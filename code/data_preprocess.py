import scanpy as sc
import pandas as pd
from typing import List
from anndata import AnnData
import argparse
from pathlib import Path
import logging

logger = logging.getLogger(__name__)
PROJECT_FOLDER = Path(__file__).parent.parent

def clean_drug_data(raw_drug_data: Path) -> pd.DataFrame:
    """
    Clean the Drug response data for the final prediction
    1. Select the Breast Tissue
    2. Summarize the duplicates by the median
    Parameters
    ----------
    raw_drug_data : Path
        path where the drug response data is stored. The file format should be CSV format.
        Following columns are required in the file:
        ``['ANCHOR_NAME','LIBRARY_NAME','ANCHOR_CONC','LIBRARY_CONC','CELL_LINE_NAME','SYNERGY_OBS_EMAX','Tissue']``

    Returns
    -------
    pd.DataFrame
        cleaned drug response data with following columns:
        ``['ANCHOR_NAME','LIBRARY_NAME','ANCHOR_CONC','LIBRARY_CONC','CELL_LINE_NAME','SYNERGY_OBS_EMAX']``
    """
    df = pd.read_csv(raw_drug_data)
    df = df.loc[df['Tissue']=='Breast',:]
    # there are duplicates of each cell line on each drug combination
    # we take the median of them
    df_d_dup = (
        df.groupby(
            [
                "ANCHOR_NAME",
                "LIBRARY_NAME",
                "ANCHOR_CONC",
                "LIBRARY_CONC",
                "CELL_LINE_NAME",
            ]
        )["SYNERGY_OBS_EMAX"]
        .median()
        .rename("SYNERGY_OBS_EMAX")
        .reset_index()
    )
    df_d_dup["CELL_LINE_NAME"] = df_d_dup["CELL_LINE_NAME"].str.upper()
    return df_d_dup


def clean_sc_data(
    adata_ref: AnnData,
    adata_query: AnnData,
    celltype_column: str,
    celltype_list: List[str],
) -> AnnData:
    """
    clean the single cell data from reference and query dataset
    - take the overlap genes between reference and query data as the input features
    - filter out low quality cells
    - attach a new column `batch` indicate ref and test
    - select 3000 highly variable genes

    Parameters
    ----------
    adata_ref : AnnData
        single cell reference data (cell x gene)
    adata_query : AnnData
        single cell query data (cell x gene)
    celltype_column : str
        the column indicate cell type in `adata_query.obs`
    celltype_list : List[str]
        the categories indicate malignant cells in the `adata_query.obs[<celltype_column>]`

    Returns
    -------
    AnnData
        concatenated reference and test data with 3000 highly variable genes. (cell x 3000)
    """
    # subtract the malignant cells
    assert (
        celltype_column in adata_query.obs
    ), f"Cannot find {celltype_column} in the query data (adata_query)."
    adata_query = adata_query[adata_query.obs[celltype_column].isin(celltype_list), :]
    adata_query.var_names_make_unique()
    adata_ref.var_names_make_unique()
    # overlap genes between reference adata and query data
    overlap_genes = adata_ref.var.index.intersection(adata_query.var.index).unique()
    adata_ref = adata_ref[:, overlap_genes]
    adata_query = adata_query[:, overlap_genes]

    # concatenate dataset
    adata_ref.obs["batch"] = "ref"
    adata_query.obs["batch"] = "test"
    adata = adata_ref.concatenate(adata_query)

    # preprocess the single cell data
    ## filter out cells with less then 200 gene been detected
    sc.pp.filter_cells(adata, min_genes=200)
    ## filter out genes that are expressed less than 3 cells
    sc.pp.filter_genes(adata, min_cells=3)
    ## calculate the mitochondrial reads pear cells
    adata.var["mt"] = adata.var_names.str.startswith(
        "MT-"
    )  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    ## filter out cells with high number of detected genes
    adata = adata[adata.obs.n_genes_by_counts < 6000, :]
    ## filter out cells with high percentage of mitochondrial reads
    adata = adata[adata.obs.pct_counts_mt < 15, :]
    ## store the raw count laryer for the training purpose
    adata.layers["counts"] = adata.X.copy()
    ## normalize the cell library size so that the cell expression data are comparable
    sc.pp.normalize_total(adata, target_sum=1e4)
    ## select the high variable genes
    sc.pp.highly_variable_genes(
        adata, n_top_genes=3000, subset=True, flavor="seurat_v3", batch_key="batch"
    )
    adata.obs.index.name = "Barcodes"
    adata.var.index.name = "Gene"
    for c in adata.obs:
        if adata.obs[c].dtype != "float" and adata.obs[c].dtype != "int":
            adata.obs[c] = (
                adata.obs[c].astype(str).replace("nan", "Unknown").astype("category")
            )
    return adata


def main():
    """
    The CLI of data preprocess
    """
    parser = argparse.ArgumentParser(
        description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--test_h5ad",
        type=Path,
        help="The h5ad data we want to embed with reference",
        required=True,
    )
    parser.add_argument(
        "--drug_csv",
        type=Path,
        help="The drug response data",
        required=True,
    )
    parser.add_argument(
        "--celltype_column",
        type=str,
        default="Paper_CellType",
        help="The column in your h5ad data indicates your celltype column",
    )
    parser.add_argument(
        "--celltype",
        type=str,
        default="Cancer_cell",
        help="The categories in the cell type column indicate the malignant cells. Different category should be separated by the comma. e.g. cancer_1,cancer_2",
    )
    parser.add_argument(
        "--ref_folder",
        type=Path,
        default=PROJECT_FOLDER
        / "data/raw_data/Single_Cell_Breast_Cancer_cell-line_Atlas",
        help="The reference data folder",
    )
    parser.add_argument(
        "--process_folder",
        type=Path,
        default=PROJECT_FOLDER / "data/processed_data",
        help="The processed data folder",
    )

    args = parser.parse_args()
    Reference_Folder = args.ref_folder
    Processed_Folder = args.process_folder
    logger.info("Start cleaning the drug data.")
    drug = clean_drug_data(raw_drug_data=args.drug_csv)
    drug.to_csv(f"{Processed_Folder}/drug_response.csv")
    logger.info(
        f"Finished the cleaning and store the processed data to {Processed_Folder}/drug_response.csv"
    )

    logger.info("Start loading single cell data.")
    # load the patient single cell data
    adata_query = sc.read(args.test_h5ad)
    # load the cell line reference data
    adata_ref = sc.read_mtx(f"{Reference_Folder}/matrix.mtx.gz").T
    adata_ref.var.index = pd.read_csv(
        f"{Reference_Folder}/features.tsv.gz", sep="\t", header=None, index_col=0
    ).index
    adata_ref.obs.index = pd.read_csv(
        f"{Reference_Folder}/barcodes.tsv.gz", sep="\t", header=None, index_col=0
    ).index

    adata_ref.obs["CellLine"] = adata_ref.obs.index.map(lambda x: x.split("_")[0])
    # rename the gene annotation by the hugo symbol
    annot = (
        sc.queries.biomart_annotations(
            "hsapiens",
            ["ensembl_gene_id", "hgnc_symbol"],
        )
        .drop_duplicates("ensembl_gene_id")
        .set_index("ensembl_gene_id")
    )
    adata_ref.var["hgnc_symbol"] = annot.loc[
        adata_ref.var.index.intersection(annot.index), :
    ]["hgnc_symbol"]

    # remove genes without the hugo symbol
    adata_ref = adata_ref[:, ~adata_ref.var.hgnc_symbol.isna()]
    adata_ref.var.index = adata_ref.var.hgnc_symbol
    adata_ref.var_names = adata_ref.var.hgnc_symbol
    adata_ref.var_names_make_unique()

    logger.info("Start cleaning single cell data.")
    adata = clean_sc_data(
        adata_ref=adata_ref,
        adata_query=adata_query,
        celltype_column=args.celltype_column,
        celltype_list=args.celltype.split(","),
    )

    adata.write(f"{Processed_Folder}/train.h5ad", compression="gzip")
    logger.info(
        f"Finished the cleaning and store the processed data to {Processed_Folder}/train.h5ad"
    )

if __name__ == "__main__":
    main()
