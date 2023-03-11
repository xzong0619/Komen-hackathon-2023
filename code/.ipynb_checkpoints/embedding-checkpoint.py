import scanpy as sc
import pandas as pd
import argparse
import commentjson
from pathlib import Path
import scvi
from scvi.model.utils import mde

PROJECT_FOLDER = Path(__file__).parent.parent


def main():
    """
    The CLI to embed the reference data together with training data
    """
    parser = argparse.ArgumentParser(
        description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--train_h5ad",
        type=Path,
        help="The h5ad data we want to embed with reference",
        required=True,
    )
    parser.add_argument(
        "--train_parameters",
        type=Path,
        help="JSON file contains the training parameters",
        required=True,
    )
    parser.add_argument(
        "--output_folder", type=Path, help="The output folder", required=True
    )
    parser.add_argument(
        "--label_column",
        type=str,
        default="CellLine",
        help="The target label in the adata.obs",
    )
    parser.add_argument(
        "--categorical_covariate_keys",
        type=str,
        default="batch",
        help="Column names that indicate the categorical covariates you want to include in your training model. Different column names should be separated by the comma. e.g. tech,Patient",
    )
    parser.add_argument(
        "--continuous_covariate_keys",
        type=str,
        default=None,
        help="Column names that indicate the continuous covariates you want to include in your training model. Different column names should be separated by the comma. e.g. pct_counts_mt,total_umis",
    )

    args = parser.parse_args()

    # specify the output
    Output_Folder = args.output_folder / Path(args.train_parameters).name.split(".")[0]
    Path(Output_Folder).mkdir(parents=True, exist_ok=True)

    # load the train data
    adata = sc.read(args.train_h5ad)
    categorical_covariate_keys = (
        args.categorical_covariate_keys.split(",")
        if args.categorical_covariate_keys is not None
        else []
    )
    continuous_covariate_keys = (
        args.continuous_covariate_keys.split(",")
        if args.continuous_covariate_keys is not None
        else []
    )
    missing_columns = pd.Index(
        categorical_covariate_keys + categorical_covariate_keys + [args.label_column]
    ).difference(adata.obs.columns)
    # make sure each column is in the dataset
    assert (
        len(missing_columns) == 0
    ), f"cannot find {','.join(missing_columns)} in adata.obs"

    # load the train parameters
    with open(args.train_parameters, "r") as f:
        arches_params = commentjson.load(f)

    # integrate data
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
    )
    vae_ref = scvi.model.SCVI(adata, **arches_params)
    vae_ref.train()

    # transfer the cell line label into the patient
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        vae_ref, labels_key=args.label_column, unlabeled_category="Unknown"
    )
    scanvi_model.train()
    adata.obs[f"{args.label_column}_prediction"] = scanvi_model.predict()

    # visualize the target label
    adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
    adata.obsm["X_mde"] = mde(adata.obsm["X_scANVI"])
    fig = sc.pl.embedding(
        adata,
        basis="X_mde",
        color=[args.label_column, f"{args.label_column}_prediction"],
        frameon=False,
        ncols=1,
        return_fig=True
    )
    # output the embedding visualization
    fig.savefig(f"{Output_Folder}/embedding.pdf", bbox_inches="tight")

    # output the embedding feature
    X = pd.DataFrame(
        adata.obsm["X_scANVI"],
        columns=[f"Latent{i}" for i in range(adata.obsm["X_scANVI"].shape[1])],
        index=adata.obs.index,
    )
    y = adata.obs[args.label_column]

    X.to_csv(f"{Output_Folder}/train_X.csv", index=True)
    y.to_frame().to_csv(f"{Output_Folder}/train_y.csv", index=True)
    # output the prediction feature
    adata.obs[f"{args.label_column}_prediction"].to_frame().to_csv(
        f"{Output_Folder}/prediction.csv", index=True
    )


if __name__ == "__main__":
    main()
