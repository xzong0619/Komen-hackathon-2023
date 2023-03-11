import pandas as pd
import argparse
from pathlib import Path

PROJECT_FOLDER = Path(__file__).parent.parent

def main():
    """
    The CLI of train on the latent space
    """
    parser = argparse.ArgumentParser(
        description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--embedding_folder",
        type=Path,
        help="The embeddding folder that stores the single cell embedding features",
        required=True,
    )
    parser.add_argument(
        "--drug_csv",
        type=Path,
        required=True,
        help="The cleaned data response data.",
    )
    parser.add_argument(
        "--process_folder",
        type=Path,
        default=PROJECT_FOLDER / "data/processed_data",
        help="The processed data folder",
    )
    args = parser.parse_args()
    X = pd.read_csv(f"{args.embedding_folder}/train_X.csv")
    y = pd.read_csv(f"{args.embedding_folder}/train_y.csv")
    df = pd.read_csv(args.drug_csv)
    y.rename(columns={"CellLine": "CELL_LINE_NAME"}, inplace=True)
    X_y = y.merge(df, on="CELL_LINE_NAME").merge(X, on="Barcodes")

    X_pred = X.loc[y["CELL_LINE_NAME"]=='Unknown',:]
    X_y.to_csv(f"{args.process_folder}/train_X_y.csv")
    X_pred.to_csv(f"{args.process_folder}/X_pred.csv")
               

if __name__ == "__main__":
    main()