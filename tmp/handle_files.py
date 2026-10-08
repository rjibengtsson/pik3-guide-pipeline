import os, sys
import pandas as pd


def read_csv_file(file_path):
    """
    Reads a CSV file and returns a pandas DataFrame.

    Parameters:
        file_path (str): The path to the CSV file.
    """

    df = pd.read_csv(file_path, header=0)

    # Filter the DataFrame
    df_filt = df[(df["predicted_efficacy"] >= 0.95) & (df["transcript_id"] == "NM_005026.5")]

    print(df_filt["start"].tolist())


def read_tsv_file(file_path):

    df = pd.read_csv(file_path, sep="\t", header=0)

    # Filter the DataFrame
    df_filt = df[(df["n_hits_nm"] > 1)]

    print(df_filt)



def main():
    in_file = sys.argv[1]
    # read_csv_file(in_file)

    read_tsv_file(in_file)


if __name__ == "__main__":
    main()