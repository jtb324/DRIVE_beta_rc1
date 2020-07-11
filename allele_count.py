# Determining the allele count in families
import pandas as pd
import numpy as np

################################################


def allele_counts(input_path, output_path):

    raw_file = pd.read_csv(input_path[0], sep=" ")

    pedigree_list = pd.read_csv(input_path[1], sep=",")

    for tuple in pedigree_list.itertuples(index=False):

        iid_list = tuple[1].strip("[]").replace(
            "'", "").replace(" ", "").split(",")

        variant = tuple[0].strip("[]").replace("(", "").replace(
            "'", "").replace(" ", "").split(",")[0]

        network = tuple[0].strip("[]").replace(")", "").replace(
            "'", "").replace(" ", "").split(",")[1]

        print(variant)

        print(network)
