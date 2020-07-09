import os
import pandas as pd
from write_path import writePath
from check_directory import check_dir

# Using the ind_network_count.csv file
#####################################################################


def determine_network_sizes(network_count_filepath, pedigree_df, output_path):
    '''This function determines the size of the matched networks. The function takes a path to a network count files, a dataframe of the pedigree, and then takes an output path. This function is used within the network_sizes functions in the searchPedigree.py file.'''

    # Making a directory for the output files
    network_directory = check_dir(output_path, "network_sizes")

    fid_df = pd.read_csv(network_count_filepath, sep=",",
                         usecols=["FID", "IID"])

    fid_list = fid_df["FID"].values.tolist()

    matched_pedigree_df = pedigree_df[pedigree_df["FID"].isin(fid_list)]

    ################################################
    # Determining Network Size

    network_size_path = writePath(network_directory, "network_sizes.csv")

    network_size_df = matched_pedigree_df.groupby("FID").count().reset_index()

    network_size_df = network_size_df.rename(
        columns={"IID": "Length of Network"})

    network_size_df.to_csv(network_size_path, index=False)

    ################################################
    # Determining the distribution of network size
    network_dist_path = writePath(
        network_directory, "network_distribution_count.csv")

    network_dist_df = network_size_df.groupby(
        "Length of Network").count().reset_index()

    network_dist_df = network_dist_df.rename(
        columns={"IID": "Size of Network", "FID": "Number of Networks"})

    network_dist_df.to_csv(network_dist_path, index=False)
