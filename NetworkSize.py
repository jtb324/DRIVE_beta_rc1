#####################################################################
# importing modules

import os
import pandas as pd

#####################################################################
# importing necessary functions from other files

from write_path import writePath
from check_directory import check_dir

# Using the ind_network_count.csv file
#####################################################################


def network_sizes(pedigree_subset, output_path, counts_file_name, list_file_name, full_pedigree_file):
    '''This function will just output a file that gives you an idea of the size of the number of individuals matched in each network and then the distribution of the number of matched individuals in the networks.'''

    print("generating a file containing the size of each network...")

    count_directory = writePath(output_path, counts_file_name)

    network_counts = pedigree_subset.groupby(
        "FID").count()  # This is a groupby object

    network_counts.reset_index().to_csv(count_directory, index=False)

    print("generating a list of all individual carrier in each network...")

    list_directory = writePath(output_path, list_file_name)

    network_list = pedigree_subset.groupby(
        "FID")["IID"].apply(list)

    network_list.reset_index().to_csv(list_directory, index=False)

    total_network_sizes(count_directory, full_pedigree_file, output_path)

###############################################################################


def total_network_sizes(network_count_filepath, pedigree_df, output_path):
    '''This function determines the full size of the matched networks. This size includes individuals who are not carrying variants. The function takes a path to a network count files, a dataframe of the pedigree, and then takes an output path. This function is used within the network_sizes functions in the searchPedigree.py file.'''

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
