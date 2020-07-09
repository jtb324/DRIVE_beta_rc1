import os
import pandas as pd

# Using the ind_network_count.csv file
#####################################################################


def writePath(write_location, file_name):
    '''This section creates a function that creates a write path that can be used. This helps keep the code DRY'''

    # This line just joins the path of the directory to write to and the filename for a complete string.
    total_var_directory = os.path.join(
        write_location, file_name)

    return total_var_directory

####################################################################


def determine_network_sizes(network_count_filepath, pedigree_df, output_path):
    '''This function determines the size of the matched networks. The function takes a path to a network count files, a dataframe of the pedigree, and then takes an output path. This function is used within the network_sizes functions in the searchPedigree.py file.'''

    # Making a directory for this out put
    network_directory = writePath(output_path, "network_sizes")

    try:
        os.mkdir(network_directory)
        print("Successfully created the {} directory".format(network_directory))

    except FileExistsError:
        print("{} is already an existing directory".format(network_directory))

    fid_df = pd.read_csv(network_count_filepath, sep=",",
                         usecols=["FID", "IID"])

    fid_list = fid_df["FID"].values.tolist()

    matched_pedigree_df = pedigree_df[pedigree_df["FID"].isin(fid_list)]

    # Determining network size###################################

    network_size_path = writePath(network_directory, "network_sizes.csv")

    network_size_df = matched_pedigree_df.groupby("FID").count().reset_index()

    network_size_df = network_size_df.rename(
        columns={"IID": "Length of Network"})

    network_size_df.to_csv(network_size_path, index=False)

    # Determining the distribution of network size ############
    network_dist_path = writePath(
        network_directory, "network_distribution_count.csv")

    network_dist_df = network_size_df.groupby(
        "Length of Network").count().reset_index()

    network_dist_df = network_dist_df.rename(
        columns={"IID": "Size of Network", "FID": "Number of Networks"})

    network_dist_df.to_csv(network_dist_path, index=False)
