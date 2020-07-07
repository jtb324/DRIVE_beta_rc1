import os
import pandas as pd

# Using the ind_network_count.csv file


def determine_network_sizes(input_path, output_path, file_name):

    fid_df = pd.read_csv(input_path[0], sep=",",
                         usecols=["FID", "IID"])

    fid_list = fid_df["FID"].values.tolist()

    for fid in fid_list:

        network_directory = os.path.join(input_path[1], fid)

        for files in os.listdir(network_directory):

            if files.endswith(".genome_networks"):

                network_file = os.path.join(network_directory, files)

                with open(network_file, "r") as network:

                    next(network)

                    for row in network:
                        print(row)
