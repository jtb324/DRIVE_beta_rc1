import logging
import re
import os
import glob
import pandas as pd

from network_creator_class import Network_Img_Maker


def get_chr_num(file: str) -> str:
    '''This function will get the chr_num from the file name'''

    match = re.search(r'chr\d\d_', file)

    # find chromosome number
    if match:

        chr_num = match.group(0)

    else:
        match = re.search(r'chr\d_', file)

        chr_num = match.group(0)

    chr_num = chr_num.strip(".").strip("_")

    print(chr_num)

    return chr_num


def create_carrier_dict(file) -> dict:
    '''This function will create a dictionary where the key is the variant 
    and the value is a list of carriers'''

    # creating the carrier dictionary
    carrier_dict: dict = dict()

    # loading the file into a dataframe
    carrier_df: pd.DataFrame = pd.read_csv(file, sep=",")

    # getting a unique list of variants
    variant_list: set = set(carrier_df["Variant ID"].tolist())

    print(variant_list)
    print(len(variant_list))

    # iterating through each variant in the list
    for variant in variant_list:

        # subsetting the dataframe for the specific variant
        carrier_df_subset: pd.DataFrame = carrier_df[carrier_df["Variant ID"] == variant]

        # Getting a list of carriers for this variant
        carrier_list: list = carrier_df_subset["IID"].values.tolist()

        carrier_dict[variant] = carrier_list

    print(f"The dictionary size is {str(carrier_dict.__sizeof__())}")

    return carrier_dict


def create_networks(segments_file_dir: str, variant_file_dir: str, ind_in_network_dict: dict, output_path: str) -> dict:

    logger = logging.getLogger(output_path+'/drawing_networks.log')

    network_drawer: object = Network_Img_Maker(
        segments_file_dir, variant_file_dir, output_path, logger)

    # Getting all of the allpair files into a list
    allpair_file_list: list = network_drawer.gather_files(
        segments_file_dir, "*.allpair.txt")

    # Getting a list of all the carrier files
    carrier_file_list: list = network_drawer.gather_files(
        variant_file_dir, "*.single_var_list_reformat.csv")

    print(allpair_file_list)
    print(carrier_file_list)

    # iterating through the carrier_file list for each file
    for file in carrier_file_list:

        # This function will get the chromosome number of the file
        chr_num: str = get_chr_num(file)

        carrier_dict: dict = create_carrier_dict(file)

        # iterating through each variant and getting the list of carriers
        for items in carrier_dict.items():

            # getting the variant from the first element of the tuple
            variant_id: str = items[0]

            carrier_list: list = items[1]

            print(chr_num)

            try:

                allpair_file_path: str = [
                    file for file in allpair_file_list if chr_num in file and variant_id in file][0]

            except IndexError:

                print(
                    f"There was no allpair.txt file found for the variant, {variant_id}")

                # TODO: add a function to write these variants to a file so that you can identify them

                continue

            # Loading the dataframe
            allpair_df: pd.DataFrame = network_drawer.load_allpair_file(
                allpair_file_path)

            # filtering the allpair_df for only rows were both individuals are carriers
            filtered_allpair_df: pd.DataFrame = network_drawer.filter_for_carriers(
                allpair_df, carrier_list)

            # This just drops any empty rows in the dataframe
            if filtered_allpair_df["pair_1"].isnull().any() or filtered_allpair_df["pair_2"].isnull().any():
                filtered_allpair_df = network_drawer.drop_empty_rows(
                    filtered_allpair_df)

            # This class method will determine the percentage of carriers in each network for each variant
            carrier_in_network_dict = network_drawer.carriers_in_network(
                carrier_list, filtered_allpair_df, ind_in_network_dict, variant_id)

            network_drawer.draw_networks(filtered_allpair_df, variant_id)

    return carrier_in_network_dict
