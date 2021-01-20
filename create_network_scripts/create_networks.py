import logging
import re
import os
import pandas as pd

import create_network_scripts


def get_chr_num(file: str) -> str:
    '''This function will get the chr_num from the file name'''

    match = re.search(r'chr\d\d.', file)

    # find chromosome number
    chr_num: str = match.group(0)
    chr_num = chr_num.strip(".")

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

    # iterating through each variant in the list
    for variant in variant_list:

        # subsetting the dataframe for the specific variant
        carrier_df_subset: pd.DataFrame = carrier_df[carrier_df["Variant ID"]
                                                     == variant]

        # Getting a list of carriers for this variant
        carrier_list: list = carrier_df_subset["IID"].values.tolist()

        carrier_dict[variant] = carrier_list


    return carrier_dict


def check_file_exist(file_path: str):
    '''This function will check if the provided file exist and then 
    delete it if it does'''

    # checking if the file exist
    if os.path.isfile(file_path):

        # removing the file if it exist
        os.remove(file_path)


def add_header_row(output_path: str):
    '''This function takes the output file path and then adds a header row to it'''

    # Loading the data into a dataframe
    output_df: pd.DataFrame = pd.read_csv(output_path)

    # writing the file to the same csv with a new header
    output_df.to_csv(output_path,
                     header=[
                         "grid", "in_network_status", "network_id",
                         " variant_id", "chr_num"
                     ])


def get_size(network_groups_file: str):
    '''This function will group the provided file so that size of each network per variant
    can be determined'''

    # Load the file into a dataframe
    network_groups_df: pd.DataFrame = pd.read_csv(network_groups_file)

    # group the dataframe by two columns and then count the values in the group
    # this block catches an error where for some reason there is a space in teh variant id column name
    try:
        grouped_df: pd.DataFrame = network_groups_df.groupby(
            ["Network ID", "variant_id"]).sum()

    except KeyError:
        # This removes the space and then tries to group it again
        network_groups_df = network_groups_df.rename(
            columns={" variant_id": "variant_id"})

        grouped_df: pd.DataFrame = network_groups_df.groupby(
            ["Network ID", "variant_id"]).sum()
    # dropping the chromosome column because it is unnecessary
    final_df: pd.DataFrame = grouped_df.drop("chr_num", 1)

    final_df.to_csv("network_groups_summary.csv", index=None)


def write_missing_variants(variant: str, output_path: str):
    '''This function will create an output file that documents which variants there was no allpair.txt file for'''

    output_file_path: str = "".join([output_path, "missing_allpairs.txt"])

    # opening the output_file
    with open(output_file_path, "a+") as missing_files:

        # checking to see if the files are empty
        if os.path.getsize(missing_files) == 0:

            missing_files.write(
                "Showing a list of variants for which an allpairs.txt file was not found"
            )
        else:
            missing_files.write(variant)
            missing_files.write("\n")


def create_networks(segments_file_dir: str, variant_file_dir: str,
                    ind_in_network_dict: dict, output_path: str) -> dict:
    check_file_exist("".join([output_path, "/network_groups", ".csv"]))

    check_file_exist("".join([output_path, "/pairs_in_networks", ".csv"]))

    # checking if the file that documnets whether or not there is an allpair file exist
    check_file_exist("".join([output_path, "missing_allpairs.txt"]))

    # Creating an empty dataframe that the pairs will be written to at the end in the
    # draw networks function
    pairs_df: pd.DataFrame = pd.DataFrame()

    network_drawer: object = create_network_scripts.Network_Img_Maker(
        segments_file_dir, variant_file_dir, output_path)
    # Getting all of the allpair files into a list
    allpair_file_list: list = network_drawer.gather_files(
        segments_file_dir, "*.allpair.txt")
    # Getting a list of all the carrier files
    carrier_file_list: list = network_drawer.gather_files(
        "".join([variant_file_dir, "reformated/"]),
        "*.single_var_list_reformat.csv")

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

            try:

                allpair_file_path: str = [
                    file for file in allpair_file_list
                    if chr_num in file and variant_id in file
                ][0]

            except IndexError:

                print(
                    f"There was no allpair.txt file found for the variant, {variant_id}"
                )

                # TODO: add a function to write these variants to a file so that you can identify them

                continue
            # Loading the dataframe
            allpair_df: pd.DataFrame = network_drawer.load_allpair_file(
                allpair_file_path)

            # filtering the allpair_df for only rows were both individuals are carriers
            filtered_allpair_df: pd.DataFrame = network_drawer.filter_for_carriers(
                allpair_df, carrier_list)

            # if the dataframe is empty then the loop needs to skip that variant
            if filtered_allpair_df.empty:
                print(
                    f"There were no pairs found where both individuals were carriers for the variant {variant_id}"
                )

                # prints the message and then moves onto the next iteration
                continue

            # This just drops any empty rows in the dataframe
            if filtered_allpair_df["pair_1"].isnull().any(
            ) or filtered_allpair_df["pair_2"].isnull().any():
                filtered_allpair_df = network_drawer.drop_empty_rows(
                    filtered_allpair_df)

            # This makes a dataframe that will list the pairs
            # network_drawer.making_pairs_df(filtered_allpair_df)

            # This class method will determine the percentage of carriers in each network for each variant
            carrier_in_network_dict = network_drawer.carriers_in_network(
                carrier_list, filtered_allpair_df, ind_in_network_dict,
                variant_id)

            output_path, pairs_df = network_drawer.draw_networks(
                filtered_allpair_df, variant_id, chr_num, pairs_df)

            # add_header_row(output_path)

            get_size(output_path)

    return carrier_in_network_dict
