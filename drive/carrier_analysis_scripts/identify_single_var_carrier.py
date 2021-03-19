# This file contains the functions used to determine how many total individuals contain some variant and how many individuals contain multiple variants

###################################################################################
# importing modules
import pandas as pd
import logging
import glob
import os
import re
import csv

###################################################################################
# importing necessary functions from other files

import file_creator_scripts
import population_filter_scripts
import utility_scripts

###################################################################################
# Function to find the total number of variants


def totalVariantIDList(iid_list: set, writeLocation: str, file_name_head: str):
    """this is a function for the inner loop that will search through each position in the row and when it encouters a one or a two it will add that to the idlist and then return so that the outer loop in the main script moves on to the next row."""

    file_name = "".join([file_name_head, ".total_variant_ID_list.txt"])

    writeDirectory = file_creator_scripts.writePath(writeLocation, file_name)

    MyFile = open(writeDirectory, "w")

    for element in iid_list:
        MyFile.write(element)
        MyFile.write("\n")
    MyFile.close()


###########################################################################################


def find_all_files(input_file_path: str):
    """This function will find a function of all files within a specified directory"""
    cur_dir = os.getcwd()

    os.chdir(input_file_path)

    recode_file_list = []

    for file in glob.glob("*.raw"):

        full_file_path = "".join([input_file_path, file])

        recode_file_list.append((full_file_path, file))

    # print(f"{len(recode_file_list)} recoded raw file found with the {os.getcwd()}")

    os.chdir(cur_dir)

    return recode_file_list


###########################################################################################
# This function determines all the individuals who have a specific variant


def get_chr_num(file_str: str) -> str:
    """Function that can identify the chromosome number in the provided file string
    Parameters
    __________
    file_str : str
        this is a file string that has a chromosome number within it

    Returns
    _______
    str
        returns a chromsome number of the form chrXX where XX are digits
    """
    match = re.search(r".chr\d\d_", file_str)

    chr_num: str = match.group(0)

    chr_num: str = chr_num.strip(".")

    file_prefix: str = chr_num.strip("_")

    return file_prefix


def run_pop_filter(pop_info: str, raw_file: str,
                   pop_code: str) -> pd.DataFrame:
    """Function to filter the raw files to individuals in a specific population
    Parameters
    __________
    pop_info : str
        file that contains information about what ancestry each grid is from

    raw_file : str
        file path to the raw file that was output from PLINK

    pop_code : str
        specified population code of interest from 1000 genomes

    Returns
    _______
    pd.DataFrame
        this is the filtered raw file loaded into a dataframe
    """
    dataset_filter = population_filter_scripts.Pop_Filter(pop_info, raw_file)

    pop_info_df, recode_df = dataset_filter.load_files()

    pop_info_subset_df = dataset_filter.get_pop_info_subset(
        pop_info_df, pop_code)

    raw_file = dataset_filter.filter_recode_df(pop_info_subset_df, recode_df)

    return raw_file


def reformatting_file(var_dict: dict, output_path: str, file_prefix: str):
    """Function to write the var dict to a reformated file
    Parameters
    __________
    var_dict : dict
        dictionary containing which single variant each grid carriers

    output_path : str
        string listing which directory to output the reformatted file to

    file_prefix : str
        string containing the chromosome number. Format should be chrXX, where
        XX are digits
    """
    var_reformat_df = pd.DataFrame(var_dict, columns=["IID", "Variant ID"])

    reformat_directory = file_creator_scripts.check_dir(
        output_path, "reformatted")

    var_reformat_df.to_csv(
        "".join([
            reformat_directory,
            "/",
            file_prefix,
            ".single_var_list_reformat.csv",
        ]),
        index=False,
    )


@utility_scripts.func_readme_generator
def single_variant_analysis(*args, **kwargs):
    """Function that identifies grids that carry at least one variant
    Parameters
    __________
    **kwargs : dict
        dictionary that contains a list of parameters. For
        this function there will be the keywords: 'recode_filepath', 'output', 'pop_info', 'pop_code'
    """
    # getting the main logger
    logger = logging.getLogger(__name__)
    # expanding parameters from the kwargs dictionary
    recodeFile: str = kwargs.get("recode_filepath")

    write_path: str = kwargs.get("output")

    pop_info: str = kwargs.get("pop_info")

    pop_code: str = kwargs.get("pop_code")

    try:

        # making a directory to put the files from this step in

        os.mkdir("".join([write_path, "carrier_analysis_output/"]))

    except FileExistsError:
        pass

    output_path: str = "".join([write_path, "carrier_analysis_output/"])

    # creating the readme
    # create_readme(output_path)

    recode_file_list = utility_scripts.get_file_list(recodeFile, "*raw")

    for file_tuple in recode_file_list:
        file_prefix: str = get_chr_num(file_tuple[1])

        output_file_name = "".join(
            [file_prefix, ".", "single_variant_carrier.csv"])

        full_output_file: str = "".join([output_path, output_file_name])

        recodeFile = file_tuple[0]

        file_checker = file_creator_scripts.Check_File_Exist(recodeFile)

        raw_file = file_checker.check_file_exist(separator=" ")

        # subsetting the raw_file for a specific population if the population code, pop_code, is provided

        if pop_code:
            raw_file: pd.DataFrame = run_pop_filter(pop_info, raw_file,
                                                    pop_code)

        column_list = raw_file.columns[6:].values.tolist()

        # var_dict = dict()

        # var_dict_reformat = dict()
        carrier_df: pd.DataFrame = pd.DataFrame()
        # total_id_set = set()
        # iid_list = []

        for column in column_list:

            IID_series: pd.Series = raw_file[raw_file[column].isin([1.0,
                                                                    2.0])].IID
            # If the IID_series has a length of 0 then there are
            # no carriers of the variant. The program will then make
            # a pandas series and insert N/A into it and then this
            # will be added to the dataframe
            if len(IID_series) == 0:
                IID_series = pd.Series("N/A")

                iid_dataframe: pd.DataFrame = IID_series.to_frame()
            else:
                iid_dataframe: pd.DataFrame = IID_series.to_frame()

            iid_dataframe["Variant ID"] = column

            iid_dataframe.columns = ["IID", "Variant ID"]

            carrier_df = pd.concat([carrier_df, iid_dataframe])

        # getting a list of IIDs
        iid_list: list = list(set(iid_dataframe.IID.values.tolist()))

        if carrier_df.empty:

            print(
                "There were no individuals found so there dictionary was not written to a csv file."
            )
            logger.info("No individuals identified as carrying variants")

        totalVariantIDList(iid_list, output_path, file_prefix)

        iid_dataframe.to_csv(full_output_file, index=False)
