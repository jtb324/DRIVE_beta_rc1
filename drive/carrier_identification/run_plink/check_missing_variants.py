# This file will contain a function that will check to see if any of the
# variants in the original variant input file were missing after the plink
# analysis file
# This function will only be used if the analysis method is gene

import pandas as pd
import os
from os import path
import glob
import re
from typing import List, Dict


def get_raw_files(plink_output_files: str) -> list:
    '''This function will get the raw files and return a list'''

    current_dir: str = os.getcwd()

    raw_file_list: list = []

    os.chdir(plink_output_files)

    for file in glob.glob("*.raw"):

        full_path = "".join([plink_output_files, file])

        raw_file_list.append(full_path)

    os.chdir(current_dir)

    return raw_file_list


def get_plink_var(raw_filepath: str) -> list:
    '''This function will load the raw file into a pandas
    dataframe and will get a list of all the variants in the dataframe'''

    # load file into pandas dataframe
    raw_df: pd.DataFrame = pd.read_csv(raw_filepath, sep=" ")

    #  getting the column titles from the dataframe
    variant_list: list = raw_df.columns.tolist()[6:]
    print(len(variant_list))
    return variant_list



def get_full_var_list(input_var_filepath: str) -> List[str]:
    '''This function will get a list of all the variants in the
    input variant file'''

    if input_var_filepath[-4:] == "xlsx":

        var_df: pd.DataFrame = pd.read_excel(input_var_filepath)

    elif input_var_filepath[-3:] == "csv":

        var_df: pd.DataFrame = pd.read_csv(input_var_filepath)

    variant_list: list = var_df["SNP"].values.tolist()

    
    return variant_list


def find_missing_variants(raw_var_list: list, input_var_list: list) -> list:
    '''This function will return a list of all variants that are in the
    input_var_list but not in the raw_var_list'''

    adjusted_raw_var_list = [variant[:-2] for variant in raw_var_list]

    missing_var_list: list = [
        variant for variant in input_var_list
        if variant not in adjusted_raw_var_list
    ]

    return missing_var_list


def write_to_file(missing_var_list: list, output: str):
    '''This function will write the missing variants to a file'''

    full_fileName: str = "".join([output, "plink_missing_variants.txt"])

    if os.path.exists(full_fileName):

        os.remove(full_fileName)

    with open(full_fileName, "w+") as missing_var_file:

        for variant in missing_var_list:

            missing_var_file.write(f"{variant}\n")


def check_for_missing_var(plink_file_path: str, input_var_filepath: str) -> int:
    '''This is the main function that will be run to generate a file
    of all variants that were missing after the initial plink run'''
    output_path: str = plink_file_path
    # gather all the raw files
    raw_file_list: list = get_raw_files(plink_file_path)
    # generating a list of variants in the raw_file

    full_raw_variant_list: list = []

    for raw_file in raw_file_list:

        raw_variant_list: list = get_plink_var(raw_file)

        full_raw_variant_list += raw_variant_list

    # generating a list of variants from the initial input file
    input_var_list: list = get_full_var_list(input_var_filepath)

    # finding the missing variants
    missing_var_list: list = find_missing_variants(full_raw_variant_list,
                                                   input_var_list)

    write_to_file(missing_var_list, output_path)

    return len(missing_var_list)
