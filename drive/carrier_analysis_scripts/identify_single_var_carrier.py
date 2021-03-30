# This file contains the functions used to determine how many total individuals contain some variant and how many individuals contain multiple variants

###################################################################################
# importing modules
import pandas as pd
import logging
import glob
import os
import re
import sys

###################################################################################
# importing necessary functions from other files

import population_filter_scripts
import utility_scripts

###################################################################################
# Function to find the total number of variants


# def totalVariantIDList(iid_list: set, writeLocation: str, file_name_head: str):
#     """this is a function for the inner loop that will search through each position in the row and when it encouters a one or a two it will add that to the idlist and then return so that the outer loop in the main script moves on to the next row."""

#     file_name = "".join([file_name_head, ".total_variant_ID_list.txt"])

#     writeDirectory = file_creator_scripts.writePath(writeLocation, file_name)

#     MyFile = open(writeDirectory, "w")

#     for element in iid_list:
#         MyFile.write(element)
#         MyFile.write("\n")
#     MyFile.close()


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


# create a decorator that can be used to check if a 
# directory exist

@utility_scripts.func_readme_generator
@utility_scripts.check_dir_decorator("carrier_analysis_output/")
def single_variant_analysis(*args, parameter_dict: dict):
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
    recodeFile: str = parameter_dict.get("recode_filepath")

    write_path: str = parameter_dict.get("output")
    
    pop_info: str = parameter_dict.get("pop_info")

    pop_code: str = parameter_dict.get("pop_code")

    output_path: str = "".join([write_path, "carrier_analysis_output/"])

    recode_file_list: list = utility_scripts.get_file_list(recodeFile, "*raw")

    # iterating through each file in the recode file list
    for recodefile in recode_file_list:   
        # getting the chromosome number to use as a file prefix
        # with the function get_chr_num and forming the output
        # file name
        output_file_name = "".join(
            [get_chr_num(recodefile), ".", "single_variant_carrier.csv"])

        # forming the full path of the output file
        full_output_file: str = os.path.join(output_path, output_file_name)

        # TODO: refactor these next two lines
        # load the raw_file into a dataframe
        raw_file: pd.DataFrame = pd.read_csv(recodefile, sep=" ")

        # file_checker = file_creator_scripts.Check_File_Exist(recodeFile)

        # raw_file = file_checker.check_file_exist(separator=" ")

        # subsetting the raw_file for a specific population if the population code, pop_code, is provided

        if pop_code:
            raw_file: pd.DataFrame = population_filter_scripts.run_pop_filter(pop_info, raw_file,
                                                    pop_code)

        column_list = raw_file.columns[6:].values.tolist()

    
        carrier_df: pd.DataFrame = pd.DataFrame()


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

        # counting how many total unique carriers
        print(f"Number of unique IIDs who are identified as carrying a variant of interest is {len(list(set(iid_dataframe.IID.values.tolist())))}")

        if carrier_df.empty:

            print(
                "There were no individuals found so there dictionary was not written to a csv file."
            )
            logger.info("No individuals identified as carrying variants")

            sys.exit(0)

        # totalVariantIDList(iid_list, output_path, file_prefix)

        iid_dataframe.to_csv(full_output_file, index=False)
