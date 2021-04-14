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

###########################################################################################
# This function determines all the individuals who have a specific variant


# @utility_scripts.func_readme_generator
def single_variant_analysis(parameter_dict: dict):
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

    # checking if the output path exists and making it if it doesn't
    full_output_dir: str =utility_scripts.check_dir(write_path, "carrier_analysis_output/")

    recode_file_list: list = utility_scripts.get_file_list(recodeFile, "*raw")

    # iterating through each file in the recode file list
    for recodefile in recode_file_list:   

        # getting the chromosome number to use as a file prefix
        # with the function get_chr_num and forming the output
        # file name
        output_file_name = "".join(
            [utility_scripts.get_chr_num(recodefile, r".chr\d\d_"), ".", "single_variant_carrier.csv"])

        # forming the full path of the output file
        full_output_file: str = os.path.join(full_output_dir, output_file_name)

        # TODO: refactor these next two lines
        # load the raw_file into a dataframe
        raw_file: pd.DataFrame = pd.read_csv(recodefile, sep=" ")

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
        # print(f"Number of unique IIDs who are identified as carrying a variant of interest is {len(list(set(iid_dataframe.IID.values.tolist())))}")

        if carrier_df.empty:

            print(
                "There were no individuals found so there dictionary was not written to a csv file."
            )
            logger.info("No individuals identified as carrying variants")

            sys.exit(0)

        # totalVariantIDList(iid_list, output_path, file_prefix)

        carrier_df.to_csv(full_output_file, index=False)
