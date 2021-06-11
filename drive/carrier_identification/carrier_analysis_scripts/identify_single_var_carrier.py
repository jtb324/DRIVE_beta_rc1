# This file contains the functions used to determine how many total individuals contain some variant and how many individuals contain multiple variants

###################################################################################
# importing modules
import pandas as pd
import logging
import os
import sys
from typing import List, Dict
from tqdm import tqdm

###################################################################################
# importing necessary functions from other files

import population_filter_scripts
import utility_scripts

###########################################################################################
def write_to_file(output_filename: str, carrier_dict: Dict) -> None:
    """Function to write the carrier dictionary to a file
    Parameters
    __________
    output_filename : str
        string that contains the full filepath to write the 
        output file to
        
    carrier_dict : Dict[str, Dict[str, List[str]]]
        dictionary that has a list of carriers for each variant 
        for each chromosome 
    """

    
    with open(output_filename, "w+") as output_file:

        # writing a header line to the file
        output_file.write("iid\tvariant_id\tchr\n")

        # iterate over each chromosome, then over each variant
        # to get a list of iids that carry the variant
        for chromo, value in carrier_dict.items():
            
            for variant, carrier_list in value.items():

                for iid in carrier_list:

                    output_file.write(f"{iid}\t{variant}\t{chromo}\n")

# This function determines all the individuals who have a specific variant
def single_variant_analysis(parameter_dict: dict) -> Dict:
    """Function that identifies grids that carry at least one variant
    Parameters
    __________
    parameter_dict : dict
        dictionary that contains a list of parameters. For
        this function there will be the keywords: 'recode_filepath', 'output', 'pop_info', 'pop_code'
    
    Return
    ______
    Dict[str, Dict[str, List[str]]]
        returns a dictionary of dictionaries that contains the carriers of each variant for each chromosome
    """

    
    # getting the main logger
    logger = logging.getLogger(__name__)
    # expanding parameters from the kwargs dictionary
    recodeFile: str = parameter_dict.get("recode_filepath")

    write_path: str = parameter_dict.get("output")

    pop_info: str = parameter_dict.get("pop_info")

    pop_code: str = parameter_dict.get("pop_code")

    # checking if the output path exists and making it if it doesn't
    full_output_dir: str = utility_scripts.check_dir(
        write_path, "carrier_analysis_output/"
    )

    recode_file_list: list = utility_scripts.get_file_list(recodeFile, "*raw")

    
    # forming the full path of the output file
    full_output_file: str = os.path.join(full_output_dir, "single_variant_carriers.csv")

    # Creating the data structure that will hold the 
    # variants and the list of carriers
    carrier_dict: Dict[str, Dict[str, List[str]]] = {}

    # creating a progress bar to show how the program is iterating through the program
    for increment in tqdm(range(len(recode_file_list))):

        # getting the file string from the list based on 
        # the increment
        recodefile: str = recode_file_list[increment]
    # for recodefile in recode_file_list:

        # getting the chromosome number 
        chr_num: str = utility_scripts.get_chr_num(recodefile, r".chr\d\d_")

        # output_file_name = "".join(
        #     [
        #         chr_num,
        #         ".",
        #         "single_variant_carrier.csv",
        #     ]
        # )

        # # forming the full path of the output file
        # full_output_file: str = os.path.join(full_output_dir, output_file_name)

        # TODO: refactor these next two lines
        # load the raw_file into a dataframe
        raw_file: pd.DataFrame = pd.read_csv(recodefile, sep=" ")

        if pop_code:
            raw_file: pd.DataFrame = population_filter_scripts.run_pop_filter(
                pop_info, raw_file, pop_code
            )

        column_list = raw_file.columns[6:].values.tolist()

        carrier_df: pd.DataFrame = pd.DataFrame()

        carrier_dict.setdefault(chr_num, {})

        for column in column_list:

            IID_series: List[str] = raw_file[raw_file[column].isin([1.0, 2.0])].IID.tolist()

            # If the IID_series has a length of 0 then there are
            # no carriers of the variant. The program will then make
            # a pandas series and insert N/A into it and then this
            # will be added to the dataframe
            if len(IID_series) == 0:
                IID_series: List[str] = list(pd.Series("N/A"))

            carrier_dict[chr_num].setdefault(column, IID_series)
            #     iid_dataframe: pd.DataFrame = IID_series.to_frame()
            # else:
            #     iid_dataframe: pd.DataFrame = IID_series.to_frame()

            # iid_dataframe["Variant ID"] = column

            # iid_dataframe.columns = ["IID", "Variant ID"]

            # carrier_df = pd.concat([carrier_df, iid_dataframe])

        # counting how many total unique carriers
        # print(f"Number of unique IIDs who are identified as carrying a variant of interest is {len(list(set(iid_dataframe.IID.values.tolist())))}")

        if len(carrier_dict) == 0:

            print(
                "There were no individuals found so there dictionary was not written to a csv file."
            )
            logger.info("No individuals identified as carrying variants")

            sys.exit(0)

        # totalVariantIDList(iid_list, output_path, file_prefix)

        # carrier_df.to_csv(full_output_file, index=False)
    
    write_to_file(full_output_file, carrier_dict)

    return carrier_dict