#!/usr/bin/python

# THese are modules used
import re
import os
from os import path
import pandas as pd

import pre_shared_segments_analysis_scripts
import utility_scripts

####################################################################################################
def form_chr_num(number: int) -> str:
    """function to form the chromosome number for the key value in the get_file_dict dictionary
    Parameters
    __________
    number : int
        an integer number from the range function 
        
    Returns
    _______
    str
        returns a string of the form chrXX where X is a number
    """
    if number < 10:
        return "chr0"+str(number)
    else:
        return "chr"+str(number)

def find_match(file_list: list, chr_num: str, file_type: str) -> str:
    """Function to find the file that matches the chromosome number
    Parameters
    __________
    file_list : list 
        list of files from the get_file_dict function

    chr_num : str
        chromosome number that will be used to find the match in the 
        above files   

    file_type : str
        string that indicates if it is the map file, the carrier file, or the iid file
    Returns
    _______
    str
        returns the filepath of the matched file
    """

    # if the chromosome number is less than 10 you have to remove the
    # zero for the ibd files
    if int(chr_num[-2:]) < 10:
        alt_chr_num: str = chr_num.replace("0", "")
    else: 
        alt_chr_num = chr_num
    
    # using a dictionary to do the appropriate conditioning during the 
    # list comprehension for each file type
    file_handler = {
            "map": [file for file in file_list if "".join([chr_num, "_"]) in file],
            "carrier":[file for file in file_list if "".join([chr_num, "."]) in file],
            "ibd": [file for file in file_list if "".join([alt_chr_num, "."]) in file]
            }

    try:
        
        file: str = file_handler[file_type][0]

    except IndexError:
        return "None"
    
    return file

def get_file_dict(map_file_list: list, carrier_file_list: list, ibd_file_list: list) -> dict:
    """Function to match up all of the files for the correct chromosome into a dictionary
    Parameters
    __________
    map_file_list : list
        a list of all the map files per chromosome as a result of running Plink

    carrier_file_list : list
        a list pf all the files that have the identified carriers per chromosome
        
    ibd_file_list : list
        a list of all the files output by the specified ibd programs
        
    Returns
    _______
    dictionary
        returns a dictionary where the keys are the chromsomes and the values are dictionaries with the corresponding files
    """
    # forming a dictionary of dictionaries where the key is the 
    # chromosome number
    file_dict: dict = {form_chr_num(i):{} for i in range(1,22)}

    # get the files from each list
    for key in file_dict:

        file_dict[key].setdefault("carrier", find_match(carrier_file_list, key, "carrier"))
        file_dict[key].setdefault("map", find_match(map_file_list, key, "map"))
        file_dict[key].setdefault("ibd", find_match(ibd_file_list, key, "ibd"))
        
    return file_dict
    
def create_no_carriers_file(iid_dict: dict):
    """Function to write to a file if the provided iid has no""" 
    pass

def create_iid_dict(variant: str, iid_dict: dict, carrier_df: pd.DataFrame) -> dict:
    """Function to create a dictionary containing a list of IIDs that are identified as carrying the variants
    Parameters
    __________
    variant : str
        string containing the variant id
    
    iid_dict : dict
        empty dictionary where the keys will be the variants and the values will be 
        the list of carriers that were identified as carriers on the MEGA array
    
    carrier_df : pd.DataFrame
        dataframe formed from the single_var_carriers.csv file
    
    Returns
    _______
    dict
        dictionary that contains the variants as a key and the values are the list of 
        carriers for that variant
    """
    iid_list: list = carrier_df[carrier_df["Variant ID"] == variant]["IID"].values.tolist()

    iid_dict[variant] = iid_list

    return iid_dict

def create_dict_with_var_pos(variant: str,  var_pos_dict: dict, map_file_df: pd.DataFrame) -> dict:
    """Function to create a dictionary where the keys are the variants and the values are the base positions
    Parameters
    __________
    variant : str
        string containing the variant id

    var_pos_dict : dict
        empty dictionary where the keys will be the variant and the 
        values will be the base positions  

    map_file_df : pd.DataFrame
        dataframe of the map file that was output by plink

    Returns 
    _______
    dict
        returns the above dictionary but fillled in
    """
    var_pos_dict[variant] = map_file_df[map_file_df["variant id"] == variant]["site"].values.tolist()[0]

    return var_pos_dict

#TODO: Need to write the unit test for this section
def create_var_dict(chromo_var_file: str, map_file: str) -> tuple:
    """function that will create a dictionary of the variants and all of 
    the corresponding information 
    Parameters
    __________
    chromo_var_file : str
        file path to the single_var_carriers.csv file 
    
    map_file : str
        filepath to the map file that was output by plink
        
    Returns
    _______
    tuple
        returns a tuple dictionaries where in the first the keys are the variant and values are the base position of the variant and in the second, the keys are the variants and the values are the IID's that carry the variants
    """
    # Return a dictionary of IIDs that carry the variant
    iid_var_dict: dict = {}
    var_pos_dict: dict = {}

    # load in the csv file that list the IIDs of grids per variant on a 
    # specific chromosome
    carrier_df = pd.read_csv(chromo_var_file, sep=",")

    # reading in the map file
    map_file_df = pd.read_csv(map_file,
                                sep="\t",
                                header=None,
                                names=["chr", "variant id", "cM", "site"])
    
    # get the specific list of variants
    carrier_df_var_list: list = list(set(carrier_df["Variant ID"].values.tolist()))

    # These next two lines create a dictionary that has all the variants as keys and then all the iids that carry that variant as values
    iid_dict_list: list = map(create_iid_dict, carrier_df_var_list, iid_var_dict, carrier_df)

    iid_dict: dict = {key:value for element in iid_dict_list for key, value in element.items()}

    var_pos_dict_list: list = map(create_dict_with_var_pos, carrier_df_var_list, var_pos_dict, map_file_df)

    var_dict: dict = {key:value for element in var_pos_dict_list for key,value in element.items()}

    # TODO: create a function that forms the dictionary of base pairs
    return var_dict, iid_dict


def filter_empty_dictionaries(file_dict: dict) -> dict:
    """Function to filter the file dictionary to only those keys that contain a value
    Parameters
    __________
    file_dict : dict
        dictionary contain keys for each chromosome and dictionaries 
        with filepaths as values
    
    Returns
    _______
    dict
        returns the same dictionary as inputted but only the keys with 
        values not empty dictionaries
    """
    return {key:file_dict[key] for key in file_dict if len(file_dict[key]) != 0}
     
# Initializing the parser
@utility_scripts.func_readme_generator
@utility_scripts.check_dir_decorator("formatted_ibd_output/")
def convert_ibd(*args, parameter_dict: dict):

    # pulling the dictionary out of the args tuple
    # need to refactor all the inputs here
    function_param_dict: dict = args[0]

    # checking the kwargs dictionary for specific
    # parameters
    ibd_files: str = function_param_dict.get("ibd_file_path")
    carrier_file: str = function_param_dict.get("carrier_file")
    ibd_program: str = function_param_dict.get("ibd_program")
    output: str = function_param_dict.get("output")
    map_file_dir: str = function_param_dict.get("map_files")
    file_suffix: str = function_param_dict.get("ibd_file_suffix")
    min_CM: str = function_param_dict.get("min_CM_threshold")
    threads: str = function_param_dict.get("threads")
    ###########################################################
    # This first section will be used to get the shared segment files for each chromosome
    # getting rid of the class
    # creating a directory
    
    # getting the segment file list
    segment_file_list: list = utility_scripts.get_file_list(ibd_files, file_suffix)


    # also need to get all of the possible variant files per chromosome
    chr_var_file_list: list = utility_scripts.get_file_list(carrier_file, "*.single_variant_carrier.csv")

    # getting the list of map files
    map_file_list: list = utility_scripts.get_file_list(map_file_dir, "*.map")

    # creating a dictionary that contains the correct map, carrier, and ibd files for the correct chromosome
    file_dict: dict = get_file_dict(map_file_list, chr_var_file_list, segment_file_list)
    
    # filter the file dict to only the keys with values for each chromosome key
    file_dict = filter_empty_dictionaries(file_dict)

    # Iterating through the chromosomes that have a value
    for key in file_dict:

        if len(file_dict[key]) != 1:

            chromo_file: str = file_dict[key]["carrier"]
            map_file: str = file_dict[key]["map"]
            
            create_var_dict()
            #TODO: This is where we are in refactoring the code
            var_info_df, variant_directory = preformater.create_variant_lists(
                chromo_file, map_file)

            # This checks if the var_info_file_path and the variant_directory are empty string because
            # this would mean that the the chromo_file only has one variant and it has no carriers
            if variant_directory == "":
                print(f"There were no carriers in the file {chromo_file}")
                continue

            iid_file_list: list = utility_scripts.get_file_list(variant_directory, "*.txt")

            variant_bp_list = var_info_df.site.values.tolist()

            variant_id_list = var_info_df.variant_id.values.tolist()

            # creating a list the base pairs with the variant id
            var_info_list: list = [
                (var_bp, var_id)
                for var_bp, var_id in zip(variant_bp_list, variant_id_list)
            ]

            # creating a dictionary to couple the iid_file_list and the var_info_list
            # into a single data structure
            file_list_dict: dict = {
                "iid_files": iid_file_list,
                "var_info_files": var_info_list
            }
            parallel_runner: object = utility_scripts.Segment_Parallel_Runner(
                int(threads), output, ibd_program, min_CM, file_list_dict,
                segment_file)

            error_filename: str = "nopairs-identified.txt"
            header_str: str = "variant_id"

            parallel_runner.run_segments_parallel(error_filename, run_main,
                                                  header_str)


def run_main(segment_file: str, output_path: str, ibd_format: list,
             min_CM: str, iid_file_list: list, que_object, var_info_tuple):

    variant_position = int(var_info_tuple[0])

    variant_id = str(var_info_tuple[1])
    print(f"running the variant {variant_id}")

    iid_file = [
        iid_file for iid_file in iid_file_list if variant_id in iid_file
    ][0]

    pheno_file = iid_file

    ibd_file_converter = pre_shared_segments_analysis_scripts.Shared_Segment_Convert(
        segment_file, pheno_file, output_path, ibd_format, min_CM, 1,
        variant_position, variant_id)

    parameter_dict = ibd_file_converter.generate_parameters()

    uniqID, dupID = ibd_file_converter.build_id_pairs()

    IBDdata, IBDindex = ibd_file_converter.create_ibd_arrays()

    ibd_file_converter.run(IBDdata, IBDindex, parameter_dict, uniqID,
                           que_object)
