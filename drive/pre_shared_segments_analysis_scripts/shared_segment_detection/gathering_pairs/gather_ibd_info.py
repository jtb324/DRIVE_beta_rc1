#!/usr/bin/python

# THese are modules 
import re
import os
from os import path
import pandas as pd

import utility_scripts
import pre_shared_segments_analysis_scripts.file_dict_creator as file_dict_creator
from .collect_shared_segments import gather_pairs, generate_parameters, build_unique_id_dict, create_ibd_arrays, write_to_file

####################################################################################################




def create_no_carriers_file(variant_id: str, chr_num: str, output_path: str):
    """Function to write to a file if the provided variant has no carriers
    Parameters
    __________
    variant_id : str 
        string containing the variant id
    
    chr_num : str 
        string containing the chromosome number taking fro the key of the file_dict

    output_path : str 
        string listing the file path that the error file will be written to 
    """ 
    
    with open("".join([output_path, "no_carriers_in_file.txt"]), "a+") as no_carrier_file:
        no_carrier_file.write(f"{variant_id}\t")
        no_carrier_file.write(f"{chr_num[-2:]}\n")

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

    return var_dict, iid_dict


# Need to check the inputs of this function
# @utility_scripts.func_readme_generator
def collect_files(parameter_dict: dict) -> dict:
    
    
    # checking the kwargs dictionary for specific
    # parameters
    ibd_files: str = parameter_dict.get("ibd_file_path")
    carrier_file: str = parameter_dict.get("carrier_file")
    ibd_program: str = parameter_dict.get("ibd_program")
    output: str = parameter_dict.get("output")
    map_file_dir: str = parameter_dict.get("map_files")
    file_suffix: str = parameter_dict.get("ibd_file_suffix")
    min_CM: str = parameter_dict.get("min_CM_threshold")
    threads: str = parameter_dict.get("threads")
    ###########################################################

    utility_scripts.check_dir(output, "formatted_ibd_output/")
    # getting the segment file list
    segment_file_list: list = utility_scripts.get_file_list(ibd_files, file_suffix)

    chr_var_file_list: list = utility_scripts.get_file_list(carrier_file, "*.single_variant_carrier.csv")

    # getting the list of map files
    map_file_list: list = utility_scripts.get_file_list(map_file_dir, "*.map")

    # creating a dictionary that contains the correct map, carrier, and ibd files for the correct chromosome
    file_dict: dict = file_dict_creator.get_file_dict(map_file_list, chr_var_file_list, segment_file_list)
    
    # filter the file dict to only the keys with values for each chromosome key
    file_dict = file_dict_creator.filter_empty_dictionaries(file_dict)

    return file_dict

def create_var_info_dict(var_info_dict: dict, var_iid_dict: dict, variant: str, bp: str) -> dict:
    """Function that will add the variants and the corresponding info about iids that carry and the bp to the var_info_dict
    Parameters
    __________
    var_info_dict : dict
        empty dictionary that will have keys of variant ids and 
        values are dictionaries containing the list of iids and the
        base position
    
    var_iid_dict : dict
        dictionary where keys are the variant id and the values are 
        the list of iids that carry based on mega array

    variant : str
        string containing the variant id

    bp : str
        this is a string for the base pair position of the variant

    Returns
    _______
    dict
        returns the filled in var_info_dict
    """

    iid_list: list = var_iid_dict[variant]

    var_info_dict[variant] = {"base_pos": bp, "iid_list": iid_list}

    return var_info_dict

# unit test this function
def filter_no_carriers(var_iid_dict: dict, output_path: str, chr_num: str) -> dict:
    """Function to filter all variants that have no carriers out 
    of the var_iid_dict
    Parameters
    __________
    var_iid_dict : dict
        dictionary containing variant ids as the keys and list
        of iids who carry the variant as values
    
    output_path : str
        string listing the directory to output the file to

    chr_num : str
        string that tells which chromosome the variant should be 
        on. This will be of the format chrXX where X is a digit
    """  
    # creating a dictionary of filtered variants to return
    filtered_var_dict: dict = {}

    # iterating through all the variants to see if there is a 
    # variant with no carriers     
    for variant in var_iid_dict:
        # if the variant has no carriers write this to a file
        if "None" in var_iid_dict[variant]:

            create_no_carriers_file(variant, chr_num, output_path)

        else:
            filtered_var_dict[variant] = var_iid_dict[variant]

    return filtered_var_dict

#TODO: refactor to make this function testable
# This function is not testable at the moment
def iterate_file_dict(file_dict: dict, output: str, threads: str, ibd_program: str, min_CM: str):
    """Function will iterate through the file dictionary which has paired the chromosome number with the appropriate files 
    Parameters
    __________
    file_dict : dict
        dictionary of dictionaries where outer dictionaries keys 
        are the chromosome number and the values are dictionaries
        where the inner key is the file type and the value is the 
        filepath 
    
    output : str
        string that list directory to output files at
    """
    # Iterating through the chromosomes that have a value
    for key in file_dict:

        if len(file_dict[key]) != 1:

            chromo_file: str = file_dict[key]["carrier"]
            map_file: str = file_dict[key]["map"]
            ibd_file: str = file_dict[key]["ibd"]

            # In the var_pos_dict the keys are the variant and values 
            # are the base position of the variant and in the 
            # var_iid_dict, the keys are the variants and the values 
            # are the IID's that carry the variants
            var_pos_dict, var_iid_dict = create_var_dict(chromo_file, map_file)
            
            # need to create a function that will check if there are no carriers for a 
            # specific file. If there are no carriers it will write that to a file. It
            # returns a dictionary that only contains variant that have carriers
            filter_no_carriers(var_iid_dict, output, key)

            #iterating through
            # forming the dictionary to have all the variant info
            var_info_dict: dict = {}

            # iterating through the above 
            for variant, bp in var_pos_dict.items():

                var_info_dict = create_var_info_dict( var_info_dict,var_iid_dict, variant, bp)

            # need to fix this part for the new function
            parallel_runner: object = utility_scripts.Segment_Parallel_Runner(
                int(threads), output, ibd_program, min_CM, var_info_dict,
                ibd_file)

            error_filename: str = "nopairs-identified.txt"
            header_str: str = "variant_id"

            parallel_runner.run_segments_parallel(error_filename, gather_shared_segments,
                                                  header_str)

# TODO: rename function
def gather_shared_segments(segment_file: str, output_path: str, ibd_format: list,
             min_CM: str, var_info_dict: list, que_object, variant):

    variant_position: int = int(var_info_dict[variant]["base_pos"])

    print(f"running the variant {variant}")

    carrier_list: list = var_info_dict[variant]["iid_list"]

    parameter_dict: dict = generate_parameters(ibd_format)

    uniqID: dict = build_unique_id_dict(carrier_list)

    IBDdata, IBDindex = create_ibd_arrays()

    chr_num: str = gather_pairs(IBDdata, IBDdata, parameter_dict, segment_file, uniqID, min_CM, variant_position) 

    write_to_file(que_object)
