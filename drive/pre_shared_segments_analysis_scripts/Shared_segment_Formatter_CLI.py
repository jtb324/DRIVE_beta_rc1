#!/usr/bin/python

# THese are modules used
import re
import os
from os import path
import pandas as pd

import pre_shared_segments_analysis_scripts
import utility_scripts
import pre_shared_segments_analysis_scripts.file_dict_creator as file_dict_creator

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

    # TODO: create a function that forms the dictionary of base pairs
    return var_dict, iid_dict


# Need to check the inputs of this function
@utility_scripts.func_readme_generator
@utility_scripts.check_dir_decorator("formatted_ibd_output/")
def collect_files(*args, parameter_dict: dict) -> dict:

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
                
def iterate_file_dict(file_dict: dict, output: str):
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

            # In the var_pos_dict the keys are the variant and values 
            # are the base position of the variant and in the 
            # var_iid_dict, the keys are the variants and the values 
            # are the IID's that carry the variants
            var_pos_dict, var_iid_dict = create_var_dict(chromo_file, map_file)
            
            # need to create a function that will check if there are no variants for a 
            # specific file
            #TODO: creating a function to do the above mentioned thing
            

            # This checks if the var_info_file_path and the variant_directory are empty string because
            # this would mean that the the chromo_file only has one variant and it has no carriers
            # if variant_directory == "":
            #     print(f"There were no carriers in the file {chromo_file}")
            #     continue
            #iterating through
            # forming the dictionary to have all the variant info
            var_info_dict: dict = {}

            # iterating through the above 
            for variant, bp in var_pos_dict.items():

                var_info_dict = create_var_info_dict( var_info_dict,var_iid_dict, variant, bp)
             
            
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
