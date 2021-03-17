#!/usr/bin/python

# THese are modules used
import re
import os
from os import path
import multiprocessing as mp
from functools import partial

import pre_shared_segments_analysis_scripts
import utility_scripts

####################################################################################################

# Initializing the parser


def remove_previous_file(file_path: str):
    """Function to remove the file from a previous run
    Parameters
    __________
    file_path : str
        string listing the filepath to the output from previously running this program
    """
    # This section will check if the output file exist from a previous run and if it does then it will delete it
    if path.exists(file_path):

        os.remove(file_path)


def alternate_chr_num_format(chr_num: str) -> str:
    '''This removes the string in some of the chromosome files'''

    zero_indx: int = chr_num.find("0")

    if zero_indx and zero_indx == 3:

        chr_list: list = chr_num.split("0")

        chr_num = "".join(chr_list)

    return chr_num


@utility_scripts.func_readme_generator
def convert_ibd(*args, **kwargs):

    # pulling the dictionary out of the args tuple
    function_param_dict: dict = args[0]

    # checking the kwargs dictionary for specific
    # parameters
    ibd_files: str = function_param_dict.get("ibd_file_path")
    carrier_file: str = function_param_dict.get("carrier_file")
    ibd_program: str = function_param_dict.get("ibd_program")
    output: str = function_param_dict.get("output_path")
    map_file_dir: str = function_param_dict.get("map_files")
    file_suffix: str = function_param_dict.get("ibd_file_suffix")
    min_CM: str = function_param_dict.get("min_CM_threshold")
    threads: str = function_param_dict.get("threads")
    ###########################################################
    # This first section will be used to get the shared segment files for each
    # chromosome

    # creating a directory
    preformater = pre_shared_segments_analysis_scripts.Pre_Shared_Segment_Converter(
        ibd_files, carrier_file, ibd_program, output, map_file_dir)

    segment_file_list = preformater.gather_segment_files(file_suffix)

    # also need to get all of the possible variant files per chromosome
    chr_var_file_list = preformater.gather_chromosome_files()

    # getting the list of map files
    map_file_list = preformater.get_map_files()
    # Need to make sure that the proper segment file is passed with the proper
    # chromosome file
    for chromo_file in chr_var_file_list:
        match = re.search(r'chr\d\d\.', chromo_file)

        chr_num = match.group(0)

        hypen_chr_num: str = "".join([chr_num.strip("."), "_"])

        # Creating a tuple that gets the proper segment_file and the proper map file that corresponds to that chr
        segment_map_tuple = [
            (segment_file, map_file) for segment_file in segment_file_list
            for map_file in map_file_list
            if chr_num in segment_file and hypen_chr_num in map_file
        ]

        # This checks to see if the tuple is empty or not
        if segment_map_tuple:

            # This gets the values out of the tuple
            # The first element is the segment file
            segment_file = segment_map_tuple[0][0]
            # The second element is the map file
            map_file = segment_map_tuple[0][1]

            var_info_df, variant_directory = preformater.create_variant_lists(
                chromo_file, map_file)

            # This checks if the var_info_file_path and the variant_directory are empty string because
            # this would mean that the the chromo_file only has one variant and it has no carriers
            if variant_directory == "":
                print(f"There were no carriers in the file {chromo_file}")
                continue

            iid_file_list = preformater.get_iid_files(variant_directory)

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
