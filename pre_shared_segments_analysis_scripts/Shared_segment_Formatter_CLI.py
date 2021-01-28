#!/usr/bin/python

# THese are modules used
import re
from numpy.lib.function_base import msort
import pandas as pd
import sys
import os
from os import path
import shutil
import multiprocessing as mp
from functools import partial

import pre_shared_segments_analysis_scripts
import utility_scripts

####################################################################################################

# Initializing the parser


def listener(que_object, output: str, header: str):
    '''This function will listen to the que and then write the element of the que to a file'''

    # opening the output file to write to
    with open(output, "a+") as output_file:

        # checking if the file size is zero
        if os.path.getsize(output) == 0:

            output_file.write(header)

        while 1:

            m = que_object.get()

            if m == "kill":

                break

            output_file.write(m)
            output_file.flush()


def remove_previous_file(file_path: str):
    '''This function will remove previous output files from previous runs'''
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


def create_readme(output: str):
    '''This function creates a readme file in the specified directory'''
    readme = utility_scripts.Readme("_README.md", output)
    readme.rm_previous_file()
    readme.write_header("formatted_ibd_output/")
    readme.create_date_info()
    readme.add_line(utility_scripts.formatted_ibd_dir_body_text_1)
    readme.add_line(utility_scripts.formatted_ibd_dir_body_text_2)


def convert_ibd(ibd_files: str, carrier_file: str, ibd_program: str,
                output: str, map_file_dir: str, file_suffix: str, min_CM: str,
                threads: str):

    ###########################################################
    # This first section will be used to get the shared segment files for each chromosome
    create_readme(output)
    # creating a directory
    preformater = pre_shared_segments_analysis_scripts.Pre_Shared_Segment_Converter(
        ibd_files, carrier_file, ibd_program, output, map_file_dir)

    segment_file_list = preformater.gather_segment_files(file_suffix)

    # also need to get all of the possible variant files per chromosome
    chr_var_file_list = preformater.gather_chromosome_files()

    # getting the list of map files
    map_file_list = preformater.get_map_files()
    # Need to make sure that the proper segment file is passed with the proper chromosome file
    for chromo_file in chr_var_file_list:
        match = re.search(r'chr\d\d\.', chromo_file)

        chr_num = match.group(0)
        if (chr_num.find("0") < len(chr_num) - 2) and chr_num.find("0") != -1:

            zero_indx: str = chr_num.find("0")

            alt_chr_num: str = "".join(
                [chr_num[:zero_indx], chr_num[zero_indx + 1:]])

        hypen_chr_num: str = "".join([chr_num.strip("."), "_"])

        # Creating a tuple that gets the proper segment_file and the proper map file that corresponds to that chr
        segment_map_tuple = [
            (segment_file, map_file) for segment_file in segment_file_list
            for map_file in map_file_list
            if chr_num in segment_file and hypen_chr_num in map_file
        ]

        if not segment_map_tuple:
            segment_map_tuple = [
                (segment_file, map_file) for segment_file in segment_file_list
                for map_file in map_file_list
                if (chr_num in segment_file or alt_chr_num in segment_file)
                and hypen_chr_num in map_file
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

            run_parallel(segment_file, output, ibd_program, min_CM,
                         iid_file_list, var_info_list, threads)


def run_parallel(segment_file: str, output_path: str, ibd_format: str,
                 min_CM: str, iid_file_list: list, variant_info_list: list,
                 threads: int):

    # starting a que to create a file that keeps track of errors
    manager = mp.Manager()

    que = manager.Queue()

    header: str = "variant_id"

    pool = mp.Pool(int(threads))

    watcher = pool.apply_async(
        listener,
        (que, "".join([output_path, "nopairs_identified.txt"]), header))

    func = partial(run_main, segment_file, output_path, ibd_format, min_CM,
                   iid_file_list, que)

    pool.map(func, variant_info_list)

    que.put("kill")

    pool.close()

    pool.join()


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
