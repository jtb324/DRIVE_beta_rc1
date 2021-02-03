# importing modules
import pandas as pd
import numpy as np
import multiprocessing as mp
import os
import sys
from functools import partial

# importing module from another file
import analysis_haplotypes


def get_variant_list(haplotype_df: pd.DataFrame) -> list:
    '''This function will return a list of all variants that can be written to a file'''

    variant_list: list = list(set(haplotype_df.variant_id.values.tolist()))

    return variant_list


def subset_df_for_variant(haplotype_df: pd.DataFrame, variant_id: str):
    '''This function will subset the haplotype just for a specific variant'''

    haplotype_subset: pd.DataFrame = haplotype_df[haplotype_df.variant_id ==
                                                  variant_id]

    return haplotype_subset


def write_ids_to_file(pair1_id: str, pair2_id: str, output_path: str) -> str:
    '''This function will write the pair ids to a file that has to be passed
    to the plink subprocess "keep" flag'''

    full_output_path: str = "".join(
        [output_path, pair1_id, "_", pair2_id, "_keep.txt"])

    with open(full_output_path, "w") as keep_file:

        keep_file.write("".join([pair1_id, " ", pair1_id, "\n"]))
        keep_file.write("".join([pair2_id, " ", pair2_id, "\n"]))

    return full_output_path


def read_in_lines(file_path: str) -> tuple:
    '''This function will read in the first and second line of the file to 
    memory from the plink output ped file'''

    with open(file_path, "r") as plink_file:

        line_1: str = plink_file.readline()
        line_2: str = plink_file.readline()

    return (line_1, line_2)


def remove_plink_files(ped_file_path: str):
    """
    This function will delete the plink output files. This includes the ped,
    the map, and the log files.
    """

    suffix_indx: int = ped_file_path.index(".ped")

    file_path: str = ped_file_path[:suffix_indx]

    # generating a suffix list of all the file outputs from running plink
    suffix_list: list = [".ped", ".map", ".log"]

    # removing the files iteratively
    for file_path in ["".join([file_path, suffix]) for suffix in suffix_list]:

        os.remove(file_path)


def run_plink(binary_file: str, output_path: str, haplotype_row_tuple: tuple):
    '''This function will be responsible for actually running plink'''
    print("This is a thread")
    # pulling relevant values out of the haplotype tuple
    pair_1_id: str = haplotype_row_tuple[1]
    pair_2_id: str = haplotype_row_tuple[2]
    chr_num: int = int(haplotype_row_tuple[3])
    hapibd_start_bp: str = haplotype_row_tuple[6]
    hapibd_end_bp: str = haplotype_row_tuple[7]
    ilash_start_bp: str = haplotype_row_tuple[9]
    ilash_end_bp: str = haplotype_row_tuple[10]

    # writing the pairs to a keep file that can then be passed to plink
    keep_file_str: str = write_ids_to_file(
        pair_1_id, pair_2_id, "".join([output_path, "keep_id_files/"]))

    if not np.isnan(hapibd_start_bp) or not np.isnan(hapibd_end_bp):
        hapibd_path: str = analysis_haplotypes.get_plink_haplotype_str(
            binary_file, keep_file_str,
            str(hapibd_start_bp), str(hapibd_end_bp), output_path, "hapibd",
            str(chr_num), pair_1_id, pair_2_id)

    if not np.isnan(ilash_start_bp) or not np.isnan(ilash_end_bp):
        ilash_path: str = analysis_haplotypes.get_plink_haplotype_str(
            binary_file, keep_file_str, str(ilash_start_bp), str(ilash_end_bp),
            output_path, "ilash", str(chr_num), pair_1_id, pair_2_id)

    # remove the keep file since it is only used by plink
    os.remove(keep_file_str)

    # Next section reads in lines from the two plink files if each file was made
    if hapibd_path:

        hapibd_tuple: tuple = read_in_lines(hapibd_path)

        analysis_haplotypes.check_tuple_size(hapibd_tuple)

        remove_plink_files(hapibd_path)

    if ilash_path:

        ilash_tuple: tuple = read_in_lines(ilash_path)

        analysis_haplotypes.check_tuple_size(ilash_tuple)

        remove_plink_files(ilash_path)

    # getting the differences in the strings for hapibd
    string_diff: int = analysis_haplotypes.compare_haplotype_str(hapibd_tuple)

    sys.exit(1)


def parallelize(variant_list: list, workers: int, haplotype_df: pd.DataFrame,
                binary_file: str, output_dir: str):
    '''This function will run the script in parallel'''
    # manager = mp.Manager()

    keep_file_output_path: str = "".join([output_dir, "keep_id_files/"])

    try:
        os.mkdir(keep_file_output_path)

    except FileExistsError:
        pass
    pool = mp.Pool(int(workers))

    # This next line allows me to pass multiple arguments into the map function
    # below
    func = partial(run_plink, binary_file, output_dir)
    # This next line will iterate through the variant list so that the program
    # can subset the haplotype dataframe
    for variant in variant_list:

        haplotype_df_subset: pd.DataFrame = subset_df_for_variant(
            haplotype_df, variant)

        pool.map(func, haplotype_df_subset.itertuples(name=None))

    pool.close()

    pool.join()

    # TODO: work on finishing parallelizing this script


def main_run(binary_file: str, output_path: str, haplotype_info_file: str,
             workers: int):
    '''This is the main function responsible for running all the functions
    associated with analysing the haplotypes. It will accept arguments from the
    MEGA_ID.py file'''

    # making the output directory that forms files containing a list of ids to keep
    # reading in the haplotype file into a dataframe
    haplotype_df: pd.DataFrame = pd.read_csv(haplotype_info_file, sep="\t")

    variant_list: list = get_variant_list(haplotype_df)

    print(output_path)

    parallelize(variant_list, workers, haplotype_df, binary_file, output_path)
