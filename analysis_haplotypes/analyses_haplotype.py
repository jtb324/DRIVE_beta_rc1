# importing modules
import pandas as pd
import multiprocessing as mp
import os
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
                                                  "variant_id"]

    return haplotype_subset


def write_ids_to_file(pair1_id: str, pair2_id: str, output_path: str) -> str:
    '''This function will write the pair ids to a file that has to be passed
    to the plink subprocess "keep" flag'''

    full_output_path: str = "".join(
        [output_path, pair1_id, "_", pair2_id, "_keep.txt"])

    with open(output_path, "w") as keep_file:

        keep_file.write("".join([pair1_id, " ", pair1_id, "\n"])
        keep_file.write("".join([pair2_id, " ", pair2_id, "\n"])

    return full_output_path

def run_plink(binary_file: str, output_path: str, haplotype_row_tuple: tuple):
    '''This function will be responsible for actually running plink'''

    # pulling relevant values out of the haplotype tuple
    pair_1_id: str=haplotype_row_tuple[0]
    pair_2_id: str=haplotype_row_tuple[1]
    chr_num: int=int(haplotype_row_tuple[2])
    hapibd_start_bp: str=haplotype_row_tuple[5]
    hapibd_end_bp: str=haplotype_row_tuple[6]
    ilash_start_bp: str=haplotype_row_tuple[8]
    ilash_end_bp: str=haplotype_row_tuple[9]

    # writing the pairs to a keep file that can then be passed to plink
    keep_file_str: str=write_ids_to_file(pair_1_id, pair_2_id, output_path)

    if hapibd_start_bp != "N/A" or hapibd_end_bp != "N/A":
        analysis_haplotypes.get_pair_haplotype_str(
            binary_file, keep_file_str, hapibd_start_bp, hapibd_end_bp, output_path, "hapibd", chr_num)

    if ilash_start_bp != "N/A" or ilash_end_bp != "N/A":
                        analysis_haplotypes.get_pair_haplotype_str(
            binary_file, keep_file_str, ilash_start_bp, ilash_end_bp, output_path, "ilash", chr_num)


def parallelize(variant_list: list, workers: int, haplotype_df: pd.DataFrame,
                binary_file: str, output_dir: str):
    '''This function will run the script in parallel'''
    # manager = mp.Manager()

    keep_file_output_path: str="".join([output_path, "keep_id_files/"])

    try:
        os.mkdir(keep_file_output_path)

    except FileExistsError:
        pass

    pool=mp.Pool(workers)

    # This next line allows me to pass multiple arguments into the map function below
    func=partial(run_plink, binary_file, output_dir)
    # This next line will iterate through the variant list so that the program
    # can subset the haplotype dataframe
    for variant in variant_list:

        haplotype_df_subset: pd.DataFrame=subset_df_for_variant(
            haplotype_df, variant)

        pool.map(func, haplotype_df_subset.itertuples(name=False))

        pool.close()

        pool.join()

    # TODO: work on finishing parallelizing this script


            def main_run(binary_file: str, output_path: str, haplotype_info_file: str, workers: int):
    '''This is the main function responsible for running all the functions
    associated with analysing the haplotypes. It will accept arguments from the
    MEGA_ID.py file'''

    # making the output directory that forms files containing a list of ids to keep
    # reading in the haplotype file into a dataframe
    haplotype_df: pd.DataFrame=pd.read_csv(haplotype_info_file, sep="\t")

    variant_list: list=get_variant_list(haplotype_df)

    parallelize(variant_list, workers, haplotype_df, binary_file, output_dir)
