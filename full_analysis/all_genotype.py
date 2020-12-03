import pandas as pd
import glob
import os
from os import path

import population_filter_scripts

# functions to gather all necessary files


def gather_files(file_directory: str, file_tag: str) -> list:
    '''This function will gather all of the allpair.new.txt files which contain information about pairs. It will also be used to get the 'chr#_list.single_variant.csv' files.'''
    curr_dir: str = os.getcwd()

    os.chdir(file_directory)

    file_list = []

    for file in glob.glob(file_tag):

        full_file_path = "".join([file_directory, file])

        file_list.append(full_file_path)

    os.chdir(curr_dir)

    return file_list


def get_variants(confirmed_carriers_df: pd.DataFrame) -> list:
    '''This function gets the variants from the provided confirmed_carriers_df'''

    variants_list: list = confirmed_carriers_df.variant_id.tolist()

    return list(set(variants_list))


def get_chr_num(confirmed_carrier_df: pd.DataFrame, variant_id: str) -> str:
    '''This returns the chromosome number for a specific variant'''

    chr_num: str = confirmed_carrier_df[confirmed_carrier_df.variant_id == variant_id]["chr"].values.tolist()[
        0]

    return chr_num


def identify_file(file_list: list, chr_num: str) -> str:
    '''This function will identify the files of interest'''
    print(file_list)
    print(chr_num)
    file_name: str = [file for file in file_list if chr_num in file][0]

    return file_name


def get_variant_map_index(map_df: pd.DataFrame, variant: str) -> int:
    '''This function will return the row index of the variant within the map file'''

    # subset map_df for just the specific variant
    map_df_subset: pd.DataFrame = map_df[map_df[1] == variant]

    # getting the index value
    var_indx: int = int(list(map_df_subset.index)[0])

    return var_indx


def get_genotype_df(ped_df: pd.DataFrame, var_indx: int) -> pd.DataFrame:
    '''This will return a dataframe that has all the carriers for  a specific variant'''

    # adjusting the index for the first six columns in the ped file
    first_indx: int = var_indx*2 + 6
    second_indx: int = var_indx*2 + 7

    # subsetting the ped file to only the IID, and then the two genotypes for the IID
    ped_file_subset: pd.DataFrame = ped_df[[1, first_indx, second_indx]]

    ped_file_subset = ped_file_subset.rename(
        columns={1: "IID", first_indx: "allele_1", second_indx: "allele_2"})

    ped_file_subset["genotype"] = ped_file_subset["allele_1"] + \
        ped_file_subset["allele_2"]

    # dropping the allele_columns
    ped_file_genotypes = ped_file_subset.drop(["allele_1", "allele_2"], axis=1)

    return ped_file_genotypes


def write_to_file(output: str, genotype_df: pd.DataFrame):

    genotype_df.to_csv(
        output, index=False, mode="a", sep="\t")


def getting_confirmed_carrier_df_subset(confirmed_carrier_df: pd.DataFrame, variant: str) -> pd.DataFrame:
    '''This function will return the subset of the dataframe for the specified variant'''

    df_subset: pd.DataFrame = confirmed_carrier_df[confirmed_carrier_df.variant_id == variant]

    return df_subset


def get_IIDs(confirmed_carriers_df_subset: pd.DataFrame) -> list:
    '''This function will return a list of all the IIDs in the confirmed_carriers_df_subset'''

    return confirmed_carriers_df_subset.IID.tolist()


def getting_the_other_iids(confirmed_IID_list: list, ped_file_genotype_df: pd.DataFrame) -> pd.DataFrame:
    '''This function returns a dataframe that has all the IIDs that are not in the confirmed_carriers.txt file'''

    filter_df: pd.DataFrame = ped_file_genotype_df[~ped_file_genotype_df.IID.isin(
        confirmed_IID_list)]

    print(confirmed_IID_list)
    print(filter_df)

    return filter_df


def merge_df(confirmed_carrier_subset: pd.DataFrame, genotype_df: str) -> pd.DataFrame:
    '''This function merges the two dataframes to get all the iids for each variant'''

    return pd.concat([confirmed_carrier_subset, genotype_df])


def reformat_chr_num(chr_num) -> str:
    '''This function reformats the chromosome number so that it is only the digits'''

    return chr_num[4:6]


def get_all_genotypes(ped_file_path: str, confirmed_carrier_file: str, population_info_file: str, output_path: str, pop_code: str = None):
    '''This function will be the main run script for this script'''

    full_output_path: str = "".join([output_path, "all_genotypes.txt"])

    # checking if the file exist from a previous run
    if path.exists(full_output_path):

        os.remove(full_output_path)

    map_file_list: list = gather_files(ped_file_path, "*.map")

    ped_file_list: list = gather_files(ped_file_path, "*.ped")

    # load the confirmed carriers file
    confirmed_carriers_df: pd.DataFrame = pd.read_csv(
        confirmed_carrier_file, sep="\t")

    # getting the list of variants from the confirmed_carriers_df
    confirmed_variants_list: list = get_variants(confirmed_carriers_df)

    for variant in confirmed_variants_list:

        chr_num: str = get_chr_num(confirmed_carriers_df, variant)

        # format the chr_num to match the map files
        chr_num = "".join([".", "chr", str(chr_num).zfill(2), "_"])

        # getting the correct map file for the variant and the correct ped file
        map_file: str = identify_file(map_file_list, chr_num)

        ped_file: str = identify_file(ped_file_list, chr_num)

        # load the files into dataframes

        map_df: pd.DataFrame = pd.read_csv(map_file, sep="\t", header=None)

        ped_df: pd.DataFrame = pd.read_csv(ped_file, sep=" ", header=None)

        var_indx: str = get_variant_map_index(map_df, variant)

        ped_file_genotypes: pd.DataFrame = get_genotype_df(ped_df, var_indx)

        # filtering the genotypes down to only the specified pop_code
        if pop_code:

            print(f"this is the pop code: {population_info_file}")

            dataset_filter = population_filter_scripts.Pop_Filter(
                population_info_file, ped_file_genotypes)

            pop_info_df, recode_df = dataset_filter.load_files()

            pop_info_subset_df = dataset_filter.get_pop_info_subset(
                pop_info_df, pop_code)

            ped_file_genotypes = dataset_filter.filter_recode_df(
                pop_info_subset_df, recode_df)

        # getting the subset of the confirmed_carrier_df for the specific variant
        confirmed_carrier_subset: pd.DataFrame = getting_confirmed_carrier_df_subset(
            confirmed_carriers_df, variant)

        # getting the list of iids in the above dataframe
        iid_list: list = get_IIDs(confirmed_carrier_subset)

        # subsetting the ped_file_genotypes to the iids not in the above list
        ped_file_genotypes = getting_the_other_iids(
            iid_list, ped_file_genotypes)

        ped_file_genotypes["confirmed_status"] = "N/A"

        ped_file_genotypes["variant_id"] = variant

        # reformat the chromosome number
        chr_num = reformat_chr_num(chr_num)

        ped_file_genotypes["chr"] = chr_num

        # merging the dataframe

        merged_df: pd.DataFrame = merge_df(
            confirmed_carrier_subset, ped_file_genotypes)

        # writing the genotype to a file
        write_to_file(full_output_path, merged_df)
