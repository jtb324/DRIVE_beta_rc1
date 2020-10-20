import pandas as pd
import glob
import os
from os import path
import argparse


def get_files(file_dir: str, file_extension) -> list:
    '''This function get the files for the directory'''

    cur_dir: str = os.getcwd()

    os.chdir(file_dir)

    file_list = []

    for file in glob.glob(file_extension):

        full_file_path = "".join([file_dir, file])

        file_list.append(full_file_path)

    os.chdir(cur_dir)

    return file_list


# def get_unique_variants(pairs_df: pd.DataFrame) -> list:
#     '''This function will return a list of all the unique variants'''

#     # the following line gets a list of all the unique variants within the dataframe
#     variant_list: list = pairs_df.variant_id.values.unique().tolist()

#     return variant_list

# def get_unique_network_ids(carriers_df_subset: pd.DataFrame) -> list:
#     '''This function will take a subset of the carriers_df only for a specific variant and
#     return a list of all the unique network ids'''

#     # the following line gets a list of all the unique variants within the dataframe
#     id_list: list = carriers_df_subset["Network ID"].values.unique(
#     ).tolist()

#     return id_list


def get_index_positions(file_path: str) -> list:

    # This dictionary contains the index of values such as IID1, IID2, Start point
    # end point, and the shared segment length. True will contain indixes for an iLash
    # file and false will contain the indices for a hapibd file
    indx_dict: dict = {
        True: [0, 2, 5, 6, 9],
        False: [0, 2, 5, 6, 7]

    }

    indx_list: list = indx_dict["match" in file_path]

    return indx_list


# def get_pair_list(df_subset: pd.DataFrame, col_num: int) -> list:
#     '''This function will take the df subset and and specific column such as 1 or 2 and return
#     a list of all the pairs for that column'''

#     col_handler: dict = {
#         1: "pair_1",
#         2: "pair_2"
#     }

#     # This next line gets the name of the column
#     col_name: str = col_handler[col_num]

#     # This next line gets a list of either all the pair 1 or pair 2 values
#     pair_list: list = df_subset[col_name].values.tolist()

#     return pair_list


# def get_chr_num(dataframe: pd.DataFrame) -> str:
#     '''This function gets the chromosome number for the specific variant and network id'''

#     # This gets the unique values of the chromosome numbers
#     chr_num: str = str(dataframe.chr_num.values.unique().tolist())

#     return chr_num


def get_ibd_file(file_list: list, chr_str: str) -> str:
    '''This function gets the ibd files using list comprehension'''

    ibd_file: str = [file for file in file_list if chr_str in file][0]

    return ibd_file


def get_haplotype_info(ibd_file: str, pair_1: str, pair_2: str) -> dict:
    '''This function will get the haplotype information for the pair in the ibd file'''
    indx_list: list = get_index_positions(ibd_file)
    # Making a dictionary that will hold the positions of interest for each variant
    info_dict: dict = dict()

    # Getting the indices for the positions of interest
    pair_1_indx: int = indx_list[0]
    pair_2_indx: int = indx_list[1]
    start_indx: int = indx_list[2]
    end_indx: int = indx_list[3]
    segment_length_indx: int = indx_list[4]

    # iterating through chunks in the ibd file to avoid memory issues
    for chunk in pd.read_csv(ibd_file, sep="\t", chunksize=100000, header=None):

        # getting the chunk where the pair one and pair2 are in the same row
        pair_row_df: pd.DataFrame = chunk[(chunk[pair_1_indx] == pair_1) & (
            chunk[pair_2_indx] == pair_2)]

        # getting the start end and segment length
        start: int = pair_row_df[start_indx].values
        end: int = pair_row_df[end_indx].values
        segment_length: float = pair_row_df[segment_length_indx].values

        # putting the values of interest into a dictionary for that pair
        info_dict["start"] = start
        info_dict["end"] = end
        info_dict["length"] = segment_length

        # these lines break out of the for loop if it has found the row that contains the pair
        if not pair_row_df.empty:
            print(pair_row_df)
            break
    print(info_dict)
    return info_dict


# def add_values_to_dict(pairs_df: pd.DataFrame, pair_1: str, pair_2: str, hapibd_dict: dict, ilash_dict: dict) -> dict:
#     '''This function will add the values of interest to the correct row in the dataframe'''

#     # Getting the values of interest from each dictionary
#     hapibd_values_list: list = hapibd_dict[(pair_1, pair_2)]

#     # getting the values out of the list
#     hapibd_start: int = hapibd_values_list[0]
#     hapibd_end: int = hapibd_values_list[1]
#     hapibd_length: int = hapibd_values_list[2]

#     # Getting the values of interest from each dictionary
#     ilash_values_list: list = ilash_dict[(pair_1, pair_2)]

#     # getting the values out of the list
#     ilash_start: int = ilash_values_list[0]
#     ilash_end: int = ilash_values_list[1]
#     ilash_length: int = ilash_values_list[2]

#     return pairs_dict
def write_to_file(pairs_dict: dict, hapibd_dict: dict, ilash_dict: dict, output: str):

    # This creates the full output path by combining the provided output argument with the name of the file
    full_output_path = "".join([output, "haplotype_info.txt"])

    file_header = "pair_1\tpair_2\tnetwork_id\tvariant_id\tchr_num\thapibd_start\thapibd_end\thapibd_length\tilash_start\tilash_end\tilash_length\n"

    with open(full_output_path, "a+") as output_file:

        if os.path.getsize(full_output_path) == 0:

            output_file.write(file_header)

        # check if the hapibd dictionary is empty
        if hapibd_dict["start"].size == 0 and ilash_dict["start"].size > 0:
            output_file.write(
                f"{pairs_dict['pair_1']}\t{pairs_dict['pair_2']}\t{pairs_dict['network_id']}\t{pairs_dict['variant_id']}\t{pairs_dict['chr']}\t{'Nan'}\t{'Nan'}\t{'Nan'}\t{ilash_dict['start'][0]}\t{ilash_dict['end'][0]}\t{ilash_dict['length'][0]}\n")

        elif hapibd_dict["start"].size > 0 and ilash_dict["start"].size == 0:
            output_file.write(
                f"{pairs_dict['pair_1']}\t{pairs_dict['pair_2']}\t{pairs_dict['network_id']}\t{pairs_dict['variant_id']}\t{pairs_dict['chr']}\t{hapibd_dict['start'][0]}\t{hapibd_dict['end'][0]}\t{hapibd_dict['length'][0]}\t{'Nan'}\t{'Nan'}\t{'Nan'}\n")
        else:
            output_file.write(f"{pairs_dict['pair_1']}\t{pairs_dict['pair_2']}\t{pairs_dict['network_id']}\t{pairs_dict['variant_id']}\t{pairs_dict['chr']}\t{hapibd_dict['start'][0]}\t{hapibd_dict['end'][0]}\t{hapibd_dict['length'][0]}\t{ilash_dict['start'][0]}\t{ilash_dict['end'][0]}\t{ilash_dict['length'][0]}\n")

#This is the main function to run the script ###############################################


def run(args):
    "function to run"

    # getting the list of ilash files

    ilash_file_list: list = get_files(args.ilash, "*match.gz")

    # getting the list of hapibd files
    hapibd_file_list: list = get_files(args.hapibd, "*ibd.gz")

    # deleting the outuput file if it already exist
    if path.exists("".join([args.output, "haplotype_info.txt"])):

        os.remove("".join([args.output, "haplotype_info.txt"]))

    with open(args.pairs, "r") as pairs_file:

        # skipping the first row of the pairs_file
        next(pairs_file)

        # iterating through the file
        for pair_row in pairs_file:

            # splitting the row
            row_list: list = pair_row.split(",")

            pair_1: str = row_list[0]
            pair_2: str = row_list[1]

            network_id: str = row_list[2]
            variant_id: str = row_list[3]
            # # Getting a list of all the unique variants
            # unique_var_list: list = get_unique_variants(
            #     carrier_chunk)

            # # iterating through each variant
            # for variant in unique_var_list:

            #     # Restricting the carrier dataframe to only those variants
            #     carrier_df_subset: pd.DataFrame = carrier_chunk[carrier_chunk.variant_id == variant]

            #     # Getting a list of all the network ids
            #     unique_network_ids: list = get_unique_network_ids(
            #         carrier_df_subset)

            # for network_id in unique_network_ids:

            #     # subsetting the dataframe further just for individuals with a certain network id
            #     in_network_df: pd.DataFrame = carrier_df_subset[carrier_df_subset["Network ID"] == network_id]

            # getting all of the first pairs in a list
            # pair_1_list: list = get_pair_list(in_network_df, 1)

            # # getting all of the second pairs in a list
            # pair_2_list: list = get_pair_list(in_network_df, 2)

            # getting the chromosome number
            chr_num: str = row_list[4].strip("\n")

            # adding chr and . to the chr_num so that it matches the files
            chr_str: str = "".join(["_chr", chr_num, "."])
            print(chr_str)
            # getting the specific ilash file for that chromosome
            ilash_ibd_file: str = get_ibd_file(ilash_file_list, chr_str)
            print(ilash_ibd_file)
            # getting the specific hapibd file for that chromosome
            hapibd_ibd_file: str = get_ibd_file(hapibd_file_list, chr_str)
            print(hapibd_ibd_file)
            # # iterating through the pairs list
            # for pair_1, pair_2 in zip(pair_1_list, pair_2_list):

            # getting values of interest from the hapibd file for the row
            hapibd_values_dict: dict = get_haplotype_info(
                hapibd_ibd_file, pair_1, pair_2)

            # getting the values of interest from the iLash files
            ilash_values_dict: dict = get_haplotype_info(
                ilash_ibd_file, pair_1, pair_2)

            # Putting the values from the original file into a dictionary so that it can be easily passed to the function
            pair_values_dict = {
                "pair_1": pair_1,
                "pair_2": pair_2,
                "network_id": network_id,
                "variant_id": variant_id,
                "chr": chr_num
            }
            # write these values to a new file
            write_to_file(pair_values_dict, hapibd_values_dict,
                          ilash_values_dict, args.output)


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument("--ilash", help="This argument list the path to the directory of ilash files",
                        dest="ilash", type=str, required=True)

    parser.add_argument("--hapibd", help="This argument list the path to the directory of hapibd files",
                        dest="hapibd", type=str, required=True)

    parser.add_argument("--pairs", help="This argument list the path to the file containing a list of pairs",
                        dest="pairs", type=str, required=True)

    parser.add_argument("--output", help="This argument list the output directory",
                        dest="output", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
