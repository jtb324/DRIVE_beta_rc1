# This script will get the haplotype information
import pandas as pd
import argparse
import os
from os import path
import glob
import re
import sys
import gzip


def identify_unique_variants(confirmed_carrier_file: str) -> list:
    '''This function will get all the unique variants that we have a confirmed carrier for'''

    # loads the variant_id column from the file into dataframe
    confirmed_carrier_df: pd.DataFrame = pd.read_csv(
        confirmed_carrier_file, sep="\t", usecols=["variant_id"])

    # gets a list of all the unique variants in the file
    unique_var_list: list = confirmed_carrier_df.variant_id.unique().tolist()
    print(len(unique_var_list))
    return unique_var_list


def get_file_list(file_dir: str, file_tag: str) -> list:

    # getting teh current directory
    cur_dir: str = os.getcwd()

    # changing to the file directory
    os.chdir(file_dir)

    # creating a list for the files to be returned in
    file_list = []

    for file in glob.glob(file_tag):

        full_file_path = "".join([file_dir, file])

        file_list.append(full_file_path)

    os.chdir(cur_dir)

    return file_list


def get_chr_id(allpair_file: str) -> str:
    '''This will get the chromosome number from the allpair file name'''

    match = re.search(r'_chr\d\d\.', allpair_file)

    # find chromosome number
    if match:

        chr_num: str = match.group(0)

    else:
        match = re.search(r'_chr\d.', allpair_file)

        chr_num: str = match.group(0)

    chr_num = chr_num[1:-1]

    return chr_num


def get_file(file_list: list, identifier: str) -> str:
    '''This function gets the file that matches a condition from a list of files'''

    file_str: str = [file for file in file_list if identifier in file][0]

    return file_str


def get_var_pos(map_file: str, variant_id: str) -> str:
    '''This function find the position of the variant id and return that as a string'''

    # This reads in the columns from the map file that contains the variant probe id and the pos
    map_df: pd.DataFrame = pd.read_csv(
        map_file, sep="\t", header=None, usecols=[1, 3])

    # getting the value
    pos: str = str(map_df[map_df[1] == variant_id][3].values[0])

    return pos


def check_ibd_program(file_str: str) -> int:
    '''This function will return a zero or a 1 depending on whether or not the file is a hapibd file or an ilash file respectively'''

    # create a dictionary to handle the situation
    ibd_handler: dict = {
        True: 1,
        False: 0
    }

    # Checks to see if the file is a ilash file
    handler_int: int = ibd_handler["match" in file_str]

    return handler_int


def get_index_positions(ibd_program_int: int) -> list:

    # This dictionary contains the index of values such as IID1, IID2, Start point
    # end point, and the shared segment length. True will contain indixes for an iLash
    # file and false will contain the indices for a hapibd file
    indx_dict: dict = {
        1: [0, 2, 5, 6, 9],
        0: [0, 2, 5, 6, 7]

    }

    indx_list: list = indx_dict[ibd_program_int]

    return indx_list


def filter_file(file: str, variant_pos: str) -> list:
    '''Filtering the ilash or hapibd file for the specific variant'''

    pair_list: list = []

    # figure out if the file is an ilash file or hapibd
    ibd_indicator: int = check_ibd_program(file)

    # Getting the index positions of the information based off of the ibd_indicator integer
    indx_list = get_index_positions(ibd_indicator)

    # expanding the index list
    pair_1_indx: int = indx_list[0]
    pair_2_indx: int = indx_list[1]
    start_indx: int = indx_list[2]
    end_indx: int = indx_list[3]
    segment_length_indx: int = indx_list[4]

    # figure out if file is ilash or hapibd
    with gzip.open(file, "rt") as ibd_file:

        for row in ibd_file:

            split_row: list = row.split()

            # This line checks to see if the variant position falls inbetween the start and end position of the row
            if int(split_row[start_indx]) <= int(variant_pos) and int(split_row[end_indx]) >= int(variant_pos):
                # if it matches then it will append that value to the pair_list
                pair_list.append(row)

            # If not it moves onto the next row
            else:
                continue

    return pair_list


def get_carriers(carrier_file: str, variant_id: str) -> list:
    '''This function will get a list of the carriers identified for a specific variant'''
    # REading the carrier file into a dataframe
    carrier_df: pd.DataFrame = pd.read_csv(carrier_file)

    # subsetting to get only the variant of interest carriers
    carrier_list: list = carrier_df[carrier_df["Variant ID"]
                                    == variant_id]["IID"].values.tolist()

    return carrier_list


def filter_for_carriers(ibd_pair_list: list, carrier_list: list) -> list:
    '''This function will filter the pair list for only pairs where both individual are carriers'''
    print(carrier_list)

    filtered_list: list = []
    # iterating through each pair string in the pair_list
    for pair_str in ibd_pair_list:

        # removing the newline
        pair_str = pair_str.strip("\n")

        # splitting the row by tab
        pair_list: list = pair_str.split("\t")

        # These next two lines get te iid1 and iid2
        iid1: str = pair_list[0]

        iid2: str = pair_list[2]

        if iid1 in carrier_list and iid2 in carrier_list:

            filtered_list.append(pair_str)

    return filtered_list


def write_to_file(hapibd_filtered_list: str, ilash_filtered_list: str, output_dir: str, chr_num: str, variant_id: str):
    '''This function will write the pairs and the segment length and the segment start and end to a file'''

    # get the chromosome number digit to be used later
    chr_digit: str = re.findall(r'[0-9]+', chr_num)[0]

    output: str = "".join([output_dir, "haplotype_info.txt"])
    # Opening the file
    with open(output, "a+") as output_file:

        # This will check if the file has any data written in it
        if os.path.getsize(output) == 0:

            output_file.write(
                f"pair_1\tpair_2\tchr\tvariant_id\thapibd_start\thapibd_end\thapibd_len\tilash_start\tilash_end\tilash_len\n")

        # iterating through each string in the filtered list and then matching that with the other filtered list
        for pair_string in hapibd_filtered_list:

            # splitting the row to get values of interest
            split_hapibd_list: list = pair_string.split("\t")
            # getting the pairs from the split hapibd_string
            pair_1: str = split_hapibd_list[0]
            pair_2: str = split_hapibd_list[2]
            print(pair_1, pair_2)
            # using list comprehension to get the string from the ilash_filtered_list that includes the same rows
            try:
                ilash_string: str = list(filter(lambda string: (pair_1
                                                                in string and pair_2 in string), ilash_filtered_list))[0]
            except IndexError:
                continue

            split_ilash_list: list = ilash_string.split("\t")

            output_file.write(
                f"{split_hapibd_list[0]}\t{split_hapibd_list[2]}\t{chr_digit}\t{variant_id}\t{split_hapibd_list[5]}\t{split_hapibd_list[6]}\t{split_hapibd_list[7]}\t{split_ilash_list[5]}\t{split_ilash_list[6]}\t{split_ilash_list[9]}\n")


def get_full_var_name(allpair_file: str, variant_id: str) -> str:
    '''This function will get the full variant name out of the allpair file name'''

    variant_indx: int = allpair_file.index(variant_id)

    # This is the ending index of the slice
    second_indx: int = variant_indx+len(variant_id)+2

    full_variant_id: str = allpair_file[variant_indx:second_indx]

    print(full_variant_id)

    return full_variant_id


def run(args):
    "function to run"
    # This section will check if the output file exist from a previous run and if it does then it will delete it
    if path.exists("".join([args.output, "haplotype_info.txt"])):

        os.remove("".join([args.output, "haplotype_info.txt"]))

    # getting a list of all the variants that have a confirmed carrier
    variant_list: list = identify_unique_variants(args.variant_file)

    # Getting all the carrier files
    carrier_file_list: list = get_file_list(args.carrier, "*.csv")

    # getting the map_files
    map_file_list: list = get_file_list(args.map_dir, "*.map")

    # getting the list of allpair files
    allpair_file_list: list = get_file_list(args.allpair_dir, "*.allpair.txt")

    for variant in variant_list:
        # getting the correct allpair file
        allpair_file: str = get_file(allpair_file_list, variant)

        print(allpair_file)

        chr_num: str = get_chr_id(allpair_file)

        full_variant_id: str = get_full_var_name(allpair_file, variant)

        map_file: str = get_file(map_file_list, "".join([".", chr_num, "_"]))
        print(map_file)
        # getting the specific carrier file for the chromosome
        carrier_file: str = get_file(
            carrier_file_list, "".join([chr_num, "_"]))

        carriers_list: list = get_carriers(carrier_file, full_variant_id)

        print(f"The map file is {map_file}")

        var_pos: str = get_var_pos(map_file, variant)

        # The next four lines get the hapibd file and the ilash file for the specific chromosome
        ilash_file_list: list = get_file_list(args.ilash_dir, "*.match.gz")

        hapibd_file_list: list = get_file_list(args.hapibd_dir, "*.ibd.gz")

        ilash_file: str = get_file(
            ilash_file_list, "".join(["_", chr_num, "."]))

        hapibd_file: str = get_file(
            hapibd_file_list, "".join(["_", chr_num, "."]))

        # getting a list of each pair
        ilash_list: list = filter_file(ilash_file, var_pos)

        ilash_carriers_list: list = filter_for_carriers(
            ilash_list, carriers_list)

        hapibd_list: list = filter_file(hapibd_file, var_pos)

        hapibd_carriers_list: list = filter_for_carriers(
            hapibd_list, carriers_list)

        write_to_file(hapibd_carriers_list, ilash_carriers_list,
                      args.output, chr_num, full_variant_id)


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument("-v", help="This argument takes the directory of the file that list all the confirmed carriers",
                        dest="variant_file", type=str, required=True)

    parser.add_argument("-m", help="This argument list the directory for the map files",
                        dest="map_dir", type=str, required=True)

    parser.add_argument("-a", help="This argument list the directory for all of the allpair files",
                        dest="allpair_dir", type=str, required=True)

    parser.add_argument("-hap", help="This argument list the directory for all the hapibd files",
                        dest="hapibd_dir", type=str, required=True)

    parser.add_argument("-ilash", help="This argument list the directory for all of the ilash files",
                        dest="ilash_dir", type=str, required=True)

    parser.add_argument("-c", help="This argument list the directory for the reformated carrier files",
                        dest="carrier", type=str, required=True)

    parser.add_argument("-o", help="This argument list the output file path for the final file",
                        dest="output", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
