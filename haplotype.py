# This script will get the haplotype information
import pandas as pd
import argparse
import os
import glob
import re
import sys
import gzip

# # creating a class for the for the binary trees to keep the first filtered rows
# class Node:
#     '''This class will help to create the binary tree that contains the filtered file strings'''

#     def __init__(self, data:str):
#         self.left = None
#         self.right = None
#         self.data = data
#         self.left_visited = False
#         self.right_visited = False

#     def insert(self, data):
#         if self.data:
#             if self.left == None:
#                 self.left = Node(data)
#                 self.left_visited = True
#             else:
#                 self.left.insert(data)
#             if self.right == None and self.left_visited == False


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

    print(file_list)
    print(identifier)
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
            print(split_row)
            print(int(split_row[start_indx]))
            # This line checks to see if the variant position falls inbetween the start and end position of the row
            if int(split_row[start_indx]) <= int(variant_pos) and int(split_row[end_indx]) >= int(variant_pos):
                # if it matches then it will append that value to the pair_list
                pair_list.append(row)

            # If not it moves onto the next row
            else:
                continue

    return pair_list

# TODO: Need to create another function that can further filter the list so that it contains only pairs were both grids are carriers


def run(args):
    "function to run"

    # getting the map_files
    map_file_list: list = get_file_list(args.map_dir, "*.map")

    # getting the list of allpair files
    allpair_file_list: list = get_file_list(args.allpair_dir, "*.allpair.txt")

    # getting the correct allpair file
    allpair_file: str = get_file(allpair_file_list, args.variant)

    print(allpair_file)

    chr_num: str = get_chr_id(allpair_file)

    map_file: str = get_file(map_file_list, "".join([".", chr_num, "_"]))

    print(f"The map file is {map_file}")

    var_pos: str = get_var_pos(map_file, args.variant)

    # The next four lines get the hapibd file and the ilash file for the specific chromosome
    ilash_file_list: list = get_file_list(args.ilash_dir, "*.match.gz")

    hapibd_file_list: list = get_file_list(args.hapibd_dir, "*.ibd.gz")

    ilash_file: str = get_file(ilash_file_list, "".join(["_", chr_num, "."]))

    hapibd_file: str = get_file(hapibd_file_list, "".join(["_", chr_num, "."]))

    # getting a list of each pair
    ilash_list: list = filter_file(ilash_file, var_pos)

    print(ilash_list)
    print(len(ilash_list))
    print(sys.getsizeof(ilash_list))

    hapibd_list: list = filter_file(hapibd_file, var_pos)

    print(hapibd_list)
    print(len(hapibd_list))
    print(sys.getsizeof(hapibd_list))


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument("-v", help="This argument takes the id of the variant of interest",
                        dest="variant", type=str, required=True)

    parser.add_argument("-m", help="This argument list the directory for the map files",
                        dest="map_dir", type=str, required=True)

    parser.add_argument("-a", help="This argument list the directory for all of the allpair files",
                        dest="allpair_dir", type=str, required=True)

    parser.add_argument("-hap", help="This argument list the directory for all the hapibd files",
                        dest="hapibd_dir", type=str, required=True)

    parser.add_argument("-ilash", help="This argument list the directory for all of the ilash files",
                        dest="ilash_dir", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
