# This script will get the haplotype information
import pandas as pd
import argparse
import os
from os import path
import glob
import re
import gzip
from functools import partial
import multiprocessing as mp

from haplotype_segments_analysis.network_ids import filter_df, filter_for_pairs
import utility_scripts


def identify_unique_variants(confirmed_carrier_file: str) -> list:
    '''This function will get all the unique variants that we have a confirmed carrier for'''

    # loads the variant_id column from the file into dataframe
    confirmed_carrier_df: pd.DataFrame = pd.read_csv(confirmed_carrier_file,
                                                     sep="\t",
                                                     usecols=["variant_id"])

    # gets a list of all the unique variants in the file
    unique_var_list: list = confirmed_carrier_df.variant_id.unique().tolist()
    # print(f"{len(unique_var_list)} unique variants identified")

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

    match = re.search(r'.chr\d\d\.', allpair_file)

    # find chromosome number
    if match:

        chr_num: str = match.group(0)

    else:
        match = re.search(r'_chr\d.', allpair_file)

        chr_num: str = match.group(0)

    chr_num = chr_num[1:-1]

    return chr_num


def alternate_chr_num_format(chr_num: str) -> str:
    '''This removes the string in some of the chromosome files'''

    zero_indx: int = chr_num.find("0")

    if zero_indx and zero_indx == 4:

        chr_list: list = chr_num.split("0")

        chr_num = "".join(chr_list)

    return chr_num


def get_file(file_list: list, identifier: str = None, chr_num=None) -> str:
    '''This function gets the file that matches a condition from a list of files'''

    # generate alternate chr number incase the formatting does not contain a zero
    if identifier:
        file_str: str = [file for file in file_list if identifier in file][0]

    alt_chr_num = None
    if chr_num:
        alt_chr_num: str = alternate_chr_num_format(chr_num)
        file_str: str = [
            file for file in file_list
            if chr_num in file or alt_chr_num in file
        ][0]

    return file_str


def get_var_pos(map_file: str, variant_id: str) -> str:
    '''This function find the position of the variant id and return that as a string'''

    # This reads in the columns from the map file that contains the variant probe id and the pos
    map_df: pd.DataFrame = pd.read_csv(map_file,
                                       sep="\t",
                                       header=None,
                                       usecols=[1, 3])

    # getting the value
    pos: str = str(map_df[map_df[1] == variant_id][3].values[0])

    return pos


def check_ibd_program(file_str: str) -> int:
    '''This function will return a zero or a 1 depending on whether or not the file is a hapibd file or an ilash file respectively'''

    # create a dictionary to handle the situation
    ibd_handler: dict = {True: 1, False: 0}

    # Checks to see if the file is a ilash file
    handler_int: int = ibd_handler["match" in file_str]

    return handler_int


def get_index_positions(ibd_program_int: int) -> list:

    # This dictionary contains the index of values such as IID1, IID2, Start point
    # end point, and the shared segment length. True will contain indixes for an iLash
    # file and false will contain the indices for a hapibd file
    indx_dict: dict = {1: [0, 2, 5, 6, 9], 0: [0, 2, 5, 6, 7]}

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

    start_indx: int = indx_list[2]
    end_indx: int = indx_list[3]

    # figure out if file is ilash or hapibd
    with gzip.open(file, "rt") as ibd_file:

        for row in ibd_file:

            split_row: list = row.split()

            # This line checks to see if the variant position falls inbetween the start and end position of the row
            if int(split_row[start_indx]) <= int(variant_pos) and int(
                    split_row[end_indx]) >= int(variant_pos):
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
    carrier_list: list = carrier_df[carrier_df["Variant ID"] ==
                                    variant_id]["IID"].values.tolist()

    return carrier_list


def filter_for_carriers(ibd_pair_list: list, carrier_on_mega_list: list,
                        confirmed_carrier_list: list) -> list:
    '''This function will filter the pair list for only pairs where both individual are carriers'''

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

        if iid1 in carrier_on_mega_list and iid2 in carrier_on_mega_list:

            if iid1 in confirmed_carrier_list and iid2 in confirmed_carrier_list:
                filtered_list.append(pair_str)

    return filtered_list


def get_confirmed_carriers(confirmed_file: str, var_id: str) -> list:
    '''This function will generate a list of carriers who are confirmed through shared segment'''

    confirmed_carriers_df: pd.DataFrame = pd.read_csv(confirmed_file, sep="\t")

    filter_for_var_df: pd.DataFrame = confirmed_carriers_df[
        confirmed_carriers_df.variant_id == var_id]

    filtered_for_confirmed_status_df: pd.DataFrame = filter_for_var_df[
        filter_for_var_df.confirmed_status == 1]

    confirmed_carrier_list: list = filtered_for_confirmed_status_df.IID.values.tolist(
    )

    return confirmed_carrier_list


def create_output_str(hapibd_filtered_list: str, ilash_filtered_list: str,
                      chr_num: str, variant_id: str, que_object,
                      var_que_object, network_file_path: str):
    '''This function will create a string containing the start and end points for each shared segment between two pairs.
    The function will place this string in the que of a listener that will then write it to a file'''

    # get the chromosome number digit to be used later
    chr_digit: str = re.findall(r'[0-9]+', chr_num)[0]

    # creatinga an empty list that will contain the output from the files
    # creating a row of N/A's incase the output string should fail

    output_str: str = f"N/A\tN/A\t{chr_digit}\t{variant_id}\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\n"

    # getting the matched values
    matched_list: list = [
        (string, pair) for string in hapibd_filtered_list
        for pair in ilash_filtered_list
        if string.split("\t")[0] in pair and string.split("\t")[2] in pair
    ]

    # getting the values that are unique in the hapibd_filtered_list

    # creating a list with only the hapibd values that are matched
    hapibd_matched: list = [pair[0] for pair in matched_list]

    hapibd_only_list: list = [
        string for string in hapibd_filtered_list
        if string not in hapibd_matched
    ]

    # getting all the values that are only in the ilash set

    ilash_matched: list = [pair[1] for pair in matched_list]

    ilash_only_list: list = [
        string for string in ilash_filtered_list if string not in ilash_matched
    ]

    # Combining the previous list to make a total of three list
    total_pair_list: list = matched_list + hapibd_only_list + ilash_only_list

    # Accounting for the situation where there are no pairs
    if len(total_pair_list) == 0:
        failed_var_str: str = f"{variant_id}\t{chr_num}\n"

        var_que_object.put(failed_var_str)

    # If the length is 2 then this means there is info from hapibd and ilash
    for pair_str in total_pair_list:

        if len(pair_str) == 2:

            split_hapibd_list: list = pair_str[0].split("\t")
            split_ilash_list: list = pair_str[1].split("\t")

            pair_1: str = split_hapibd_list[0]
            pair_2: str = split_hapibd_list[2]

            # getting the network ids for this set of pairs. Have to remove final two characters from teh variant id name
            filtered_network_df: pd.DataFrame = filter_df(
                network_file_path, variant_id)

            network_id: str = filter_for_pairs(filtered_network_df, pair_1,
                                               pair_2)

            output_str: str = f"{split_hapibd_list[0]}\t{split_hapibd_list[2]}\t{chr_digit}\t{variant_id}\t{network_id}\t{split_hapibd_list[5]}\t{split_hapibd_list[6]}\t{split_hapibd_list[7]}\t{split_ilash_list[5]}\t{split_ilash_list[6]}\t{split_ilash_list[9]}\n"

        else:

            split_str_list: list = pair_str.split("\t")

            pair_1: str = split_str_list[0]
            pair_2: str = split_str_list[2]

            # getting the network ids for this set of pairs. Have to remove final two characters from teh variant id name
            filtered_network_df: pd.DataFrame = filter_df(
                network_file_path, variant_id)

            network_id: str = filter_for_pairs(filtered_network_df, pair_1,
                                               pair_2)

            if len(split_str_list) == 11:

                output_str: str = f"{split_str_list[0]}\t{split_str_list[2]}\t{chr_digit}\t{variant_id}\t{network_id}\tN/A\tN/A\tN/A\t{split_str_list[5]}\t{split_str_list[6]}\t{split_str_list[9]}\n"

            elif len(split_str_list) == 8:

                output_str: str = f"{split_str_list[0]}\t{split_str_list[2]}\t{chr_digit}\t{variant_id}\t{network_id}\t{split_str_list[5]}\t{split_str_list[6]}\t{split_str_list[7]}\tN/A\tN/A\tN/A\n"

        que_object.put(output_str)


def get_full_var_name(allpair_file: str, variant_id: str) -> str:
    '''This function will get the full variant name out of the allpair file name'''

    variant_indx: int = allpair_file.index(variant_id)

    # This is the ending index of the slice
    second_indx: int = variant_indx + len(variant_id) + 2

    full_variant_id: str = allpair_file[variant_indx:second_indx]

    return full_variant_id


def get_haplotype(allpair_file_list: list, carrier_file_list: list,
                  map_file_list: list, ilash_file_list: list,
                  hapibd_file_list: list, que_object, var_que_object,
                  network_file_path: str, confirmed_carrier_file: str,
                  variant: str):
    '''This function will contain the main segments of code that will run in the run function.
    It will be used in the parallel_map funcion which is an attempt to parallelize the function.'''

    allpair_file: str = get_file(allpair_file_list, identifier=variant)

    # print(f"using allpair file {allpair_file}")

    chr_num: str = get_chr_id(allpair_file)

    full_variant_id: str = get_full_var_name(allpair_file, variant)

    map_file: str = get_file(map_file_list,
                             chr_num="".join([".", chr_num, "_"]))

    # getting the specific carrier file for the chromosome
    carrier_file: str = get_file(carrier_file_list,
                                 chr_num="".join([chr_num, "."]))

    carriers_list: list = get_carriers(carrier_file, full_variant_id)

    confirmed_carriers_list: list = get_confirmed_carriers(
        confirmed_carrier_file, variant)

    var_pos: str = get_var_pos(map_file, variant)

    # The next four lines get the hapibd file and the ilash file for the specific chromosome

    ilash_file: str = get_file(ilash_file_list,
                               chr_num="".join(["_", chr_num, "."]))

    hapibd_file: str = get_file(hapibd_file_list,
                                chr_num="".join(["_", chr_num, "."]))

    # getting a list of each pair
    ilash_list: list = filter_file(ilash_file, var_pos)

    ilash_carriers_list: list = filter_for_carriers(ilash_list, carriers_list,
                                                    confirmed_carriers_list)

    hapibd_list: list = filter_file(hapibd_file, var_pos)

    hapibd_carriers_list: list = filter_for_carriers(hapibd_list,
                                                     carriers_list,
                                                     confirmed_carriers_list)

    create_output_str(hapibd_carriers_list, ilash_carriers_list, chr_num,
                      full_variant_id, que_object, var_que_object,
                      network_file_path)


def parallel_map(workers: int, allpair_file_list: list, variant_list: list,
                 carrier_file_list: list, map_file_list: list, output: str,
                 ilash_file_list: list, hapibd_file_list: list,
                 network_file_path: str, confirmed_carrier_file: str):
    print("attempting to run in parallel...")

    # creating a manager que that can be used to line up how it write to the file
    manager = mp.Manager()

    # creating the que object for the haplotypes
    que = manager.Queue()

    # create a que object for the failed variants
    variant_que = manager.Queue()

    pool = mp.Pool(workers)

    header: str = f"pair_1\tpair_2\tchr\tvariant_id\tnetwork_id\thapibd_start\thapibd_end\thapibd_len\tilash_start\tilash_end\tilash_len\n"
    # activate the listener function so that it can write from the que as it is going
    watcher = pool.apply_async(
        utility_scripts.listener,
        (que, "".join([output, "haplotype_lengths.txt"]), header))

    variant_header: str = f"variant\tchr\n"
    var_watcher = pool.apply_async(
        utility_scripts.listener, (variant_que, "".join(
            [output, "nopairs_haplotype_analysis.txt"]), variant_header))

    # creating a partial function so that we can pass the necessary parameters to the get_haplotype function
    func = partial(get_haplotype, allpair_file_list, carrier_file_list,
                   map_file_list, ilash_file_list, hapibd_file_list, que,
                   variant_que, network_file_path, confirmed_carrier_file)

    pool.map(func, variant_list)

    que.put("kill")
    variant_que.put("kill")
    pool.close()

    pool.join()


def remove_previous_file(file_path: str):
    '''This function will remove previous output files from previous runs'''
    # This section will check if the output file exist from a previous run and if it does then it will delete it
    if path.exists(file_path):

        os.remove(file_path)


def sort_file(output_file_path: str):
    '''This function will load the haplotype information into a datframe and then it will sort the dataframe and write that to the file'''

    haplotype_df: pd.DataFrame = pd.read_csv(output_file_path, sep="\t")

    haplotype_df = haplotype_df.sort_values("variant_id")

    haplotype_df.to_csv(output_file_path, na_rep="N/A", index=False, sep="\t")


def create_readme(output_path: str):
    '''This function creates a readme for the specified output path'''
    readme = utility_scripts.Readme("_README.md", output_path)
    readme.rm_previous_file()
    readme.write_header("haplotype_analysis/")
    readme.create_date_info()
    readme.add_line(utility_scripts.haplotype_analysis_body_text)


def get_segment_lengths(confirmed_carrier_file: str, output_path: str,
                        ilash_dir: str, hapibd_dir: str, threads: int,
                        map_file_dir: str, reformated_carrier_files: str,
                        network_file: str, allpair_files: str) -> str:
    "function to run"

    create_readme(output_path)
    # removing output files from previous runs
    remove_previous_file("".join([output_path, "haplotype_lengths.txt"]))

    remove_previous_file("".join(
        [output_path, "nopairs_haplotype_analysis.txt"]))

    # getting a list of all the variants that have a confirmed carrier
    variant_list: list = identify_unique_variants(confirmed_carrier_file)
    # Getting all the carrier files

    carrier_file_list: list = get_file_list(reformated_carrier_files, "*.csv")
    # getting the map_files
    map_file_list: list = get_file_list(map_file_dir, "*.map")

    # getting the list of allpair files
    allpair_file_list: list = get_file_list(allpair_files, "*.allpair.txt")

    # getting a list of the ilash and hapibd files
    ilash_file_list: list = get_file_list(ilash_dir, "*.match.gz")

    hapibd_file_list: list = get_file_list(hapibd_dir, "*.ibd.gz")

    file_handler_dict: dict = {
        "allpair_files": allpair_file_list,
        "map_files": map_file_list,
        "carrier_files": carrier_file_list,
        "ilash_files": ilash_file_list,
        "hapibd_files": hapibd_file_list
    }

    parallel_runner: object = utility_scripts.Haplotype_Parallel_Runner(
        int(threads), output_path, file_handler_dict, network_file,
        variant_list, confirmed_carrier_file)

    file_name: str = "haplotype_lengths.txt"

    header_str: str = "pair_1\tpair_2\tchr\tvariant_id\tnetwork_id\thapibd_start\thapibd_end\thapibd_len\tilash_start\tilash_end\tilash_len\n"

    parallel_runner.run_haplotypes_parallel(file_name, get_haplotype,
                                            header_str)

    sort_file("".join([output_path, "haplotype_lengths.txt"]))

    return "".join([output_path, "haplotype_lengths.txt"])
