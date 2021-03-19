import sys
import gzip
import pandas as pd
import itertools
import glob
import os
import re

# Getting all the ibd files that end in .small.txt.gz
import pre_shared_segments_analysis_scripts
import utility_scripts
import pair_info


def get_variant_id(ibd_file: str) -> str:
    '''This function will get the proper variant_id'''

    dot_indx_list: list = [
        indx for indx, letter in enumerate(ibd_file) if letter == "."
    ]

    variant_id_handler = {
        4: ibd_file[:dot_indx_list[0]],
        5: ibd_file[:dot_indx_list[1]]
    }

    variant_id: str = variant_id_handler[len(dot_indx_list)]

    return variant_id


def build_file_dict(ibd_file_list: list, program_list: list) -> dict:
    '''This function aligns all the files from each ibd program together'''
    file_dict = dict()

    # iterate through the files to build the dictionary
    for ibd_file in ibd_file_list:

        match = re.search(r'.chr\d\d\.', ibd_file)

        # find chromosome number
        if match:

            chr_num = match.group(0)

        else:
            match = re.search(r'.chr\d_', ibd_file)

            if not match:
                match = re.search(r'.chr\d\.', ibd_file)

            chr_num = match.group(0)

        # Finding the variant id of the file. file names are built so that the
        # variant id sits between the first "_" and the first "."
        for ibd_program in program_list:

            if ibd_program in ibd_file:
                underscore_indx = ibd_file.find("".join([ibd_program, "_"]))

                shorten_ibd_file_string: str = ibd_file[underscore_indx +
                                                        len(ibd_program) + 1:]

        variant_id = get_variant_id(shorten_ibd_file_string)

        # using list comprehension to get all the files that contain that
        # variant and chromosome
        filter_ibd_file_list = [
            file for file in ibd_file_list
            if variant_id in file and chr_num in file
        ]

        # match up the variants with the IBD program
        for ibd_program in program_list:

            # This goes through the three files in the filter_ibd_file_list
            for file in filter_ibd_file_list:

                # This checks to see if the variant id is in the dictionary
                # and that the ibd program is in the file
                if ibd_program in file and (
                        chr_num, variant_id) not in file_dict.keys():

                    # If it is not then the
                    file_dict[(chr_num, variant_id)] = set()

                    file_dict[(chr_num, variant_id)].add("".join(
                        [ibd_program, ":", file]))

                elif ibd_program in file and (chr_num,
                                              variant_id) in file_dict.keys():

                    file_dict[(chr_num, variant_id)].add("".join(
                        [ibd_program, ":", file]))

    return file_dict


# Checking the system arguments


def findkey(i, mydict):
    result = []
    offset = -1
    while True:
        try:
            offset = list(mydict.values()).index(i, offset + 1)
        except ValueError:
            return list(map(lambda i: list(mydict.keys())[i], result))
        result.append(offset)


def allinter(mylist,
             curr_pair) -> int:  # Finds the pair that intersects for all files
    intu = curr_pair[mylist[0]]
    for f in mylist[1:]:
        intu = intu & curr_pair[f]
    return intu


def get_uniqrow(i, allcomb, curr_pair, combtab) -> list:
    uniqdic = {}
    for comb in allcomb.keys():
        raw_n = len(allinter(allcomb[comb], curr_pair))
        octab = pd.DataFrame(
            map(lambda ff: combtab[ff] > combtab.loc[comb, ff],
                allcomb[comb])).all()
        overcount = octab.index[octab == True].tolist()
        uniqdic[comb] = raw_n - \
            sum(list(map(lambda oc: uniqdic[oc], overcount)))
    return list(uniqdic.values())


# This is pair that has union for all files
def all_agree_pair(pair_list: dict) -> list:
    unionpair = list(pair_list.values())[0]
    for f in pair_list.keys():
        unionpair = unionpair.union(pair_list[f])
    return unionpair


def get_pair_info(pair_info_object, input_file_object) -> dict:
    """Function to get a dictionary of The pairs and other pair information
    Parameters
    __________

    Returns
    _______
    dict
        returns a dictionary where the key is the tuple of (pair1, pair2) and there 
        are values "carrier_status", "connected_carriers", "missed_carriers"
    """

    pair_info_dict: dict = form_pair_dict(allpair_file_name, pair_list, variant_id, carrier_file_dir,
                                          chr_num, hapibd_file, ilash_file, map_file)

    return pair_info_dict


def get_carrier_list(file: str, variant_id: str) -> list:
    """Function to get a list of IIDs who are carrier from the carrier file
    Parameters
    __________
    file : str
        string that list the path to a file called chr**.single_variant_carrier.csv. 

    variant_id : str
        string containing the illumina mega probe id for the variant of
        interest

    Returns
    _______
    list
        returns a list of IIDs who carry the variant of interest
    """

    carrier_df = pd.read_csv(file, sep=",")

    # filter dataframe for those individuals carrying the variant
    carrier_list = carrier_df[carrier_df["Variant ID"] ==
                              variant_id]["IID"].values.tolist()
    return carrier_list


# def alternate_chr_num_format(chr_num: str) -> str:
#     '''This function will add a zero between the number if it is a single digit chromosome'''

#     if len(chr_num.strip(".")) == 4:

#         chr_num = chr_num.strip(".")

#         chr_num = "".join([".", chr_num[:-1], "0", chr_num[3], "."])

#     return chr_num


def check_for_missed_carriers(pair_2: str, pair_list: list,
                              carrier_list: list) -> int:
    '''This function will check for variants that may be missed carriers'''

    # get a list of all strings that contain the pair
    pair2_in_str: list = [string for string in pair_list if pair_2 in string]

    connected_carriers_list: list = []

    for pair in pair2_in_str:

        pair_iids: list = pair.split(":")[1]

        pair_1, pair_2 = pair_iids.split("-")

        if pair_1 in carrier_list or pair_2 in carrier_list:

            connected_carriers_list.append(pair)

    return int(len(connected_carriers_list))


def form_pair_dict(pair_info_object, input_file_object,write_path: str, pair_list: list, variant_id: str,
                   carrier_file_dir: str, chr_num: str, hapibd_file: str,
                   ilash_file: str, map_file: str) -> dict:
    '''this function takes the original allpair.txt file and reformats it to four columns.
    The first columnn is the ibd program that identified the pairs, the second column is the
    first pair, the third column is the second pair, and then the fourth column tells if the
    second pair is a carrier'''

    # Removing the allpair.txt file if it already exists from a previous run
    if os.path.isfile(pair_info_object.output_path):

        # removing the file
        os.remove(pair_info_object.output_path`)


    # use list comprehension to find the file with that chr_num
    carrier_list: list = input_file_object.get_carrier_file_list()
    # alt_chr_num: str = alternate_chr_num_format(chr_num)
    carrier_file = [
        file for file in carrier_list
        if chr_num.strip(".") in file][0]

    # getting the list of carriers' iids for the specific variant

    carrier_iid_list = get_carrier_list(carrier_file, variant_id)
    pairs_dict: dict = {}
    for pair in pair_list:

        # with open(write_path, "a+") as allpair_new_file:
        #     # Checks to see if the file is empty. If it is empty then it write the header
        #     if os.path.getsize(write_path) == 0:

        #         allpair_new_file.write(
        #             "IBD_programs\tpair_1\tpair_2\tchr\tvariant_id\tcarrier_status\tpotential_missed_carrier\tconnected_carriers\thapibd_phase1\thapibd_phase2\tilash_phase1\tilash_phase2\thapibd_start\thapibd_end\thapibd_len\tilash_start\tilash_end\tilash_len\n"
        #         )

        # getting the IBD programs that found the output
        programs = pair.split(":")[0]

        # getting pair 1
        pair1 = pair.split(":")[1].split("-")[0]

        # Removing the extra newline off of the final pair
        pair1 = pair1.strip("\n")

        # getting pair 2
        pair2 = pair.split(":")[1].split("-")[1]

        # stripping the extra newline off of the final pair
        pair2 = pair2.strip("\n")

        pairs_dict[(pair1, pair2)] = {}

        # determining carrier status of second iid pair
        if pair2 in carrier_iid_list:
            # the carrier status is a one if the second iid pair is in the list of carriers

            pairs_dict[(pair1, pair2)]["carrier_status"] = 1
            print("car_status = 1")
            print(pair1, pair2)

            # hapibd_segment_dict, ilash_segment_dict = pre_shared_segments_analysis_scripts.get_segment_lens(
            #     pair1, pair2, map_file, variant_id, hapibd_file,
            #     ilash_file)
            # print(hapibd_segment_dict)

        else:
            # If the second iid pair is not in the list of carriers then the status is a 0
            pairs_dict[(pair1, pair2)]["carrier_status"] = 0
        # This function will check if the second pair is a potential missed carrierq

            print("car_status = 0")
            print(pair1, pair2)

            # hapibd_segment_dict, ilash_segment_dict = pre_shared_segments_analysis_scripts.get_segment_lens(
            #     pair1, pair2, map_file, variant_id, hapibd_file,
            #     ilash_file)

            # print("hapibd_segment_dict")
            # print(hapibd_segment_dict)

            connected_carriers = check_for_missed_carriers(
                pair2, pair_list, carrier_iid_list)

            pairs_dict[(pair1, pair2)
                       ]["connected_carriers"] = connected_carriers

            if connected_carriers >= 2:

                pairs_dict[(pair1, pair2)]["missed_carrier"] = 1

            else:
                pairs_dict[(pair1, pair2)]["missed_carrier"] = 0

        if pairs_dict[(pair1, pair2)]["carrier_status"] == 1:

            connected_carriers = "N/A"

            pairs_dict[(pair1, pair2)
                       ]["connected_carriers"] = connected_carriers

            pairs_dict[(pair1, pair2)]["missed_carrier"] = 0

            print(pair1, pair2)
            # hapibd_segment_dict, ilash_segment_dict = pre_shared_segments_analysis_scripts.get_segment_lens(
            #     pair1, pair2, map_file, variant_id, hapibd_file,
            #     ilash_file)
            # print(hapibd_segment_dict)

            # if not hapibd_segment_dict:
            #     hapibd_segment_dict = {
            #         "start": "N/A",
            #         "end": "N/A",
            #         "length": "N/A",
            #         "phase1": "N/A",
            #         "phase2": "N/A"
            #     }
            # if not ilash_segment_dict:
            #     ilash_segment_dict = {
            #         "start": "N/A",
            #         "end": "N/A",
            #         "length": "N/A",
            #         "phase1": "N/A",
            #         "phase2": "N/A"
            #     }

            # # Need to call the haplotype function to get the lengths
            # # writing the pairs to different columns
            # allpair_new_file.write(
            #     f"{programs}\t{pair1}\t{pair2}\t{chr_num.strip('.')}\t{variant_id}\t{car_status}\t{str(potential_missed_carrier)}\t{str(connected_carriers)}\t{hapibd_segment_dict['phase1']}\t{hapibd_segment_dict['phase2']}\t{ilash_segment_dict['phase1']}\t{ilash_segment_dict['phase2']}\t{hapibd_segment_dict['start']}\t{hapibd_segment_dict['end']}\t{hapibd_segment_dict['length']}\t{ilash_segment_dict['start']}\t{ilash_segment_dict['end']}\t{ilash_segment_dict['length']}\n"
            # )
    return pairs_dict


def is_max_pairs_found(curr_max_pairs: int, new_max_pairs: int) -> int:

    pair_handler_dict = {False: 0, True: 1}

    # This function will return either 1 or zero from the pair handler_dict based on whether or not curr_max_pair from the previous row is less than or greater than the max pairs from the current row
    return pair_handler_dict[curr_max_pairs >= new_max_pairs]


def after_max_pair_found(curr_max_pair: int, new_max_pair: int) -> int:
    '''This function will break out of the loop if a max pair is found and then the size of the pair list starts decreasing'''

    # If the previous rows max number of pairs is higher than the current row then this function will return a one
    if curr_max_pair > new_max_pair:
        return 1
    elif curr_max_pair == new_max_pair:  # If the the above is not true than it returns a 0
        return 0


# TODO: incorporate this function into the main combine_output function
def get_file(file_list: list, identifier: str = None, chr_num=None) -> str:
    '''This function gets the file that matches a condition from a list of files'''

    # generate alternate chr number incase the formatting does not contain a zero
    if identifier:
        file_str: str = [file for file in file_list if identifier in file][0]

    # alt_chr_num = None
    if chr_num:
        # alt_chr_num: str = alternate_chr_num_format(chr_num)
        file_str: str = [
            file for file in file_list
            if chr_num in file
        ][0]

    return file_str


def find_ibd_file(ibd_dir: list, chr_num: str) -> str:
    """Function to return the correct ibd file for the chromosome"""

    ibd_file: str = [
        file for file in ibd_dir if "".join(["_", chr_num[1:]]) in file
    ][0]

    return ibd_file


def combine_output(segment_dir: str, ibd_programs: list, output: str,
                   car_file: str, ibd_files: dict, map_file_dir: str):
    ilash_dir_str: str = ibd_files["ilash"]
    hapibd_dir_str: str = ibd_files["hapibd"]
    # using the above two directories to extract the hapibd and iLash files:
    ilash_file_list: list = utility_scripts.get_file_list(
        ilash_dir_str, "*match.gz")
    hapibd_file_list: list = utility_scripts.get_file_list(
        hapibd_dir_str, "*ibd.gz")

    # getting all the map files from the
    map_file_list: list = utility_scripts.get_file_list(map_file_dir, "*.map")

    output: str = "".join([output, ""])

    ibd_file_list: list = utility_scripts.get_file_list(
        segment_dir, "*.small.txt.gz")

    file_dict: dict = build_file_dict(ibd_file_list, ibd_programs)

    for chr_num, variant_id in file_dict.keys():
        print(chr_num)
        # Setting a max_number of pairs parameter ot use for comparision so that it only keeps one line
        max_pairs: int = 0

        file_list: list = file_dict[(chr_num, variant_id)]

        # alt_chr_num: str = alternate_chr_num_format(chr_num)
        # getting the hapibd file that corresponds to the correct chromosome number
        hapibd_file: str = find_ibd_file(hapibd_file_list, chr_num)
        # getting the ilash file that corresponds to the correct chromsome number
        ilash_file: str = find_ibd_file(ilash_file_list, chr_num)

        # getting the map file that corresponds to the correct chromosome number
        map_file: str = [
            file for file in map_file_list if "".join([chr_num[:-1], "_"])
        ][0]

        if len(file_list) == 0:  # Checking length of system arguments
            sys.exit(f"no files found for this variant {variant_id}")

        # Making the first argument the output variable

        out = "".join([output, "IBD_", variant_id, chr_num[:-1]])
        # next three lines write the files to a dictionary
        files = {}

        for f in file_list:

            files[f.split(':', 1)[0]] = f.split(':', 1)[1]

        # print('input {0} files: {1}'.format(len(files),
        #                                  ' '.join(files.keys())))

        curr_pair = {}
        curr_pos = {}
        curr_ibd = {}
        newpos = {}
        newline = {}
        openfile = {}
        endtest = {}
        endline = {}

        # read first line for all input files
        for f in files:
            openfile[f] = gzip.open(files[f], 'rt')
            openfile[f].seek(0, 2)
            endline[f] = openfile[f].tell()
            openfile[f].seek(0)
            line0 = openfile[f].readline()
            curr_pos[f] = 0
            curr_ibd[f] = set([])
            curr_pair[f] = set([])
            line1 = openfile[f].readline()
            line1 = line1.strip()
            newpos[f] = int(line1.split('\t')[1])
            newline[f] = line1
            endtest[f] = 0

        allcomb = {}
        for i in range(len(files.keys()), 0, -1):
            for item in list(itertools.combinations(files.keys(), i)):
                allcomb['+'.join(item)] = item

        combtab = pd.DataFrame(0, columns=files.keys(), index=allcomb.keys())

        for item in allcomb:
            combtab.loc[item, allcomb[item]] = len(allcomb[item])

        oldallpair = set([])

        # this sets a counter to determine how many times the max_pair_int branch is entered
        count: int = 0

        while sum(list(map(lambda f: endtest[f],
                           endtest.keys()))) < len(endtest):
            pos = min(newpos.values())
            nowf = findkey(pos, newpos)

            #    print('{0} from {1}'.format(str(pos), ' '.join(nowf)))
            for f in nowf:

                CHR = str(newline[f].split('\t')[0])
                if newline[f].split('\t')[4] != 'NA':
                    addibd = set(newline[f].split('\t')[4].split(' '))
                else:
                    addibd = set([])
                if newline[f].split('\t')[5] != 'NA':
                    delibd = set(newline[f].split('\t')[5].split(' '))
                else:
                    delibd = set([])
                curr_ibd[f] = set(set(curr_ibd[f] | addibd) - delibd)
                if len(curr_ibd[f]) > 0:
                    curr_pair[f] = set(
                        map(lambda x: x.split(':')[1], curr_ibd[f]))
                else:
                    curr_pair[f] = set([])
                nextline = openfile[f].readline()
                nextline = nextline.strip()
                if nextline == '' and openfile[f].tell() == endline[f]:
                    endtest[f] = 1
                    newpos[f] = float('inf')
                else:
                    #            nextline = nextline.strip()
                    newpos[f] = int(nextline.split('\t')[1])
                    newline[f] = nextline

            # print('chr{0}:{1} from {2}'.format(CHR, str(pos), ' '.join(nowf)))

            # sumrow: list = list(
            #     map(lambda comb: len(allinter(comb, curr_pair)),
            #         allcomb.values()))
            # uniqrow = get_uniqrow(1, allcomb, curr_pair, combtab)

            newallpair: list = all_agree_pair(curr_pair)

            max_pairs_int: int = is_max_pairs_found(max_pairs, len(newallpair))

            outpair: list = []

            for pp in newallpair:
                tool = []
                for f in curr_pair.keys():
                    if pp in curr_pair[f]:
                        tool.append(f)
                outpair.append('{0}:{1}'.format(','.join(tool), pp))

            if len(outpair) == 0:
                outpair = ['NA']

            if max_pairs_int == 1:

                # update the counter
                count += 1

                # This will return a 0 if the max pairs == the len(newallpair) and a one if max_pairs is greater
                after_max_pair: int = after_max_pair_found(
                    max_pairs, len(newallpair))

                if after_max_pair == 0 and count == 1:

                    max_pairs_str: str = previous_row_str

                    # get the start base position which will be used later in the allpair path
                    start_bp: str = previous_row_bp

                end_bp: str = str(pos)

                if after_max_pair == 1:
                    # creating an output path that just has the allpair.txt files and

                    allagree_path = "".join(
                        [out, ".", start_bp, "-", end_bp, ".allpair.txt"])

                    # Entering into the get_max_pairs function
                    # TODO: Make a pairs object that can contain some of this information
                    pair_info_object = pair_info.Pair_Info_Class(
                        max_pairs_str, variant_id, chr_num, allagree_path)

                    input_file_object = pair_info.Input_Files(
                        car_file, hapibd_file, ilash_file, map_file)

                    get_pair_info(pair_info_object, input_file_object)

                    break

            if max_pairs_int == 0:
                # This if statement is made to reset the counter if the max_pair integer goes from being 1 back to 0

                # Reseting the counter
                count = 0

            max_pairs = len(newallpair)
            # keeping track of the previous row so that it can be used if necessary
            previous_row_str: str = f"{str(CHR)}\t{str(pos)}\tNA\t{len(newallpair)}\t{' '.join(outpair)}\n"

            # Also keeping track of the base position
            previous_row_bp: str = str(pos)
