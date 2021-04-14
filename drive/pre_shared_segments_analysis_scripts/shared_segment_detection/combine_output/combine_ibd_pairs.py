import sys
import gzip
import pandas as pd
import itertools
import glob
import os
import re



# Getting all the ibd files that end in .small.txt.gz
import pre_shared_segments_analysis_scripts
from .pair_functions import is_max_pairs_found, after_max_pair_found, Pair_Info_Class
from .build_analysis_dict import get_analysis_files
import utility_scripts


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


def find_ibd_file(ibd_file_list: list, chr_num: str) -> str:
    """Function to return the correct ibd file for the chromosome"""
    
    chr_num = chr_num.strip(".")
    ibd_file: str = [
        file for file in ibd_file_list if "".join(["_", chr_num, "."]) in file
    ][0]

    return ibd_file

def form_file_dict(file_list: list) -> dict:
    """function to form a dictionary of file
    Parameters
    __________
    file_list : list
        list of all the files ibd files"""
    files = {}

    for f in file_list:

        files[f.split(':', 1)[0]] = f.split(':', 1)[1]
    
    return files

def gather_files(ibd_files: dict, segment_dir: str, map_file_dir: str = None) -> dict:
    """Function to gather the necessary files and then return them as a dictionary
    Parameters
    __________
    ibd_files : dict
        dictionary containing two keys: ilash and hapibd. The 
        values of these keys are the directories for each 
        respective program that contains the programs output
    
    segment_dir: str
        string listing the directory that has the .small.txt.gz 
        files
    
    map_file_dir: str
        string that list the directory that has all of the .map 
        files
    
    Returns
    _______
    dict
        returns a dictionary that contains list to all of the output files"""

    # pulling out the hapibd and the ilash directory files               
    ilash_dir_str: str = ibd_files["ilash"]
    hapibd_dir_str: str = ibd_files["hapibd"]

    # using the above two directories to extract the hapibd and iLash files:
    ilash_file_list: list = utility_scripts.get_file_list(
        ilash_dir_str, "*match.gz")
    hapibd_file_list: list = utility_scripts.get_file_list(
        hapibd_dir_str, "*ibd.gz")

    # getting a list of the ibd files

    ibd_file_list: list = utility_scripts.get_file_list(
    segment_dir, "*.small.txt.gz")

    # getting all the map files from the
    if map_file_dir:
        map_file_list: list = utility_scripts.get_file_list(map_file_dir, "*.map")

        return {
        "ilash_file_list": ilash_file_list,
        "hapibd_file_list": hapibd_file_list,
        "map_file_list": map_file_list,
        "ibd_pair_file_list": ibd_file_list
        }

    else:
        return {
            "ilash_file_list": ilash_file_list,
            "hapibd_file_list": hapibd_file_list,
            "ibd_pair_file_list": ibd_file_list
        }
    # output: str = "".join([output, ""])
    
    

def read_first_line(files: list, openfile: dict, endline: dict, curr_pos: dict, curr_ibd: dict, curr_pair: dict, newpos: dict, newline: dict, endtest: dict) -> str:
    """This function will read the first line of the files into the appropriate dictionary
    Parameters
    __________"""

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
    
    return line1

def form_all_combinations(file_dict: dict, all_comb_dict: dict):
    """Function to form all the combinations from the file dict
    Parameters
    __________
    file_dict : dict
        dictionary that has the files per certain chromosomes and the identifier
    
    all_comb_dict : dict
        dictionary that will have the different combinations of items from
        the file dictionary
    """

    for i in range(len(file_dict.keys()), 0, -1):
        for item in list(itertools.combinations(file_dict.keys(), i)):
            all_comb_dict['+'.join(item)] = item


def write_to_file(output_path: str, pair_list: list):
    """Function to write the pair string to a file 
    Parameters
    _________
    output_path : str
        string that has the filepath to write the output to
    
    pair_list : list
        list that has the string of information for each pair in the 
        pair iid list
    """

    with open(output_path, "a+") as allpair_new_file:
            # Checks to see if the file is empty. If it is empty then it write the header
        if os.path.getsize(output_path) == 0:
            # creating a header_line
            allpair_new_file.write(
                "IBD_programs\tpair_1\tpair_2\tchr\tvariant_id\tgene_name\tcarrier_status\tpotential_missed_carrier\tconnected_carriers\thapibd_phase1\thapibd_phase2\tilash_phase1\tilash_phase2\thapibd_start\thapibd_end\thapibd_len\tilash_start\tilash_end\tilash_len\n"
            )
        # writing the pairs strings to the output file
        for pair_str in pair_list:
            allpair_new_file.write(pair_str)

def fix_chr_num(chr_num: str) -> str:
    """Function to convert the chromosome number into a single digit if it is a single digit chromosome so chr08 == chr8
    Parameters
    __________
    chr_num : str
        string that has the two digit chromosome number
    
    Returns
    _______
    str
        returns a string that has the chromosome number formatted for comparison 
        to the ibd files
    """

    if int(chr_num[-2:]) >= 10:

        return chr_num
    
    else:

        return "".join([chr_num[:3], chr_num[-1]])


def combine_output(gathered_file_dict: dict, file_dict: dict, output: str, analysis_type: str, car_file_dir: str= None, pheno_carrier_df: pd.DataFrame=None, pheno_gmap_df: pd.DataFrame=None):
    
    
    for chr_num, identifier in file_dict.keys():

        # Setting a max_number of pairs parameter ot use for comparision so that it only keeps one line
        max_pairs: int = 0

        file_list: list = file_dict[(chr_num, identifier)]

        if len(file_list) == 0:              
            sys.exit(f"no files found for {identifier}")

        # alt_chr_num: str = alternate_chr_num_format(chr_num)
        # getting the hapibd file that corresponds to the correct chromosome number
        # and loading it into a dataframe
        hapibd_file: str = find_ibd_file(gathered_file_dict["hapibd_file_list"], fix_chr_num(chr_num))
        hapibd_df: pd.DataFrame = pd.read_csv(hapibd_file, sep="\t", header=None)
        # getting the ilash file that corresponds to the correct chromsome number and loading it into a dataframe
        ilash_file: str = find_ibd_file(gathered_file_dict["ilash_file_list"], fix_chr_num(chr_num))

        ilash_df: pd.DataFrame = pd.read_csv(ilash_file, sep="\t", header=None)

        
        # The next if/else statement will get the analysis type 
        # dictionary depending on the analysis type
        #If the map_file_list is a key in teh gathered_file_dict then
        # the analysis type is not phenotype and the if statement will return the
        # variant position along with the analysis type
        # If the analysis is the phenotype then it will get the gene 
        # positions from the pheno_gmap_df dictionary. 
        if "map_file_list" in gathered_file_dict:
            map_file: str = [
                file for file in gathered_file_dict["map_file_list"] if "".join([chr_num[:-1], "_"])
            ][0]
            analysis_type_dict: dict = get_analysis_files(analysis_type, identifier, map_file)
        else:
            analysis_type_dict: dict = get_analysis_files(analysis_type, identifier, pheno_gmap_df=pheno_gmap_df)


        output_dir: str = utility_scripts.check_dir(output, "pairs")
        # creating the output path for the file
        out = os.path.join(output_dir, "".join(["IBD_", identifier, ".",chr_num]))
        # next line writes the file
        files = form_file_dict(file_list)

        # creating dictionaries that will be used in the rest of the 
        # code
        curr_pair = {}
        curr_pos = {}
        curr_ibd = {}
        newpos = {}
        newline = {}
        openfile = {}
        endtest = {}
        endline = {}

        # read first line for all input files
        line1: str = read_first_line(files, openfile, endline, curr_pos, curr_ibd, curr_pair, newpos, newline, endtest)

        #TODO: keep working on changing this into smaller functions
        allcomb = {}

        # creating all the combinations based on the file dict
        form_all_combinations(files, allcomb)
    

        combtab = pd.DataFrame(0, columns=files.keys(), index=allcomb.keys())

        for item in allcomb:
            combtab.loc[item, allcomb[item]] = len(allcomb[item])

        oldallpair = set([])


        count: int = 0
        
        while sum(list(map(lambda f: endtest[f],
                           endtest.keys()))) < len(endtest):
            pos = min(newpos.values())
            nowf = findkey(pos, newpos)

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

                    # deleting the file if it is there from a previous run
                    utility_scripts.check_file(allagree_path)

                    # Entering into the get_max_pairs function
                    # TODO: Make a pairs object that can contain the information about the pair object such as the string of pairs, the identifier which is the variant_id or gene name, the chromosome number, the output_path, and the analysis type
                    pair_info_object = Pair_Info_Class(
                        max_pairs_str, identifier, chr_num, allagree_path, analysis_type)

                    # forming a list of iids that that are identified as carrying 
                    # the variant. This list becomes an attribute of the 
                    # pair_info_object called self.iid_list
                    # will check if the pheno_carrier_df exists and if it doesn't
                    # then it will use the car_file_dir

                    if analysis_type == "phenotype":
                        pair_info_object.iid_list_handler(pheno_carriers=pheno_carrier_df)
                    else:
                        
                        pair_info_object.iid_list_handler(carrier_dir=car_file_dir, pheno_carriers=None)
                    #
                    # Next line will actually generate a string with all the necesary information in it
                    pair_info_list: list = pair_info_object.generate_pairs_dict(hapibd_df, ilash_df, analysis_type_dict)

                    write_to_file(pair_info_object.output_path, pair_info_list)

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
