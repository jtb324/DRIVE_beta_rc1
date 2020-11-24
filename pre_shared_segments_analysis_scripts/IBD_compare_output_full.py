import sys
import gzip
import pandas as pd
import itertools
import glob
import os
import re
import argparse

# Getting all the ibd files that end in .small.txt.gz


def gather_ibd_files(segment_dir: str) -> list:
    '''This function will get all of the files for a specific chromosome'''
    cur_dir = os.getcwd()
    os.chdir(segment_dir)

    ibd_file_list = []

    for file in glob.glob("*.small.txt.gz"):

        full_file_path = "".join([segment_dir, file])

        ibd_file_list.append(full_file_path)

    os.chdir(cur_dir)

    return ibd_file_list


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

                underscore_indx = ibd_file.find(
                    "".join([ibd_program, "_"]))

                ibd_program_indx = underscore_indx+len(ibd_program)+1

        dot_indx = ibd_file.find(".")

        variant_id = ibd_file[ibd_program_indx:dot_indx]

        # There are a few variants that have a dot in the middle and therefore
        # the above code would split it.
        # Therefore if there is no number in the string than it enters this if
        # statement and then finds the next dot.
        if not any(map(str.isdigit, variant_id)):
            # This finds the second dot in the string
            dot_indx = ibd_file.find(".", dot_indx+1)

            # This finds the full variant id
            variant_id = ibd_file[ibd_program_indx:dot_indx]

        # using list comprehension to get all the files that contain that
        # variant and chromosome
        filter_ibd_file_list = [
            file for file in ibd_file_list if variant_id in file and chr_num in file]

        # match up the variants with the IBD program
        for ibd_program in program_list:

            # This goes through the three files in the filter_ibd_file_list
            for file in filter_ibd_file_list:

                # This checks to see if the variant id is in the dictionary
                # and that the ibd program is in the file
                if ibd_program in file and (chr_num, variant_id) not in file_dict.keys():

                    # If it is not then the
                    file_dict[(chr_num, variant_id)] = set()

                    file_dict[(chr_num, variant_id)].add(
                        "".join([ibd_program, ":", file]))

                elif ibd_program in file and (chr_num, variant_id) in file_dict.keys():

                    file_dict[(chr_num, variant_id)].add(
                        "".join([ibd_program, ":", file]))

    return file_dict

# Checking the system arguments


def findkey(i, mydict):
    result = []
    offset = -1
    while True:
        try:
            offset = list(mydict.values()).index(i, offset+1)
        except ValueError:
            return list(map(lambda i: list(mydict.keys())[i], result))
        result.append(offset)


def allinter(mylist, curr_pair) -> int:  # Finds the pair that intersects for all files
    intu = curr_pair[mylist[0]]
    for f in mylist[1:]:
        intu = intu & curr_pair[f]
    return intu


def get_uniqrow(i, allcomb, curr_pair, combtab) -> list:
    uniqdic = {}
    for comb in allcomb.keys():
        raw_n = len(allinter(allcomb[comb], curr_pair))
        octab = pd.DataFrame(
            map(lambda ff: combtab[ff] > combtab.loc[comb, ff], allcomb[comb])).all()
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


def get_max_pairs(allpair_file_name: str, pairs_row: str, variant_id: str, carrier_file_dir: str, chr_num: str):
    '''This function will find the highest number of pairs and writes it to a file called .allpair.txt'''

    # spliting the pair row
    # {str(CHR)}\t{str(pos)}\tNA\t{len(newallpair)}\t{' '.join(outpair)}

    pair_str: str = pairs_row.split("\t")[4]

    pair_list: list = pair_str.split(" ")

    reformat(allpair_file_name,
             pair_list, variant_id, carrier_file_dir, chr_num)


def get_carrier_file_list(carrier_files_dir: str) -> list:
    '''This function will return a list of files that describe the carriers for each variant. It will use the reformated files'''
    cur_dir = os.getcwd()
    os.chdir(carrier_files_dir)

    carrier_file_list = []

    for file in glob.glob("*.csv"):

        full_file_path = "".join([carrier_files_dir, file])

        carrier_file_list.append(full_file_path)

    os.chdir(cur_dir)

    return carrier_file_list


def get_carrier_list(file: str, variant_id: str) -> list:
    '''Getting the list of characters'''
    carrier_df = pd.read_csv(file, sep=",")

    # filter dataframe for those individuals carrying the variant
    carrier_list = carrier_df[carrier_df["Variant ID"]
                              == variant_id]["IID"].values.tolist()

    return carrier_list


def reformat(write_path: str, pair_list: list, variant_id: str, carrier_file_dir: str, chr_num: str):
    '''this function takes the original allpair.txt file and reformats it to four columns.
    The first columnn is the ibd program that identified the pairs, the second column is the 
    first pair, the third column is the second pair, and then the fourth column tells if the 
    second pair is a carrier'''

    # Removing the allpair.txt file if it already exists from a previous run
    if os.path.isfile('filename.txt'):

        # removing the file
        os.remove(write_path)

    carrier_list = get_carrier_file_list(carrier_file_dir)

    # use list comprehension to find the file with that chr_num

    carrier_file = [
        file for file in carrier_list if chr_num.strip(".") in file][0]

    # getting the list of carriers' iids for the specific variant

    carrier_iid_list = get_carrier_list(carrier_file, variant_id)

    for pair in pair_list:

        with open(write_path, "a+") as allpair_new_file:
            # Checks to see if the file is empty. If it is empty then it write the header
            if os.path.getsize(write_path) == 0:

                allpair_new_file.write(
                    "IBD_programs\tpair_1\tpair_2\tcarrier_status\n")

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

            # determining carrier status of second iid pair
            if pair2 in carrier_iid_list:
                # the carrier status is a one if the second iid pair is in the list of carriers
                car_status = 1

            else:
                # If the second iid pair is not in the list of carriers then the status is a 0
                car_status = 0

            # writing the pairs to different columns
            allpair_new_file.write(
                f"{programs}\t{pair1}\t{pair2}\t{car_status}\n")


def is_max_pairs_found(curr_max_pairs: int, new_max_pairs: int) -> int:

    pair_handler_dict = {
        False: 0,
        True: 1
    }

    # This function will return either 1 or zero from the pair handler_dict based on whether or not curr_max_pair from the previous row is less than or greater than the max pairs from the current row
    return pair_handler_dict[curr_max_pairs >= new_max_pairs]


def after_max_pair_found(curr_max_pair: int, new_max_pair: int) -> int:
    '''This function will break out of the loop if a max pair is found and then the size of the pair list starts decreasing'''

    # If the previous rows max number of pairs is higher than the current row then this function will return a one
    if curr_max_pair > new_max_pair:
        return 1
    elif curr_max_pair == new_max_pair:  # If the the above is not true than it returns a 0
        return 0


def combine_output(segment_dir: str, ibd_programs: str, output: str, car_file: str):
    ibd_file_list: list = gather_ibd_files(segment_dir)

    # Setting a max_number of pairs parameter ot use for comparision so that it only keeps one line
    max_pairs: int = 0

    file_dict: dict = build_file_dict(ibd_file_list, ibd_programs)

    for chr_num, variant_id in file_dict.keys():

        file_list: list = file_dict[(chr_num, variant_id)]

        if len(file_list) == 0:  # Checking length of system arguments
            sys.exit(f"no files found for this variant {variant_id}")

        # Making the first argument the output variable
        out = "".join([output, "IBD_", variant_id,
                       "_", chr_num[1:len(chr_num)-1]])
        # next three lines write the files to a dictionary
        files = {}

        for f in file_list:
            print(f)
            files[f.split(':', 1)[0]] = f.split(':', 1)[1]

        print('input {0} files: {1}'.format(
            len(files), ' '.join(files.keys())))

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
        print("printing output path...")
        print(out)
        # sumtab = open(out+'.sum.txt', 'w')
        # uniqtab = open(out+'.uniquq.txt', 'w')

        # allagree = open(out+'.allpair.txt', 'w')
        # # sumtab.write('chr\tpos\tsource\t{}\n'.format(
        # #     '\t'.join(allcomb.keys())))
        # # uniqtab.write('chr\tpos\tsource\t{}\n'.format(
        # #     '\t'.join(allcomb.keys())))
        # allagree.write('chr\tpos\tsegments\tnpair\tpair.list\n')

        oldallpair = set([])

        # this sets a counter to determine how many times the max_pair_int branch is entered
        count: int = 0

        while sum(list(map(lambda f: endtest[f], endtest.keys()))) < len(endtest):
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

            sumrow = list(map(lambda comb: len(
                allinter(comb, curr_pair)), allcomb.values()))
            uniqrow = get_uniqrow(1, allcomb, curr_pair, combtab)
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
                    get_max_pairs(allagree_path, max_pairs_str, variant_id,
                                  car_file, chr_num)

                    break

            if max_pairs_int == 0:
                # This if statement is made to reset the counter if the max_pair integer goes from being 1 back to 0

                # Reseting the counter
                count = 0

            # sumtab.write('{0}\t{1}\t{2}\t{3}\n'.format(str(CHR), str(
            #     pos), ",".join(nowf), '\t'.join(map(str, sumrow))))
            # uniqtab.write('{0}\t{1}\t{2}\t{3}\n'.format(str(CHR), str(
            #     pos), ",".join(nowf), '\t'.join(map(str, uniqrow))))
            max_pairs = len(newallpair)

            print(len(newallpair))

            print(outpair)

            # keeping track of the previous row so that it can be used if necessary
            previous_row_str: str = f"{str(CHR)}\t{str(pos)}\tNA\t{len(newallpair)}\t{' '.join(outpair)}\n"

            # Also keeping track of the base position
            previous_row_bp: str = str(pos)

        # sumtab.close()
        # uniqtab.close()
        # allagree.close()
