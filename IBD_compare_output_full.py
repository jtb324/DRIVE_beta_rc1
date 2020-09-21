from pydoc import cli
import sys
import gzip
from numpy.core.numeric import NaN
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

        print(chr_num)
        # Finding the variant id of the file. file names are built so that the variant id sits
        # between the first "_" and the first "."
        for ibd_program in program_list:

            if ibd_program in ibd_file:

                underscore_indx = ibd_file.find(
                    "".join([ibd_program, "_"]))

                ibd_program_indx = underscore_indx+len(ibd_program)+1

        dot_indx = ibd_file.find(".")

        variant_id = ibd_file[ibd_program_indx:dot_indx]

        # There are a few variants that have a dot in the middle and therefore the above code would split it.
        # Therefore if there is no number in the string than it enters this if statement and then finds the next dot.
        if not any(map(str.isdigit, variant_id)):
            # This finds the second dot in the string
            dot_indx = ibd_file.find(".", dot_indx+1)

            # This finds the full variant id
            variant_id = ibd_file[ibd_program_indx:dot_indx]

        # using list comprehension to get all the files that contain that variant and chromosome
        filter_ibd_file_list = [
            file for file in ibd_file_list if variant_id in file and chr_num in file]

        # match up the variants with the IBD program
        for ibd_program in program_list:

            # This goes through the three files in the filter_ibd_file_list
            for file in filter_ibd_file_list:

                # This checks to see if the variant id is in the dictionary and that the ibd program is in the file
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


def get_max_pairs(allpair_file_path: str, output_path: str):
    '''This function will find the highest number of pairs and writes it to a file called .allpair.txt'''
    # Opening the allpair file
    print(allpair_file_path)
    with open(allpair_file_path, "r") as allpair_file:
        max_pairs = 0
        count = 0
        # Going through each row
        for row in allpair_file:

            pair_count = row.split("\t")[3]

            try:

                if int(pair_count) < max_pairs:
                    count += 1
                    if count == 1:
                        end_pos = row.split("\t")[1]

                elif int(pair_count) > max_pairs:
                    max_pairs = int(row.split("\t")[3])
                    previous_row = row
                    start_pos = previous_row.split("\t")[1]

            except ValueError:
                pass

        # This creates the file path up to .allpair.txt. It is upto that point because it will be reused for the allpair.new.txt
        write_path = "".join([output_path, ".", start_pos, "_", end_pos])

        print(write_path)

        with open("".join([write_path, ".allpair.txt"]), "w") as single_allpair_file:
            single_allpair_file.write('chr\tpos\tsegments\tnpair\tpair.list\n')
            single_allpair_file.write(previous_row)

        single_allpair_file.close()
    allpair_file.close()
    reformat("".join([write_path, ".allpair.txt"]), write_path)
    # Need to make the reformat function that takes that and then writes out line by line


def reformat(single_allpair_file: str, write_path: str):
    with open(single_allpair_file, "r") as allpair_file:
        next(allpair_file)
        for row in allpair_file:
            pair_str = row.split("\t")[4]
            pair_list = pair_str.split(" ")
            header = row.split("\t")[0:4]

            with open("".join([write_path, ".allpair.new.txt"]), "w") as allpair_new_file:

                allpair_new_file.write("IBD_programs\tpair_1\tpair_2\n")

                for pair in pair_list:
                    # getting the IBD programs that found the output
                    programs = pair.split(":")[0]

                    # getting pair 1
                    pair1 = pair.split(":")[1].split("-")[0]

                    # getting pair 2
                    pair2 = pair.split(":")[1].split("-")[1]

                    # writing the pairs to different columns
                    allpair_new_file.write(f"{programs}\t{pair1}\t{pair2}\n")


def run(args):
    ibd_file_list = gather_ibd_files(args.segment_dir)

    file_dict = build_file_dict(ibd_file_list, args.ibd_programs)

    for chr_num, variant_id in file_dict.keys():

        file_list = file_dict[(chr_num, variant_id)]

        if len(file_list) == 0:  # Checking length of system arguments
            sys.exit(f"no files found for this variant {variant_id}")

        # Making the first argument the output variable
        out = "".join([args.output, "IBD_", variant_id,
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
        sumtab = open(out+'.sum.txt', 'w')
        uniqtab = open(out+'.uniquq.txt', 'w')
        allagree_path = "".join([out, ".allpair.txt"])
        allagree = open(out+'.allpair.txt', 'w')
        sumtab.write('chr\tpos\tsource\t{}\n'.format(
            '\t'.join(allcomb.keys())))
        uniqtab.write('chr\tpos\tsource\t{}\n'.format(
            '\t'.join(allcomb.keys())))
        allagree.write('chr\tpos\tsegments\tnpair\tpair.list\n')

        oldallpair = set([])

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

            print('chr{0}:{1} from {2}'.format(CHR, str(pos), ' '.join(nowf)))

            sumrow = list(map(lambda comb: len(
                allinter(comb, curr_pair)), allcomb.values()))
            uniqrow = get_uniqrow(1, allcomb, curr_pair, combtab)
            newallpair = all_agree_pair(curr_pair)
            outpair = []
            for pp in newallpair:
                tool = []
                for f in curr_pair.keys():
                    if pp in curr_pair[f]:
                        tool.append(f)
                outpair.append('{0}:{1}'.format(','.join(tool), pp))

            if len(outpair) == 0:
                outpair = ['NA']
            allagree.write('{0}\t{1}\tNA\t{2}\t{3}\n'.format(
                str(CHR), str(pos), len(newallpair), ' '.join(outpair)))

            oldallpair = set(newallpair)

            sumtab.write('{0}\t{1}\t{2}\t{3}\n'.format(str(CHR), str(
                pos), ",".join(nowf), '\t'.join(map(str, sumrow))))
            uniqtab.write('{0}\t{1}\t{2}\t{3}\n'.format(str(CHR), str(
                pos), ",".join(nowf), '\t'.join(map(str, uniqrow))))

        sumtab.close()
        uniqtab.close()
        allagree.close()

        get_max_pairs(allagree_path, out)


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument("-s", help="This argument list the filepath for the ibd files ending in small.txt.gz output",
                        dest="segment_dir", type=str, required=True)

    parser.add_argument("-p", help="This argument list the different ibd programs used",
                        dest="ibd_programs", nargs="+", type=str, required=True)

    parser.add_argument("-o", help="This argument list the output directory",
                        dest="output", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
