# import statements
import sys
import gzip
from numpy.core.numeric import NaN
import pandas as pd
import itertools


class Output_Comparer:

    def __init__(self, output, variant_position, var_of_interest, input_dir):
        self.output = output
        self.variant_pos = variant_position
        self.var_of_interest = var_of_interest
        self.input_dir = input_dir

    def check_arguments(self, args_list):
        if len(args_list) == 0:
            sys.exit('this.py output format1:file1 format2:file2 ... ...')

    def create_file_dict(self, args_list, var_of_interest):

        files = {}
        for f in range(0, len(args_list)):
            print(args_list[f])
            # making the files match the variant
            underscore_pos = args_list[f].split(":")[1].find("_")
            dot_pos = args_list[f].split(":")[1].find(".")
            software_name = args_list[f].split(":")[1][:underscore_pos]
            file_tag = args_list[f].split(":")[1][dot_pos:]

            filename = "".join([software_name, "_", var_of_interest, file_tag])

            full_filename = "".join([self.input_dir, "/", filename])

            # writing the software name and file name to a dictionary
            files[args_list[f].split(':')[0]] = full_filename

        print('input {0} files: {1}'.format(
            len(files), ' '.join(files.keys())))

        return files

    def read_first_line(self, files):
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

        first_line_dict = {
            "curr_pair": curr_pair,
            "curr_pos": curr_pos,
            "curr_ibd": curr_ibd,
            "newpos": newpos,
            "newline": newline,
            "openfile": openfile,
            "endtest": endtest,
            "endline": endline
        }

        return first_line_dict

    def define_input_combinations(self, file_dict):
        allcomb = {}
        for i in range(len(file_dict.keys()), 0, -1):

            for item in list(itertools.combinations(file_dict.keys(), i)):

                allcomb['+'.join(item)] = item

        combtab = pd.DataFrame(
            0, columns=file_dict.keys(), index=allcomb.keys())

        for item in allcomb:
            combtab.loc[item, allcomb[item]] = len(allcomb[item])

        return allcomb, combtab

    def findkey(self, i, new_pos_dict):
        result = []
        offset = -1
        while True:
            try:
                offset = list(new_pos_dict.values()).index(i, offset+1)
            except ValueError:
                return list(map(lambda i: list(new_pos_dict.keys())[i], result))
            result.append(offset)

    def allinter(self, mylist, curr_pair):
        intu = curr_pair[mylist[0]]
        for f in mylist[1:]:
            intu = intu & curr_pair[f]

        return intu

    def get_uniqrow(self, i, allcomb, combtab, curr_pair_dict):
        uniqdic = {}
        for comb in allcomb.keys():
            raw_n = len(self.allinter(allcomb[comb], curr_pair_dict))
            octab = pd.DataFrame(
                map(lambda ff: combtab[ff] > combtab.loc[comb, ff], allcomb[comb])).all()
            overcount = octab.index[octab == True].tolist()
            uniqdic[comb] = raw_n - \
                sum(list(map(lambda oc: uniqdic[oc], overcount)))

        return list(uniqdic.values())

    def all_agree_pair(self, pair_list):
        unionpair = list(pair_list.values())[0]
        for f in pair_list.keys():
            unionpair = unionpair.union(pair_list[f])
        return unionpair

    def write_to_file(self, allcomb, combtab, first_line_dict):

        # Pulling parameters from first_line_dict
        curr_pair = first_line_dict["curr_pair"]
        curr_pos = first_line_dict["curr_pos"]
        curr_ibd = first_line_dict["curr_ibd"]
        newpos = first_line_dict["newpos"]
        newline = first_line_dict["newline"]
        openfile = first_line_dict["openfile"]
        endtest = first_line_dict["endtest"]
        endline = first_line_dict["endline"]

        sumtab = open(self.output+'.sum.txt', 'w')

        uniqtab = open(self.output+'.uniquq.txt', 'w')

        allagree = open(self.output+'.allpair.txt', 'w')

        allagree_path = "".join([self.output, '.allpair.txt'])

        sumtab.write('chr\tpos\tsource\t{}\n'.format(
            '\t'.join(allcomb.keys())))
        uniqtab.write('chr\tpos\tsource\t{}\n'.format(
            '\t'.join(allcomb.keys())))
        allagree.write('chr\tpos\tsegments\tnpair\tpair.list\n')

        oldallpair = set([])
        while sum(list(map(lambda f: endtest[f], endtest.keys()))) < len(endtest):

            pos = min(newpos.values())

            nowf = self.findkey(pos, newpos)
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

            sumrow = list(map(lambda comb: len(
                self.allinter(comb, curr_pair)), allcomb.values()))
            uniqrow = self.get_uniqrow(1, allcomb, combtab, curr_pair)
            newallpair = self.all_agree_pair(curr_pair)
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

        self.identify_most_pairs(allagree_path)

    def identify_most_pairs(self, allpair_file_path):
        print(f"using file at {allpair_file_path}")
        allpairs_df = pd.read_csv(allpair_file_path, sep="\t")

        max_pairs = allpairs_df["npair"].max()

        allpair_df_max_pairs = allpairs_df[allpairs_df["npair"] == max_pairs]

        start_variant_position_df = allpair_df_max_pairs[(
            allpair_df_max_pairs["pos"] < int(self.variant_pos))]["pos"]

        bp_before_variant = start_variant_position_df.iloc[len(
            start_variant_position_df)-1]

        variant_after_position_df = allpairs_df[(
            allpairs_df["pos"] > int(self.variant_pos))]["pos"]

        bp_after_variant = variant_after_position_df.iloc[0]

        dataframe_to_write = allpair_df_max_pairs[allpair_df_max_pairs["pos"]
                                                  == bp_before_variant]

        write_path = "".join([self.output, ".", str(bp_before_variant),
                              "_", str(bp_after_variant), ".allpair.txt"])

        print(f"writing to {write_path}")

        dataframe_to_write.to_csv(
            write_path, header=False, sep=" ", index=None, mode='a', na_rep='NA')

        self.reformat_file(write_path, str(
            bp_before_variant), str(bp_after_variant))

    def reformat_file(self, allpairs_file_path, bp_before_variant: str, bp_after_variant: str):
        '''This function reformats the file so that the first line of the file contains the chromosome id, the bp, and information
        about the number of pairs. Every line after line 1 contains information about the pairs. This file is given the same name 
        as the .allpair.txt file but the ending is changed to the .allpair.new.txt'''
        # Creating a write path for the new file that includes the bp before and after the variant and adds an ending
        # .allpair.new.txt
        write_path = "".join([self.output, ".", bp_before_variant,
                              "_", bp_after_variant, ".allpair.new.txt"])

        # This try and except statement first tries to read teh file in just using a regular delimiter
        # If that fails because of a FileNotFoundError then it says that the file at the provided
        # directory was not found. If it fails for any other reason than it tries a different tab delimiter instead.

        with open(allpairs_file_path, "r") as allpair_file:
            # This section opens the .allpair.new.txt file and writes the reformated text to it.
            for row in allpair_file:

                # print(row)

                row = row.split(" ")

                top_row = row[:4]

                pair_list = row[4:]
                # print(pair_list)

                with open(write_path, "w") as file:

                    for info in top_row:

                        file.write("%s\t" % info)

                    file.write("\n")

                    for pair in pair_list:

                        # Removes extra quotation marks
                        pair = pair.replace('"', '')

                        print(pair)

                        file.write("%s\n" % pair)

                    file.close()
