#!/usr/bin/python

from Scripts_to_Understand.IBDinput_byBP import IBDdata
import sys  # THese are modules used
import getopt  # This is a commandline parser
import argparse
import gzip
import multiprocessing as mp
#import tqdm
#import os

min_cM = 3
cpu = 1
# read argument
####################################################################################################


class newPOS:
    __slots__ = 'add', 'rem'

    def __init__(self, add, rem):
        self.add = add
        self.rem = rem


class Shared_Segment_Convert(newPOS):

    def __init__(self, shared_segment_file, iid_file, output_path, ibd_program_used, min_cM_threshold, thread, base_position):
        self.segment_file = shared_segment_file
        self.iid_file = iid_file
        self.output = output_path
        self.format = ibd_program_used
        self.min_cM = min_cM_threshold
        self.thread = thread
        self.bp = base_position

        # Printing the initialized
        print('Input: {0}'.format(self.segment_file))
        print('Phenotype file: {}'.format(self.iid_file))
        print('Output: {}'.format(self.output))
        print('Input file format: {}'.format(self.format))
        print('Min output IBD length: {}'.format(self.min_cM))

    def generate_parameters(self):
        '''This will get some of the parameters used later'''

        parameter_dict = {
            "id1_indx": 0,
            "id2_indx": 2,
            "chr_indx": 4,
            "str_indx": 5,
            "end_indx": 6
        }

        # This section checks the in_format argument to determine which software program was used for IBD detection
        if self.format.lower() == 'germline':
            parameter_dict["cM_indx"] = 10
            parameter_dict["unit"] = 11

        elif self.format.lower() == 'ilash':
            parameter_dict["cM_indx"] = 9

        elif self.format.lower() in ['hap-ibd', 'hapibd']:
            parameter_dict["cM_indx"] = 7

        elif self.format.lower() == 'rapid':
            parameter_dict["id1_indx"] = 1
            parameter_dict["id2_indx"] = 2
            parameter_dict["chr_indx"] = 0
            parameter_dict["cM_indx"] = 7

        else:
            if self.format.lower().split(':')[0] in ['other', 'others']:
                indx = self.format.lower().split(':')[1].split(';')
                parameter_dict["id1_indx"] = int(indx[0])-1
                parameter_dict["id2_indx"] = int(indx[1])-1
                parameter_dict["chr_indx"] = int(indx[2])-1
                parameter_dict["str_indx"] = int(indx[3])-1
                parameter_dict["end_indx"] = int(indx[4])-1
                parameter_dict["cM_indx"] = int(indx[5])-1

            else:
                sys.exit(
                    'unrecognized or incorrect format: GREMLINE/iLASH/RaPID/hap-ibd/other:id1;id2;chr;str;start bp;end bp;cM')

        return parameter_dict

    def build_id_pairs(self):
        '''This creates a two list of unique iids and duplicate iids'''
        # read phenotype and build possible ID pairs
        uniqID = {}  # creates empty dictionary
        dupID = []  # creates empty list

        # Opens the file from for the list of IIDs to search through
        pheno_file = open(self.iid_file, 'r')

        IDnum = 0

        for line in pheno_file:  # This goes through each line and will get the id's
            line = line.strip()

            if line in uniqID:
                dupID.append(line[0])
            else:

                uniqID[line] = IDnum
                IDnum = IDnum+1

        print('identified '+str(len(uniqID))+' unique IDs')

        # Closing the file
        pheno_file.close()

        return uniqID, dupID

    def create_ibd_arrays(self):
        '''This creates two IBD arrays that will be used later'''

        # creating a dictionary with 22 key slots and 22 empty dictionaries
        # Also creating a dicitonary IBDindex with 22 dictionaries containing 'start': 999999999, 'end': 0, 'allpos': []
        # Using dictionary comprehension to make the two dictionaries. Just a little more concise than the for loop.
        # The 22 is for the different chromosomes.
        # the "allpos" is the breakpoints
        IBDdata = {str(i): {} for i in range(1, 23)}
        IBDindex = {str(i): {'start': 999999999, 'end': 0, 'allpos': []}
                    for i in range(1, 23)}

        return IBDdata, IBDindex

    def IBDsumm(self, i, IBDdata, IBDindex, parameter_dict, uniqID):
        '''This function will be used in the parallelism function'''

        # undoing the parameter_dict
        id1_indx = parameter_dict["id1_indx"]
        id2_indx = parameter_dict["id2_indx"]
        chr_indx = parameter_dict["chr_indx"]
        str_indx = parameter_dict["str_indx"]
        end_indx = parameter_dict["end_indx"]
        cM_indx = parameter_dict["cM_indx"]
        unit = parameter_dict["unit"]

        # Figuring out if the IBD estimate files are gunzipped
        if self.segment_file[-3:] == '.gz':
            data = gzip.open(self.segment_file, 'rt')
        else:
            # If it is not gunzipped then the file is just normally openned
            data = open(self.segment_file, 'rt')

        for line in data:  # Read through line by line

            line = line.strip()  # This strips all the space out
            line = line.split()  # This then splits based on the whitespace
            id1 = str(line[id1_indx])  # This index finds the first id
            # Then the second id which is at position 3
            id2 = str(line[id2_indx])
            # Then find the chromomsome atq position 5
            CHR = str(line[chr_indx])
            # This just ensures that the smallest start point is giving
            start = min(int(line[str_indx]), int(line[end_indx]))
            # This ensures the maximum endpoint is found
            end = max(int(line[str_indx]), int(line[end_indx]))
            cM = str(line[cM_indx])  # This gives the identified cM length

            if id1 not in uniqID and id2 not in uniqID:
                continue  # ignores ids that are not in the uniqID

            elif float(cM) < min_cM or ('unit' in vars() and str(line[unit]) != 'cM'):
                continue  # also ignores if the cM distance is less then the specified threshold or if the unit is not specified as cM

            elif start < self.bp and end > self.bp:  # Checks to see if the segment starts before the variant position and then if it ends after the variant position

                if id1 in uniqID and id2 in uniqID:  # Checks to see if the ids are in the uniqID list

                    if uniqID[id1] < uniqID[id2]:
                        # If both ids are in the list then it writes the pairs to a variable pair
                        pair = '{0}:{1}-{2}'.format(cM, id1, id2)

                    else:
                        # this just puts the ids in order
                        pair = '{0}:{1}-{2}'.format(cM, id2, id1)

                elif id1 in uniqID:  # If only one id is in the uniqID then it writes it this way with the matched id in
                    pair = '{0}:{1}-{2}'.format(cM, id1, id2)

                else:
                    pair = '{0}:{1}-{2}'.format(cM, id2, id1)

                print(pair)
            # start and end not in identified breakpoints
                if int(start) not in IBDindex[CHR]['allpos'] and int(end) not in IBDindex[CHR]['allpos']:
                    IBDdata[CHR][str(start)] = newPOS([pair], [])
                    IBDdata[CHR][str(end)] = newPOS([], [pair])
                    IBDindex[CHR]['allpos'].append(int(start))
                    IBDindex[CHR]['allpos'].append(int(end))
            #                print('new start: {0}; new end: {1}'.format(str(start), str(end)))
            # start is not in identified breakpoints but end is
                elif int(start) not in IBDindex[CHR]['allpos'] and int(end) in IBDindex[CHR]['allpos']:
                    IBDdata[CHR][str(start)] = newPOS([pair], [])
                    IBDdata[CHR][str(end)].rem.append(str(pair))
                    IBDindex[CHR]['allpos'].append(int(start))
            #                print('new start: {0}; existed end pairs: {1}'.format(str(start), str(len(IBDdata[CHR][str(end)].rem))))
            # start is in identified breakpoints but end not
                elif int(start) in IBDindex[CHR]['allpos'] and int(end) not in IBDindex[CHR]['allpos']:
                    IBDdata[CHR][str(start)].add.append(str(pair))
                    IBDdata[CHR][str(end)] = newPOS([], [pair])
                    IBDindex[CHR]['allpos'].append(int(end))
        #                print('new end: {0}; existed start pairs: {1}'.format(str(end), str(len(IBDdata[CHR][str(start)].add))))
        # both start and end in identified breakpoints
                elif int(start) in IBDindex[CHR]['allpos'] and int(end) in IBDindex[CHR]['allpos']:
                    IBDdata[CHR][str(start)].add.append(str(pair))
                    IBDdata[CHR][str(end)].rem.append(str(pair))

        # writing output
        print('identified ' +
              str(len(IBDindex[str(i)]['allpos']))+' breakpoints on chr'+str(CHR))

        # Opening the file .small/txt/gz file to write to
        out = gzip.open(self.output + '.chr' +
                        str(CHR) + '.small.txt.gz', 'wt')

        # Writing the header line to the file
        out.write('chr\tpos\tsegments\tpairs\tadd\tdel\n')

        allibd = set([])

        for pos in sorted(IBDindex[str(i)]['allpos']):

            allibd = allibd | set(IBDdata[str(i)][str(pos)].add)
            allibd = allibd - set(IBDdata[str(i)][str(pos)].rem)

            print(allibd)

            allibdpair = {}

            if len(IBDdata[str(i)][str(pos)].add) == 0:
                IBDdata[str(i)][str(pos)].add.append('NA')
            if len(IBDdata[str(i)][str(pos)].rem) == 0:
                IBDdata[str(i)][str(pos)].rem.append('NA')

            nseg = str(len(allibd))

            for cM_pair in allibd:
                pair = cM_pair.split(':')[1]
                cM = cM_pair.split(':')[0]

                if pair in allibdpair:

                    allibdpair[pair] = '{0};{1}'.format(allibdpair[pair], cM)

                else:

                    allibdpair[pair] = str(cM)

            npair = str(len(allibdpair))

            out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(str(CHR), str(pos), nseg, npair, ' '.join(
                IBDdata[str(i)][str(pos)].add), ' '.join(IBDdata[str(i)][str(pos)].rem)))

        IBDdata[str(i)] = []

        out.close()

    def run_parallel(self, IBDdata, IBDindex, parameter_dict, uniqID):

        pool = mp.Pool(self.thread)

        try:
            for i in [7]:
                pool.apply_async(self.IBDsumm, args=(
                    i, IBDdata, IBDindex, parameter_dict, uniqID))
        except:
            print("Something went wrong")
            pool.close()
            pool.join()

        finally:
            pool.close()
            pool.join()

        del(IBDdata)
        del(IBDindex)

####################################################################################################
# Initializing the parser


def run(args):
    print("running")

    x = Shared_Segment_Convert(
        args.input, args.pheno, args.output, args.format, args.min, args.thread, args.bp)

    parameter_dict = x.generate_parameters()

    uniqID, dupID = x.build_id_pairs()

    IBDdata, IBDindex = x.create_ibd_arrays()

    x.run_parallel(IBDdata, IBDindex, parameter_dict, uniqID)


def main():
    parser = argparse.ArgumentParser(
        description="This script converts the output of the  IBD detection programs to a human readable form")

    parser.add_argument("-i", '--input', help="This argument just list the path to the IBD detection software output",
                        dest="input", nargs="+", type=str, require=True)

    parser.add_argument("-p", "--pheno", help="This argument just list the path to a txt file of all the IIDs being focused on",
                        dest="pheno", type=str, require=True)

    parser.add_argument("-o", "--output", help="This argument list the output directory",
                        dest="output", type=str, require=True)

    parser.add_argument("-f", "--format", help="This argument specifies the IBD program used",
                        dest="format", nargs="+",  type=str, require=True)

    parser.add_argument("-m", "--min", help="This argument specifies the minimum cM threshold",
                        dest="min", type=str, require=True)

    parser.add_argument("-t", "--thread", help="This argument specifies the threadcount to be used in parallel",
                        dest="thread", type=str, require=True)

    parser.add_argument(
        "-b", "--bp", help="This argument specifies the nucleotide position of the variant of interest", dest="bp", nargs="+", type=str, require=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

#####################################################################################################
# try:
#     # Parses command line options and parameter list. The return value consist of two elements. The first element is a list of (option, value)
#     # pairs and the second is the list of program arguments left after the option list was stripped. The sys.argv[1:] just gets all elements after
#     # the first argument so it drops the leading reference to hte running program.

#     opts, args = getopt.getopt(sys.argv[1:], "hi:p:o:f:m:t:b:", [
#                                "help", "input=", "pheno=", "output=", "format=", "min=", "thread=", "bp="])

# except getopt.GetoptError:
#     # This just recognizes the error thrown if an unrecognized argument is passed.
#     print('test.py -i <inputfile: prefix($i).(.gz)> -p <ID list> -o <outputfile> -f <format: GERMLINE, iLASH, RaPID, hap-ibd, other> -t <threads> -m <minCM> -b <BP>')
#     sys.exit()

# for opt, arg in opts:  # This breaks up the key value pair
#     # This checks to see which inputs are in which argument.
#     if opt in ('-h', '--help'):
#         print('test.py -i <inputfile> -p <ID list> -o <outputfile> -f <format: GERMLINE, iLASH, RaPID, hap-ibd, other> -t <threads> -m <minCM> -b <BP>')
#         sys.exit()
#     elif opt in ('-i', '--input'):
#         #        if 'chr1' in arg:
#         #            input_pre = arg.split('chr1')[0]
#         #            input_ext = arg.split('chr1')[1]
#         #        else:
#         input_pre = arg
#         input_ext = ''
#     elif opt in ('-p', '--pheno'):
#         pheno = arg  # Pheno in this case is the IIDs that are being focused on
#     elif opt in ('-o', '--output'):
#         output = arg
#     elif opt in ('-f', '--format'):
#         # Specifies if the OBD detection program was GERMLINE, hapibd, iLASH, or RaPID
#         in_format = arg
#     elif opt in ('-m', '--min'):
#         min_cM = float(arg)  # shortest IBD segment to be used
#     elif opt in ('-t', '--thread'):
#         cpu = min(22, int(arg))
#     elif opt in ('-b', '--bp'):
#         bp = int(arg)  # This is the pos of the variant of interest

# print('Input: {0}'.format(input_pre))
# print('Phenotype file: {}'.format(pheno))
# print('Output: {}'.format(output))
# print('Input file format: {}'.format(in_format))
# print('Min output IBD length: {}'.format(min_cM))

# set input format
# id1_indx = 0
# id2_indx = 2
# chr_indx = 4
# str_indx = 5
# end_indx = 6

# # This section checks the in_format argument to determine which software program was used for IBD detection
# if in_format.lower() == 'germline':
#     cM_indx = 10
#     unit = 11
# elif in_format.lower() == 'ilash':
#     cM_indx = 9
# elif in_format.lower() in ['hap-ibd', 'hapibd']:
#     cM_indx = 7
# elif in_format.lower() == 'rapid':
#     id1_indx = 1
#     id2_indx = 2
#     chr_indx = 0
#     cM_indx = 7
# else:
#     if in_format.lower().split(':')[0] in ['other', 'others']:
#         indx = in_format.lower().split(':')[1].split(';')
#         id1_indx = int(indx[0])-1
#         id2_indx = int(indx[1])-1
#         chr_indx = int(indx[2])-1
#         str_indx = int(indx[3])-1
#         end_indx = int(indx[4])-1
#         cM_indx = int(indx[5])-1
#     else:
#         sys.exit(
#             'unrecognized or incorrect format: GREMLINE/iLASH/RaPID/hap-ibd/other:id1;id2;chr;str;start bp;end bp;cM')

# read phenotype and build possible ID pairs
# uniqID = {}  # creates empty dictionary
# dupID = []  # creates empty list
# # Opens the file from for the list of IIDs to search through
# pheno_file = open(pheno, 'r')
# IDnum = 0

# for line in pheno_file:  # This goes through each line and will get the id's
#     line = line.strip()
#     # line = line.split()
#     # if str(line[0]) in uniqID:
#     if line in uniqID:
#         dupID.append(line[0])
#     else:
#         # uniqID[str(line[0])] = IDnum
#         uniqID[line] = IDnum
#         IDnum = IDnum+1

# print('identified '+str(len(uniqID))+' unique IDs')
# pheno_file.close()

# define function

###############################################################
# This function seems to not be used


# def find_nearest(pos_int, pos_list):
#     nearest_pos = min(pos_list, key=lambda x: abs(x - pos_int))
#     if nearest_pos < pos_int:
#         return pos_list.index(nearest_pos)+1
#     else:
#         return pos_list.index(nearest_pos)

################################################################

#summ = open(output+'.summary.txt','w')
# summ.write('chr\tpos\tn_pairs\n')


# read IBD data
# IBDdata = {}
# IBDindex = {}

# ###################################################################
# # WHAT IS THIS FOR?
# # creating a dictionary with 22 key slots and 22 empty dictionaries
# # Also creating a dicitonary IBDindex with 22 dictionaries containing 'start': 999999999, 'end': 0, 'allpos': []
# # Using dictionary comprehension to make the two dictionaries. Just a little more concise than the for loop.
# # The 22 is for the different chromosomes.
# # the "allpos" is the breakpoints
# IBDdata = {str(i): {} for i in range(1, 23)}
# IBDindex = {str(i): {'start': 999999999, 'end': 0, 'allpos': []}
#             for i in range(1, 23)}

# for i in list(range(1, 23)):
#     #    filepath = input_pre + 'chr' +str(i) + input_ext
#     #    if input_ext[-3:] == '.gz':
#     #        data = gzip.open(input_pre + 'chr' +str(i) + input_ext, 'rt')
#     #        pbar = tqdm(total=os.path.getsize(filepath))
#     #    else:
#     #        data = open(input_pre + 'chr'+str(i) +input_ext, 'rt')
#     #        pbar = tqdm(total=os.path.getsize(filepath))
#     #    print('reading chr'+str(i)+'...')
#     IBDdata[str(i)] = {}
#     IBDindex[str(i)] = {'start': 999999999, 'end': 0, 'allpos': []}
####################################################################


# def IBDsumm(i):
#     if input_pre[-3:] == '.gz':  # Figuring out if the IBD estimate files are gunzipped
#         data = gzip.open(input_pre, 'rt')
#     else:
#         # If it is not gunzipped then the file is just normally openned
#         data = open(input_pre, 'rt')
#         #print('working on chr'+str(i)+'...')
#     for line in data:  # Read through line by line
#         #        pbar.update(len(line.encode('utf-8')))
#         line = line.strip()  # This strips all the space out
#         line = line.split()  # This then splits based on the whitespace
#         id1 = str(line[id1_indx])  # This index finds the first id
#         id2 = str(line[id2_indx])  # Then the second id which is at position 3
#         CHR = str(line[chr_indx])  # Then find the chromomsome atq position 5
#         # This just ensures that the smallest start point is giving
#         start = min(int(line[str_indx]), int(line[end_indx]))
#         # This ensures the maximum endpoint is found
#         end = max(int(line[str_indx]), int(line[end_indx]))
#         cM = str(line[cM_indx])  # This gives the identified cM length
#         if id1 not in uniqID and id2 not in uniqID:
#             continue  # ignores ids that are not in the uniqID
#         elif float(cM) < min_cM or ('unit' in vars() and str(line[unit]) != 'cM'):
#             continue  # also ignores if the cM distance is less then the specified threshold or if the unit is not specified as cM
#         elif start < bp and end > bp:  # Checks to see if the segment starts before the variant position and then if it ends after the variant position
#             if id1 in uniqID and id2 in uniqID:  # Checks to see if the ids are in the uniqID list
#                 if uniqID[id1] < uniqID[id2]:
#                     # If both ids are in the list then it writes the pairs to a variable pair
#                     pair = '{0}:{1}-{2}'.format(cM, id1, id2)
#                 else:
#                     # this just puts the ids in order
#                     pair = '{0}:{1}-{2}'.format(cM, id2, id1)

#             elif id1 in uniqID:  # If only one id is in the uniqID then it writes it this way with the matched id in
#                 pair = '{0}:{1}-{2}'.format(cM, id1, id2)
#             else:
#                 pair = '{0}:{1}-{2}'.format(cM, id2, id1)
#             print(pair)
#         # start and end not in identified breakpoints
#             if int(start) not in IBDindex[CHR]['allpos'] and int(end) not in IBDindex[CHR]['allpos']:
#                 IBDdata[CHR][str(start)] = newPOS([pair], [])
#                 IBDdata[CHR][str(end)] = newPOS([], [pair])
#                 IBDindex[CHR]['allpos'].append(int(start))
#                 IBDindex[CHR]['allpos'].append(int(end))
#         #                print('new start: {0}; new end: {1}'.format(str(start), str(end)))
#         # start is not in identified breakpoints but end is
#             elif int(start) not in IBDindex[CHR]['allpos'] and int(end) in IBDindex[CHR]['allpos']:
#                 IBDdata[CHR][str(start)] = newPOS([pair], [])
#                 IBDdata[CHR][str(end)].rem.append(str(pair))
#                 IBDindex[CHR]['allpos'].append(int(start))
#         #                print('new start: {0}; existed end pairs: {1}'.format(str(start), str(len(IBDdata[CHR][str(end)].rem))))
#         # start is in identified breakpoints but end not
#             elif int(start) in IBDindex[CHR]['allpos'] and int(end) not in IBDindex[CHR]['allpos']:
#                 IBDdata[CHR][str(start)].add.append(str(pair))
#                 IBDdata[CHR][str(end)] = newPOS([], [pair])
#                 IBDindex[CHR]['allpos'].append(int(end))
#     #                print('new end: {0}; existed start pairs: {1}'.format(str(end), str(len(IBDdata[CHR][str(start)].add))))
#     # both start and end in identified breakpoints
#             elif int(start) in IBDindex[CHR]['allpos'] and int(end) in IBDindex[CHR]['allpos']:
#                 IBDdata[CHR][str(start)].add.append(str(pair))
#                 IBDdata[CHR][str(end)].rem.append(str(pair))
#     # writing output
#     print('identified ' +
#           str(len(IBDindex[str(i)]['allpos']))+' breakpoints on chr'+str(CHR))
#     # Opening the file .small/txt/gz file to write to
#     out = gzip.open(output + '.chr' + str(CHR) + '.small.txt.gz', 'wt')
#     # Writing the header line to the file
#     out.write('chr\tpos\tsegments\tpairs\tadd\tdel\n')
#     #    out = open(output+'.chr'+str(i)+'.txt','w')s
#     allibd = set([])
#     #    ii = 0
#     #    tpos = str(len(IBDindex[str(i)]['allpos']))
#     for pos in sorted(IBDindex[str(i)]['allpos']):
#         #        ii = ii + 1
#         #        print(str(pos)+':'+str(ii)+'/'+tpos)
#         allibd = allibd | set(IBDdata[str(i)][str(pos)].add)
#         allibd = allibd - set(IBDdata[str(i)][str(pos)].rem)
#         print(allibd)
#         allibdpair = {}
#         if len(IBDdata[str(i)][str(pos)].add) == 0:
#             IBDdata[str(i)][str(pos)].add.append('NA')
#         if len(IBDdata[str(i)][str(pos)].rem) == 0:
#             IBDdata[str(i)][str(pos)].rem.append('NA')
#         nseg = str(len(allibd))
#         for cM_pair in allibd:
#             pair = cM_pair.split(':')[1]
#             cM = cM_pair.split(':')[0]
#             if pair in allibdpair:
#                 allibdpair[pair] = '{0};{1}'.format(allibdpair[pair], cM)
#             else:
#                 allibdpair[pair] = str(cM)
#         npair = str(len(allibdpair))
#         out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(str(CHR), str(pos), nseg, npair, ' '.join(
#             IBDdata[str(i)][str(pos)].add), ' '.join(IBDdata[str(i)][str(pos)].rem)))
#     #     summ.write(str(i)+'\t'+str(pos)+'\t'+str(len(IBDdata[str(i)][str(pos)]))+'\n')
#     IBDdata[str(i)] = []
#     out.close()


# pool = mp.Pool(cpu)
# try:
#     for i in [7]:
#         pool.apply_async(IBDsumm, args=(i,))
# except:
#     print("Something went wrong")
#     pool.close()
#     pool.join()
# finally:
#     pool.close()
#     pool.join()

# del(IBDdata)
# del(IBDindex)
# del(possible_pairs)
# summ.close()
#out = open(output,'w')
# for CHR in IBDdata.keys():
# for pos in IBDindex[CHR]['allpos']:
#out.write(CHR+'\t'+str(pos)+'\t'+" ".join(IBDdata[CHR][str(pos)])+'\n')
# out.close():
