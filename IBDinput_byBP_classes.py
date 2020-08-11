#!/usr/bin/python

import sys  # THese are modules used
import gzip
import multiprocessing as mp
import re


min_cM = 3
cpu = 1

####################################################################################################


class newPOS:
    __slots__ = 'add', 'rem'

    def __init__(self, add, rem):
        self.add = add
        self.rem = rem


class Shared_Segment_Convert(newPOS):

    def __init__(self, shared_segment_file, iid_file, output_path, ibd_program_used, min_cM_threshold, thread, base_position):
        self.segment_file = str(shared_segment_file)
        self.iid_file = str(iid_file)
        self.output = str(output_path)
        self.format = str(ibd_program_used)
        self.min_cM = int(min_cM_threshold)
        self.thread = int(thread)
        self.bp = int(base_position)
        # This gets the name of the variant of interest assuming it is input as a text file
        self.variant_name = re.findall(
            "shared_segment_variants/(.*)", self.iid_file)[0][:-4]

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
        id1_indx = int(parameter_dict["id1_indx"])
        id2_indx = int(parameter_dict["id2_indx"])
        chr_indx = int(parameter_dict["chr_indx"])
        str_indx = int(parameter_dict["str_indx"])
        end_indx = int(parameter_dict["end_indx"])
        cM_indx = int(parameter_dict["cM_indx"])

        # This catches the KeyError raised because unit is only found in GERMLINE files
        try:
            unit = parameter_dict["unit"]
        except KeyError:
            pass

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

            elif float(cM) < self.min_cM or ('unit' in vars() and str(line[unit]) != 'cM'):
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

            # start and end not in identified breakpoints
                if int(start) not in IBDindex[CHR]['allpos'] and int(end) not in IBDindex[CHR]['allpos']:
                    IBDdata[CHR][str(start)] = newPOS([pair], [])
                    IBDdata[CHR][str(end)] = newPOS([], [pair])
                    IBDindex[CHR]['allpos'].append(int(start))
                    IBDindex[CHR]['allpos'].append(int(end))
            #
            # start is not in identified breakpoints but end is
                elif int(start) not in IBDindex[CHR]['allpos'] and int(end) in IBDindex[CHR]['allpos']:
                    IBDdata[CHR][str(start)] = newPOS([pair], [])
                    IBDdata[CHR][str(end)].rem.append(str(pair))
                    IBDindex[CHR]['allpos'].append(int(start))
            #
            # start is in identified breakpoints but end not
                elif int(start) in IBDindex[CHR]['allpos'] and int(end) not in IBDindex[CHR]['allpos']:
                    IBDdata[CHR][str(start)].add.append(str(pair))
                    IBDdata[CHR][str(end)] = newPOS([], [pair])
                    IBDindex[CHR]['allpos'].append(int(end))
        #
        # both start and end in identified breakpoints
                elif int(start) in IBDindex[CHR]['allpos'] and int(end) in IBDindex[CHR]['allpos']:
                    IBDdata[CHR][str(start)].add.append(str(pair))
                    IBDdata[CHR][str(end)].rem.append(str(pair))

        data.close()

        # writing output
        print('identified ' +
              str(len(IBDindex[str(i)]['allpos']))+' breakpoints on chr'+str(CHR))

        # Opening the file .small/txt/gz file to write to
        # NEED TO FIX THIS LINE HERE
        write_path = self.output + '_' + self.variant_name + \
            '.chr' + str(CHR) + '.small.txt.gz'

        print(write_path)

        out = gzip.open(write_path, 'wt')

        # Writing the header line to the file
        out.write('chr\tpos\tsegments\tpairs\tadd\tdel\n')

        allibd = set([])

        for pos in sorted(IBDindex[str(i)]['allpos']):

            allibd = allibd | set(IBDdata[str(i)][str(pos)].add)
            allibd = allibd - set(IBDdata[str(i)][str(pos)].rem)

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

        print("entering the IBDsumm function...")

    def run_parallel(self, IBDdata, IBDindex, parameter_dict, uniqID):
        print("initializing the parallel run")

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
