#!/usr/bin/python

import sys
import getopt
import gzip
import multiprocessing as mp
#import tqdm
#import os

min_cM = 3
cpu = 1
## read argument 
try:
    opts,args = getopt.getopt(sys.argv[1:],"hi:p:o:f:m:t:b:",["help","input=","pheno=","output=","format=","min=","thread=","bp="])
except getopt.GetoptError:
    print('test.py -i <inputfile: prefix($i).(.gz)> -p <ID list> -o <outputfile> -f <format: GERMLINE, iLASH, RaPID, hap-ibd, other> -t <threads> -m <minCM> -b <BP>')
    sys.exit()
for opt,arg in opts:
    if opt in ('-h','--help'):
        print('test.py -i <inputfile> -p <ID list> -o <outputfile> -f <format: GERMLINE, iLASH, RaPID, hap-ibd, other> -t <threads> -m <minCM> -b <BP>')
        sys.exit()
    elif opt in ('-i', '--input'):
#        if 'chr1' in arg:
#            input_pre = arg.split('chr1')[0]
#            input_ext = arg.split('chr1')[1]
#        else:
        input_pre = arg
        input_ext = ''
    elif opt in ('-p', '--pheno'):
        pheno = arg
    elif opt in ('-o', '--output'):
        output = arg
    elif opt in ('-f', '--format'):
        in_format = arg
    elif opt in ('-m', '--min'):
        min_cM = float(arg)
    elif opt in ('-t', '--thread'):
        cpu = min(22, int(arg))
    elif opt in ('-b','--bp'):
        bp= int(arg)

print('Input: {0}'.format(input_pre))
print('Phenotype file: {}'.format(pheno))
print('Output: {}'.format(output))
print('Input file format: {}'.format(in_format))
print('Min output IBD length: {}'.format(min_cM))
## set input format
id1_indx = 0
id2_indx = 2
chr_indx = 4
str_indx = 5
end_indx = 6

if in_format.lower() == 'germline':
    cM_indx = 10
    unit = 11
elif in_format.lower() == 'ilash':
    cM_indx = 9
elif in_format.lower() in ['hap-ibd', 'hapibd']:
    cM_indx = 7
elif in_format.lower() == 'rapid':
    id1_indx = 1
    id2_indx = 2
    chr_indx = 0
    cM_indx = 7
else:
    if in_format.lower().split(':')[0] in ['other','others']:
        indx = in_format.lower().split(':')[1].split(';')
        id1_indx = int(indx[0])-1
        id2_indx = int(indx[1])-1
        chr_indx = int(indx[2])-1
        str_indx = int(indx[3])-1
        end_indx = int(indx[4])-1
        cM_indx = int(indx[5])-1
    else:
        sys.exit('unrecognized or incorrect format: GREMLINE/iLASH/RaPID/hap-ibd/other:id1;id2;chr;str;start bp;end bp;cM')

##read phenotype and build possibele ID pairs
uniqID = {}
dupID = []
pheno_file=open(pheno, 'r')
IDnum=0
for line in pheno_file:
    line = line.strip()
    line = line.split()
    if str(line[0]) in uniqID:
        dupID.append(line[0])
    else:
        uniqID[str(line[0])] = IDnum
        IDnum = IDnum+1

print('identified '+str(len(uniqID))+' unique IDs')
pheno_file.close()

## define function

def find_nearest(pos_int, pos_list):
 nearest_pos = min(pos_list, key=lambda x:abs(x- pos_int))
 if nearest_pos < pos_int:
  return pos_list.index(nearest_pos)+1
 else:
  return pos_list.index(nearest_pos)

#summ = open(output+'.summary.txt','w')
#summ.write('chr\tpos\tn_pairs\n')

class newPOS:
    __slots__ = 'add', 'rem'
    def __init__(self, add, rem):
        self.add = add
        self.rem = rem

##read IBD data
IBDdata = {}
IBDindex = {}
for i in list(range(1, 23)):
#    filepath = input_pre + 'chr' +str(i) + input_ext
#    if input_ext[-3:] == '.gz':
#        data = gzip.open(input_pre + 'chr' +str(i) + input_ext, 'rt')
#        pbar = tqdm(total=os.path.getsize(filepath))
#    else:
#        data = open(input_pre + 'chr'+str(i) +input_ext, 'rt')
#        pbar = tqdm(total=os.path.getsize(filepath))
#    print('reading chr'+str(i)+'...')
    IBDdata[str(i)] = {}
    IBDindex[str(i)] = {'start':999999999, 'end':0, 'allpos':[]}

def IBDsumm(i):
    if input_pre[-3:] == '.gz':
        data = gzip.open(input_pre, 'rt')
    else:
        data = open(input_pre, 'rt')
#    print('working on chr'+str(i)+'...')
    for line in data:
#        pbar.update(len(line.encode('utf-8')))
        line = line.strip()
        line = line.split()
        id1 = str(line[id1_indx])
        id2 = str(line[id2_indx])
        CHR = str(line[chr_indx])
        start = min(int(line[str_indx]), int(line[end_indx]))
        end = max(int(line[str_indx]), int(line[end_indx]))
        cM = str(line[cM_indx])
        if id1 not in uniqID and id2 not in uniqID:
            continue
        elif float(cM)< min_cM or ( 'unit' in vars() and  str(line[unit]) != 'cM'):
            continue
        elif start < bp and end > bp:
            if id1 in uniqID and id2 in uniqID:
                if uniqID[id1] < uniqID[id2]:
                    pair =  '{0}:{1}-{2}'.format(cM, id1, id2)
                else:
                    pair =  '{0}:{1}-{2}'.format(cM, id2, id1)
            
            elif id1 in uniqID:
                pair =  '{0}:{1}-{2}'.format(cM, id1, id2)
            else:
                pair =  '{0}:{1}-{2}'.format(cM, id2, id1)
            print(pair)
## start and end not in identified breakpoints
            if int(start) not in IBDindex[CHR]['allpos'] and int(end) not in IBDindex[CHR]['allpos']:
                IBDdata[CHR][str(start)] = newPOS([pair], [])
                IBDdata[CHR][str(end)] = newPOS( [], [pair])
                IBDindex[CHR]['allpos'].append(int(start))
                IBDindex[CHR]['allpos'].append(int(end))
#                print('new start: {0}; new end: {1}'.format(str(start), str(end)))
## start is not in identified breakpoints but end is
            elif int(start) not in IBDindex[CHR]['allpos'] and int(end) in IBDindex[CHR]['allpos']:
                IBDdata[CHR][str(start)] = newPOS( [pair], [])
                IBDdata[CHR][str(end)].rem.append(str(pair))
                IBDindex[CHR]['allpos'].append(int(start))
#                print('new start: {0}; existed end pairs: {1}'.format(str(start), str(len(IBDdata[CHR][str(end)].rem))))
## start is in identified breakpoints but end not
            elif int(start) in IBDindex[CHR]['allpos'] and int(end) not in IBDindex[CHR]['allpos']:
                IBDdata[CHR][str(start)].add.append(str(pair))
                IBDdata[CHR][str(end)] = newPOS( [], [pair])
                IBDindex[CHR]['allpos'].append(int(end))
#                print('new end: {0}; existed start pairs: {1}'.format(str(end), str(len(IBDdata[CHR][str(start)].add))))
## both start and end in identified breakpoints
            elif int(start) in IBDindex[CHR]['allpos'] and int(end) in IBDindex[CHR]['allpos']:
                IBDdata[CHR][str(start)].add.append(str(pair))
                IBDdata[CHR][str(end)].rem.append(str(pair))
## writing output
    print('identified '+str(len(IBDindex[str(i)]['allpos']))+' breakpoints on chr'+str(CHR))
    out = gzip.open(output +'.chr'+ str(CHR) +'.small.txt.gz','wt')
    out.write('chr\tpos\tsegments\tpairs\tadd\tdel\n')
#    out = open(output+'.chr'+str(i)+'.txt','w')
    allibd = set([])
#    ii = 0
#    tpos = str(len(IBDindex[str(i)]['allpos']))
    for pos in sorted(IBDindex[str(i)]['allpos']):
#        ii = ii + 1
#        print(str(pos)+':'+str(ii)+'/'+tpos)
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
        out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(str(CHR), str(pos), nseg, npair , ' '.join(IBDdata[str(i)][str(pos)].add), ' '.join(IBDdata[str(i)][str(pos)].rem)))
#    summ.write(str(i)+'\t'+str(pos)+'\t'+str(len(IBDdata[str(i)][str(pos)]))+'\n')
    IBDdata[str(i)]=[]
    out.close()

pool = mp.Pool(cpu)
try:
    for i in [7]:
        pool.apply_async(IBDsumm, args=(i,))
except:
    print("Something went wrong")
    pool.close()
    pool.join()
finally:
    pool.close()
    pool.join() 

del(IBDdata)
del(IBDindex)
#del(possible_pairs)
#summ.close()
#out = open(output,'w')
#for CHR in IBDdata.keys():
# for pos in IBDindex[CHR]['allpos']:
#out.write(CHR+'\t'+str(pos)+'\t'+" ".join(IBDdata[CHR][str(pos)])+'\n')
#out.close():
