import sys
import gzip
import pandas as pd
import itertools

if len(sys.argv) == 2:  # Checking length of system arguments
    sys.exit('this.py output format1:file1 format2:file2 ... ...')

out = str(sys.argv[1])  # Making the first argument the output variable
# next three lines write the files to a dictionary
files = {}
for f in range(2, len(sys.argv)):
    files[sys.argv[f].split(':')[0]] = sys.argv[f].split(':')[1]

print('input {0} files: {1}'.format(len(files), ' '.join(files.keys())))

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

# define all possible input combination
allcomb = {}
for i in range(len(files.keys()), 0, -1):
    for item in list(itertools.combinations(files.keys(), i)):
        allcomb['+'.join(item)] = item

combtab = pd.DataFrame(0, columns=files.keys(), index=allcomb.keys())
for item in allcomb:
    combtab.loc[item, allcomb[item]] = len(allcomb[item])


def findkey(i, mydict):
    result = []
    offset = -1
    while True:
        try:
            offset = list(mydict.values()).index(i, offset+1)
        except ValueError:
            return list(map(lambda i: list(mydict.keys())[i], result))
        result.append(offset)


def allinter(mylist):  # Finds the pair that intersects for all files
    intu = curr_pair[mylist[0]]
    for f in mylist[1:]:
        intu = intu & curr_pair[f]
    return intu


def get_uniqrow(i):
    uniqdic = {}
    for comb in allcomb.keys():
        raw_n = len(allinter(allcomb[comb]))
        octab = pd.DataFrame(
            map(lambda ff: combtab[ff] > combtab.loc[comb, ff], allcomb[comb])).all()
        overcount = octab.index[octab == True].tolist()
        uniqdic[comb] = raw_n - \
            sum(list(map(lambda oc: uniqdic[oc], overcount)))
    return list(uniqdic.values())


def all_agree_pair(pair_list):  # This is pair that has union for all files
    unionpair = list(pair_list.values())[0]
    for f in pair_list.keys():
        unionpair = unionpair.union(pair_list[f])
    return unionpair


sumtab = open(out+'.sum.txt', 'w')
uniqtab = open(out+'.uniquq.txt', 'w')
allagree = open(out+'.allpair.txt', 'w')
sumtab.write('chr\tpos\tsource\t{}\n'.format('\t'.join(allcomb.keys())))
uniqtab.write('chr\tpos\tsource\t{}\n'.format('\t'.join(allcomb.keys())))
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
            curr_pair[f] = set(map(lambda x: x.split(':')[1], curr_ibd[f]))
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
    sumrow = list(map(lambda comb: len(allinter(comb)), allcomb.values()))
    uniqrow = get_uniqrow(1)
    newallpair = all_agree_pair(curr_pair)
    outpair = []
    for pp in newallpair:
        tool = []
        for f in curr_pair.keys():
            if pp in curr_pair[f]:
                tool.append(f)
        outpair.append('{0}:{1}'.format(','.join(tool), pp))
#    if len(newallpair - oldallpair) > 0:
#        new_add = ['.:'+ s for s in newallpair - oldallpair]
#    else:
#        new_add = ['NA']
#    if len(oldallpair - newallpair) > 0:
#        new_del = ['.:'+ s for s in oldallpair - newallpair]
#    else:
#        new_del = ['NA']
#    if new_add != ['NA'] or new_del != ['NA']:
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
# aa
#{'2:2-5', '1:1-4', '2:1-4', '1:1-3', '2:1-5', '2:3-5', '1:1-2'}
#set(map(lambda x: x.split(':')[1], aa ))
#{'1-3', '3-5', '2-5', '1-5', '1-2', '1-4'}
#
#
