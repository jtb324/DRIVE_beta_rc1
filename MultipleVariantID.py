# This file identifies individuals with multiple variants of interest

def multiVariantID(row, variantList):
    '''This function is also for the inner loop and will find individuals that have multiple variants'''

    for i in range(6, (len(row)-1)):

        variantCount = 0

        if row[i] == '1' or row[i] == '2':

            variantCount += 1

    if variantCount > 1:
        variantList.append(row[1])
