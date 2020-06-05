# This file identifies individuals with multiple variants of interest
import os


def multiVariantID(row, variantList):
    '''This function is also for the inner loop and will find individuals that have multiple variants'''

    variantCount = 0

    for i in range(6, (len(row)-1)):

        if row[i] == '1' or row[i] == '2':

            variantCount += 1

    if variantCount > 1:
        variantList.append(row[1])

    return
