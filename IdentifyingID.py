# This is the main file that uses the VariantID.py and the MultipleVariantID.py to determine
# the IDs of individuals who contain a variant
import os
from VariantID import individualVariantID
from MultipleVariantID import multiVariantID


def IdentifyID(filepath, writeLocation):
    with open(filepath) as geno_file:

        headerLine = next(geno_file)  # This skips the 1st row

        # This next two lines create lists for the total variants and the multivariants ids
        totalVariantList = []

        multiVariantList = []

        for row in geno_file:  # This iterates through each row in the file

            row = row.split()  # This will split the row by white space
            # Next two lines run the corresponding functions
            individualVariantID(row, totalVariantList)

            multiVariantID(row, multiVariantList)

        # The next two lines print the length of the two list made which corresponds to number of ids found

        print("The number of individuals that have at least one desired variant is:")
        print(len(totalVariantList))
        print("The number of individuals carrying two or multiple desired variants is:")
        print(len(multiVariantList))

        ####################################################
        # This section writes the above two list to txt files
        writePath = writeLocation

        totalVarDirectory = os.path.join(
            writeLocation, "totalVariantIDList.txt")

        multiVarDirectory = os.path.join(
            writeLocation, "multiVariantIDList.txt")

        ##############################################
        # This section writes the total list of variants
        MyFile = open(
            totalVarDirectory, 'w')

        for element in totalVariantList:
            MyFile.write(element)
            MyFile.write('\n')
        MyFile.close()

        ###################################################
        # This creates the list of individuals with mulitple variants
        MyFile = open(
            multiVarDirectory, 'w')

        for element in multiVariantList:
            MyFile.write(element)
            MyFile.write('\n')
        MyFile.close()
