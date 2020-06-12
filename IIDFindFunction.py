# This file contains the functions used to determine how many toatl individuals contain some variant and how many individuals contain multiple variants
# BUG: This version using the CLI is 1KB bigger than the version in the scripts folder so I need to figure out why thats happening
###################################################################################
import os

###################################################################################


def totalVariantID(filepath, writeLocation):
    '''this is a function for the inner loop that will search through each position in the row and when it encouters a one or a two it will add that to the idlist and then return so that the outer loop in the main script moves on to the next row.'''

    with open(filepath) as geno_file:

        headerLine = next(geno_file)  # This skips the 1st row

        # This next two lines create lists for the total variants and the multivariants ids
        totalVariantList = []

        for row in geno_file:  # This iterates through each row in the file

            row = row.split()  # This will split the row by white space

            # This for loop loops through the length of the recode row
            i = 6
            while i <= len(row)-1:

                # Checks to see if the element at row[i] is either 1 or 2
                if row[i] == '1' or row[i] == '2':

                    # If it is 1 or 2 then it appends that id to totalVariantList
                    totalVariantList.append(row[1])

                    # This return will then break out of the loop and move on to the next row
                    i = len(row)-1

                i += 1

        writePath = writeLocation

        totalVarDirectory = os.path.join(
            writeLocation, "totalVariantIDList.txt")

        MyFile = open(
            totalVarDirectory, 'w')

        for element in totalVariantList:
            MyFile.write(element)
            MyFile.write('\n')
        MyFile.close()
