# This file contains the functions used to determine how many toatl individuals contain some variant and how many individuals contain multiple variants
# BUG: This version using the CLI is 1KB bigger than the version in the scripts folder so I need to figure out why thats happening
###################################################################################
import os
import csv
###################################################################################
# Function to find the total number of variants


def totalVariantID(filepath, writeLocation):
    '''this is a function for the inner loop that will search through each position in the row and when it encouters a one or a two it will add that to the idlist and then return so that the outer loop in the main script moves on to the next row.'''

    with open(filepath) as geno_file:

        headerLine = next(geno_file)  # This skips the 1st row

        # This next two lines create lists for the total variants and the multivariants ids
        totalVariantList = []

        for row in geno_file:  # This iterates through each row in the file

            row = row.split()  # This will split the row by white space

            genoRow = row[6:]

            if '1' in genoRow or '2' in genoRow:

                totalVariantList.append(row[1])

        print("The total number of individual carrier at least one desired variant is: {}".format(
            len(totalVariantList)))

        writeDirectory = writePath(writeLocation, "totalVariantIDList.txt")

        MyFile = open(
            writeDirectory, 'w')

        for element in totalVariantList:
            MyFile.write(element)
            MyFile.write('\n')
        MyFile.close()

############################################################################################
# This function determines the directory to write to


def writePath(writeLocation, fileName):

    totalVarDirectory = os.path.join(
        writeLocation, fileName)

    return totalVarDirectory
############################################################################################
# Function that counts how many individuals carry a set of variants


def individualCount(multiVarDict, writePath):
    '''This function will create a new multiVarDict where the keys are the index of each variant and the values are the number of individuals containing those variants'''

    individCountDict = dict()  # This creates an empty multiVarDict

    for key in multiVarDict:  # This goes through each key of the multiVarDict that was passed into the function
        # This line assigns the key to the key in the new multiVarDict and then finds the value in the old multiVarDict and uses the len() function to determine the # of individuals in the multiVarDict variable which then gets stored as the value in the individCountDict
        individCountDict[key] = len(multiVarDict[key])

    # This uses the csvDictWriter function to write the individCountDict to a csv file named IndividualCount.csv
    csvDictWriter(individCountDict,
                  writePath, "IndividualCount.csv")

################################################################################################
# Function that groups individuals by which variants they carry


def multiVariantAnalysis(filepath, writePath, fileName):
    '''This function preforms the main multiple variant analysis and will make two dictionaries. One multiVarDict contains key that are the index of each variant from the original PLINK recode file (starts at the seventh position because the first 6 values are not important info in this function) and then the values are a list of individuals who contain those variants. The second multiVarDict contains the same keys, but the values are the number of individuals which carry those variants'''

    with open(filepath) as geno_file:

        headerLine = next(geno_file)  # This skips the 1st row

        indexList = []  # This creates an empty indexList that the index of the variants will be added to. It is outside of the while loop so that it gets reset ot empty for each new row
        variantCount = 0

        multiVarDict = dict()

        for row in geno_file:  # This iterates through each row in the file

            row = row.split()  # This will split the row by white space

            i = 6  # This begins the i counter at 6 since the 7th position is the the first index of interest in the row and python is zero based indexing

            # This while loop iterates over the row from the 6th index to the end of the row
            while i <= len(row)-1:

                # Checks to see if the value at row[i] is a string one or two
                if row[i] == '1' or row[i] == '2':

                    variantCount += 1
                    # If the condition is true then it appends i to the index list
                    indexList.append(i)

                    # This then starts another counter that iterates through the rest of the of the row by beginning at index i+1
                    for j in range(i+1, len(row)):

                        # Checks to see if row[j] == string 1 or 2
                        if row[j] == '1' or row[j] == '2':

                            variantCount += 1
                            # If true then it also appends this index to the indexList and then moves to the next j
                            indexList.append(j)

                    if variantCount > 1:
                        # This converts the indexList to a tuple so that it can be used as a key in the multiVarDict
                        indexTuple = tuple(indexList)

                        if indexTuple in multiVarDict:  # This checks to see if the indexTuple is already a key in the multiVarDict

                            # If true then it just appends the IID to the value of the multiVarDict
                            multiVarDict[indexTuple].append(row[1])

                        else:

                            # If false then it creates a new multiVarDict input with that key and value
                            multiVarDict[indexTuple] = [row[1]]

                    # This then sets i to be the length of the row minus 1 so that the while loop will break and the function moves onto the next row.
                    i = len(row)-1

                # adds 1 to index i if the original if statement is false and the row[i] does not equal a string of 1 or 2
                i += 1

                # This then passes the multiVarDict to the csvDictWriter function and ouptuts a csv file
            csvDictWriter(
                multiVarDict, writePath, fileName)

            # This passes the multiVarDict to the individualCount function to determine how many individuals have each combination of variants
    individualCount(multiVarDict)


######################################################################################################
# Function that writes multiVarDict to a csv file

# This function writes the passed multiVarDict to a file whose title is given by the name argument
def csvDictWriter(multiVarDict, directoryName, fileName):

    # This line opens the csv file with write permissions
    with open(writePath(directoryName, fileName), 'w') as csvfile:

        # this line creates a writer object that the multiVarDict will be written to
        writer = csv.writer(csvfile)

        for key, value in multiVarDict.items():  # This line iterates through the keys and values in the multiVarDict

            # This line writes the keys and values to a row
            writer.writerow([key, value])

    csvfile.close()  # This closes the csv file
