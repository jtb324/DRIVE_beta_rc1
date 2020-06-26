# This file contains the functions used to determine how many toatl individuals contain some variant and how many individuals contain multiple variants

###################################################################################
import os
import csv
###################################################################################
# Function to find the total number of variants
# TODO: figure out how to convert the multiVariantAnalysis to an is in or parallel


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
# This function determines all the individuals who hace a specific variant
# TODO: Want to form a dictionary where the key is a variant and the values are the IIDs of individuals with those variants


# def singleVariant(filepath, writeLocation):
#     with open(filepath) as geno_file:

#         headerLine = next(geno_file)

#         singleVariantDict = dict()

#         singleVariantList = []

#         for row in geno_file:

#             row = row.split()

#             genoRow = row[6:]

#             row = row[0:5]

#             for i in range(1, len(genoRow)):

#                 if genoRow[i] == '1' or genoRow[i] == '2':


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

    # TODO: try to use a map function

    map_iterator = map(len, multiVarDict)
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

        multiVarDict = dict()

        for row in geno_file:  # This iterates through each row in the file

            # This establishes a counter to keep track of how many variants there are in a row
            variantCount = 0

            row = row.split()  # This will split the row by white space

            # These next two lines split the row into the rowID, which contains just the identifying info, and then the genotyping info in the row variable
            rowID = row[0:6]
            row = row[6:]

            # This for loop uses enumerate to get the elements and their index in the row
            # TODO: explore using list comprehension for this instead of loop
            indexList = [i+6 for i,
                         x in enumerate(row) if x == '1' or x == '2']

            # This line checks to see if the row had more than one variant
            if len(indexList) > 1:
                # This converts the indexList to a tuple so that it can be used as a key in the multiVarDict
                indexTuple = tuple(indexList)

                if indexTuple in multiVarDict:  # This checks to see if the indexTuple is already a key in the multiVarDict

                    # If true then it just appends the IID to the value of the multiVarDict
                    multiVarDict[indexTuple].append(rowID[1])

                else:

                    # If false then it creates a new multiVarDict input with that key and value
                    multiVarDict[indexTuple] = [rowID[1]]

                # This then passes the multiVarDict to the csvDictWriter function and ouptuts a csv file
            csvDictWriter(
                multiVarDict, writePath, fileName)

            # This passes the multiVarDict to the individualCount function to determine how many individuals have each combination of variants
    individualCount(multiVarDict, writePath)


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
