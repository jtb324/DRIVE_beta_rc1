#################################################
# Import packages
import os
import csv
import pandas as pd

###############################################################


def writePath(writeLocation, fileName):
    '''This section creates a function that creates a write path that can be used. This helps keep the code DRY'''

    totalVarDirectory = os.path.join(
        writeLocation, fileName)

    return totalVarDirectory

################################################################


def csvDictWriter(variantDict, directoryName, fileName):
    '''This function writes a dictionary to a csv file. The function takes a dictionary, a directory path, and a fileName as inputs and outputs a csv file at that directory'''

    # This line opens the csv file with write permissions
    with open(writePath(directoryName, fileName), 'w') as csvfile:

        # this line creates a writer object that the multiVarDict will be written to
        writer = csv.writer(csvfile)

        for key, value in variantDict.items():  # This line iterates through the keys and values in the multiVarDict

            # This line writes the keys and values to a row
            writer.writerow([key, value])

    csvfile.close()  # This closes the csv file

###############################################################################


def pedigreeCount(multiVarDict, writePath):
    '''This function will create a new multiVarDict where the keys are the index of each variant and the values are the number of individuals containing those variants'''

    # This uses a map function. The .items makes of tuple of key:value pairs and then
    # the lambda function takes the items as a input and updates the original dictionary by
    # by assigning the length of the second element of the tuple to the correct key
    multiVarDict = dict(map(lambda x: (x[0], len(x[1])), multiVarDict.items()))

    # This uses the csvDictWriter function to write the individCountDict to a csv file named IndividualCount.csv
    csvDictWriter(multiVarDict,
                  writePath, "pedigreeCount.csv")


###############################################################################
def searchPedigree(inputPath, outputPath, fileName):
    '''This function will search through the provided pedigree file and output two csv files. One file is a list of the variant index positions and a list of individuals with that variant. Then another file is made with the variant index and then the number of individuals that are parts of pedigrees that carry that variant.'''

    IdList = pd.read_csv(
        inputPath[0], skiprows=2, header=None)  # TODO: need to adjust the skiprows so that it is customizable depending on inputs.

    Ped_dict = dict()

    with open(inputPath[1]) as pedigree_file:
        for row in pedigree_file:

            row = row.split()

            for i, j in IdList.iterrows():

                if row[1] in j[1] and row[1] != row[0]:

                    if j[0] in Ped_dict:
                        Ped_dict[j[0]].update({row[0]: row[1]})
                    else:
                        Ped_dict[j[0]] = {row[0]: row[1]}

        csvDictWriter(
            Ped_dict, outputPath, fileName)

    pedigreeCount(Ped_dict, outputPath)
