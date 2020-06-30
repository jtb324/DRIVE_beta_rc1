#################################################
# Import packages
import os
import csv
import pandas as pd

#################################################
#################################################
# This section will create a dataframes of all individuals with a specific combination of variants


def writePath(writeLocation, fileName):

    totalVarDirectory = os.path.join(
        writeLocation, fileName)

    return totalVarDirectory

################################################################


def csvDictWriter(multiVarDict, directoryName, fileName):

    # This line opens the csv file with write permissions
    with open(writePath(directoryName, fileName), 'w') as csvfile:

        # this line creates a writer object that the multiVarDict will be written to
        writer = csv.writer(csvfile)

        for key, value in multiVarDict.items():  # This line iterates through the keys and values in the multiVarDict

            # This line writes the keys and values to a row
            writer.writerow([key, value])

    csvfile.close()  # This closes the csv file


###############################################################################
def searchPedigree(inputPath, outputPath, filename):
    IdList = pd.read_csv(
        inputPath[0], skiprows=2, header=None)

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
            Ped_dict, "outPath", "IndInPedigree.csv")
