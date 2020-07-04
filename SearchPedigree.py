#################################################
# Import packages
import os
import csv
import pandas as pd
import numpy as np

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


def drop_variant(file, drop_list):
    '''This function will drop specific rows in the dataframe containing certain variants that you specify'''

    file = file.drop(drop_list, axis=0)

    return file
#############################################################################################


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

def searchPedigree(inputPath, outputPath, drop_value, fileName):
    '''This function will search through the provided pedigree file and output two csv files. One file is a list of the variant index positions and a list of individuals with that variant. Then another file is made with the variant index and then the number of individuals that are parts of pedigrees that carry that variant.'''

    # These next lines read in the list of carriers for each variant as a dataframe and can drop any variant from the dataframe
    ind_var_carrier_df = pd.read_csv(
        inputPath[0], sep=",", header=None, names=["MEGA_ID", "IID"])  # Read in the csv will all the carriers of a variant

    if drop_value != None:

        drop_list = ind_var_carrier_df[ind_var_carrier_df["MEGA_ID"].isin(
            drop_value)].index.tolist()

        # This is the function used to drop variants
        ind_var_carrier_df = drop_variant(ind_var_carrier_df, drop_list)

    # This next part reads in the pedigree file
    pedigree_df = pd.read_csv(inputPath[1], sep="\t")

    pedigree_iid_dict = dict()

    for ind in ind_var_carrier_df.index:

        id_list = ind_var_carrier_df.loc[ind, "IID"].strip(
            "[]").replace("'", "").replace(" ", "").split(",")

        index = pedigree_df[pedigree_df["IID"].isin(id_list)].index.tolist()

        iid_in_pedigree_list = pedigree_df.loc[index, "IID"].values.tolist()

        if ind_var_carrier_df.loc[ind, "MEGA_ID"] in pedigree_iid_dict:

            pedigree_iid_dict[ind_var_carrier_df.loc[ind,
                                                     "MEGA_ID"]].append(iid_in_pedigree_list)

        else:

            pedigree_iid_dict[ind_var_carrier_df.loc[ind,
                                                     "MEGA_ID"]] = iid_in_pedigree_list

    csvDictWriter(
        pedigree_iid_dict, outputPath, fileName)

    pedigreeCount(pedigree_iid_dict, outputPath)

####################################################################


# def multiCarriers(inputPath, outputPath, fileName):
#     '''This function will output a list of individuals who carry the same variant and are in the same pedigree.'''

#     multiIndividDict = dict()

#     with open(inputPath[0]) as IndInPedigreeCSV:

#         for row in IndInPedigreeCSV:

#             row = row.split()

#             variant = row[0]

#             IIDList = row[1:]

    #         if len(IIDList) > 1:

    #             for i in range(0, len(IIDList)):

    #                 individual_count = IIDList.count(IIDList[i][0])

    #                 if individual_count > 1:

    #                     if variant in Ped_dict:
    #                         multiIndividDict[variant].update(IIDList[i])

    #                     else:
    #                         multiIndividDict[variant] = {IIDList[i]}
    # csvDictWriter(
    #     multiIndividDict, outputPath, fileName)
