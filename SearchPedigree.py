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


def pedigreeCount(pedigreeDict, writePath):
    '''This function will create a new pedigreeDict where the keys are the index of each variant and the values are the number of individuals containing those variants'''

    print("generating a file containing the number of individuals found within the pedigree for each variant...")
    # This uses a map function. The .items makes of tuple of key:value pairs and then
    # the lambda function takes the items as a input and updates the original dictionary by
    # by assigning the length of the second element of the tuple to the correct key
    pedigreeDict = dict(map(lambda x: (x[0], len(x[1])), pedigreeDict.items()))

    # This uses the csvDictWriter function to write the individCountDict to a csv file named IndividualCount.csv
    csvDictWriter(pedigreeDict,
                  writePath, "pedigreeCount.csv")

    individual_sum = sum(pedigreeDict.values())

    print("The total number of individuals found within the pedigree is {}".format(
        individual_sum))


###############################################################################

def network_sizes(pedigree_file, outputPath):
    '''This function will just output a file that gives you an idea of the size of the network. '''
    print(pedigree_file)
    print("generating a file containing the size of each network...")

    count_directory = writePath(outputPath, "network_counts.csv")

    network_counts = pedigree_file.groupby(
        "FID").count()  # This is a groupby object
    network_counts.reset_index().to_csv(count_directory)

    print("generating a list of all individual carrier in each network...")

    list_directory = writePath(outputPath, "network_list.csv")
    network_list = pedigree_file.groupby(
        "FID")["IID"].apply(list)

    network_list.reset_index().to_csv(list_directory)

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

    index_list = []

    for ind in ind_var_carrier_df.index:

        id_list = ind_var_carrier_df.loc[ind, "IID"].strip(
            "[]").replace("'", "").replace(" ", "").split(",")

        index = pedigree_df[pedigree_df["IID"].isin(id_list)].index.tolist()

        index_list.extend(index)
        iid_in_pedigree_list = pedigree_df.loc[index,
                                               ["FID", "IID"]].values.tolist()

        fid_in_pedigree_list = pedigree_df.loc[index, "FID"].values.tolist()

        if len(iid_in_pedigree_list) > 0:

            if ind_var_carrier_df.loc[ind, "MEGA_ID"] in pedigree_iid_dict:

                pedigree_iid_dict[ind_var_carrier_df.loc[ind,
                                                         "MEGA_ID"]].append(iid_in_pedigree_list)

            else:

                pedigree_iid_dict[ind_var_carrier_df.loc[ind,
                                                         "MEGA_ID"]] = iid_in_pedigree_list

    csvDictWriter(
        pedigree_iid_dict, outputPath, fileName)

    pedigreeCount(pedigree_iid_dict, outputPath)
    print(index_list)
    network_sizes(pedigree_df.loc[np.unique(index_list)], outputPath)

    multi_ind_in_pedigree(
        ind_var_carrier_df, pedigree_df.loc[np.uniqueindex_list], fileName)

####################################################################


def multi_ind_in_pedigree(carrier_df, pedigree_df, fileName):
    '''This function will output a list of individuals who carry the same variant and are in the same pedigree.'''

    print("finished.")


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
