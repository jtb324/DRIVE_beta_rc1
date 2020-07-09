#################################################
# Import packages
import os
import csv
import pandas as pd
import numpy as np
from NetworkSize import determine_network_sizes
from write_path import writePath
from check_directory import check_dir
###############################################################


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


def pedigreeCount(pedigreeDict, writePath, file_name):
    '''This function will create a new pedigreeDict where the keys are the index of each variant and the values are the number of individuals containing those variants'''

    print("generating a file containing the number of individuals found within the pedigree for each variant...")
    # This uses a map function. The .items makes of tuple of key:value pairs and then
    # the lambda function takes the items as a input and updates the original dictionary by
    # by assigning the length of the second element of the tuple to the correct key
    pedigreeDict = dict(map(lambda x: (x[0], len(x[1])), pedigreeDict.items()))

    # This uses the csvDictWriter function to write the individCountDict to a csv file named IndividualCount.csv
    csvDictWriter(pedigreeDict,
                  writePath, file_name)

    individual_sum = sum(pedigreeDict.values())

    print("The total number of individuals found within the pedigree is {}".format(
        individual_sum))


###############################################################################

def network_sizes(pedigree_subset, output_path, counts_file_name, list_file_name, full_pedigree_file):
    '''This function will just output a file that gives you an idea of the size of the network. '''

    print("generating a file containing the size of each network...")

    count_directory = writePath(output_path, counts_file_name)

    network_counts = pedigree_subset.groupby(
        "FID").count()  # This is a groupby object

    network_counts.reset_index().to_csv(count_directory, index=False)

    print("generating a list of all individual carrier in each network...")

    list_directory = writePath(output_path, list_file_name)

    network_list = pedigree_subset.groupby(
        "FID")["IID"].apply(list)

    network_list.reset_index().to_csv(list_directory, index=False)

    determine_network_sizes(count_directory, full_pedigree_file, output_path)

###############################################################################


def searchPedigree(inputPath, outputPath, drop_value, reformat, fileName):
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

    ############################
    # Creating necessary global variables
    pedigree_iid_dict = dict()

    index_list = []
    ###########################

    for ind in ind_var_carrier_df.index:

        id_list = ind_var_carrier_df.loc[ind, "IID"].strip(
            "[]").replace("'", "").replace(" ", "").split(",")

        index = pedigree_df[pedigree_df["IID"].isin(id_list)].index

        index_list.extend(index)

        iid_in_pedigree_list = pedigree_df.loc[index,
                                               ["IID"]].values
        if len(iid_in_pedigree_list) > 0:

            iid_in_pedigree_list = ' '.join(map(str, iid_in_pedigree_list)).strip(
                "[]").replace("] [", " ").replace("'", "").split(" ")

            pedigree_iid_dict[ind_var_carrier_df.loc[ind,
                                                     "MEGA_ID"]] = iid_in_pedigree_list

    ##################################################
    # This creates a subset of the network file where the FID does not equal the IID. This is because the individuals are only part of a network if the FID is different than the IID

    network_subset = pedigree_df.loc[np.unique(
        index_list)]

    network_subset = network_subset.query("FID != IID")

    ##################################################
    # This section will

    if reformat == True:

        file_directory = check_dir(outputPath, "reformated")

        print("producing a file more conducive output file format for ")

        single_var_df = pd.read_csv(
            writePath(inputPath[2], "single_var_list_reformat.csv"), sep=",")

        single_var_df = single_var_df[single_var_df["IID"].isin(
            network_subset["IID"].values.tolist())]

        network_subset_reformat = pd.merge(
            network_subset, single_var_df, on="IID")

        network_subset_reformat.to_csv(
            writePath(file_directory, "ind_in_ped-reformat.csv"), index=False)

    #####################################################
    csvDictWriter(
        pedigree_iid_dict, outputPath, fileName)

    pedigreeCount(pedigree_iid_dict, outputPath, "pedigree_count.csv")

    network_sizes(network_subset, outputPath,
                  "ind_network_counts.csv", "ind_network_list.csv",
                  pedigree_df)

    multi_ind_in_pedigree(
        pedigree_iid_dict, network_subset, outputPath)

####################################################################


def multi_ind_in_pedigree(carrier_dict, pedigree, output_path):
    '''This function will output a list of individuals who carry the same variant and are in the same pedigree.'''

    print("generating a csv file of networks containing multiple individuals who carry the same variant...")

    carrier_df = pd.DataFrame.from_dict(
        carrier_dict, orient="index").reset_index()

    carrier_dict = dict()

    for ind in carrier_df.index:

        id_list = carrier_df.loc[ind].dropna().values.tolist()

        if len(id_list) > 1:

            pedigree_df = pedigree[pedigree["IID"].isin(id_list)]

            ped_count_df = pedigree_df.groupby(
                "FID").count().reset_index()

            fid_df = ped_count_df[ped_count_df["IID"] > 1]

            if not fid_df.empty:  # This line skips fid_df where there are no FIDs that meet the above condition

                for fid_ind in fid_df.index:

                    fid = fid_df.loc[fid_ind, "FID"]

                    indexTuple = tuple(
                        [carrier_df.loc[ind, "index"],
                         fid]
                    )

                    carrier_dict[indexTuple] = pedigree_df["IID"][pedigree_df["FID"]
                                                                  == fid].values.tolist()

    csvDictWriter(carrier_dict, output_path, "multiple_ind_in_pedigree.csv")

    pedigreeCount(carrier_dict, output_path, "multi_ped_count.csv")
