#################################################
# importing packages
import os
import csv
import pandas as pd
import numpy as np
import copy
import logging
import sys

###################################################################################
# importing necessary functions from other files

from NetworkSize import network_sizes
from write_path import writePath
from check_directory import check_dir
from csv_dict_writer import csvDictWriter

###############################################################


def drop_variant(file, drop_list):
    '''This function will drop specific rows in the dataframe containing certain variants that you specify'''

    file = file.drop(drop_list, axis=0)

    return file

#############################################################################################


def pedigreeCount(pedigreeDict, writePath, file_name, logger):
    '''This function will determine the number of individuals found for each combination of variant. It will create a new pedigreeDict where the keys are the index of each variant and the values are the number of individuals containing those variants.'''

    print("generating a file containing the number of individuals found within the pedigree for each variant...")

    # This uses a map function. The .items makes of tuple of key:value pairs and then
    # the lambda function takes the items as a input and updates the original dictionary by
    # by assigning the length of the second element of the tuple to the correct key
    pedigreeDict = dict(map(lambda x: (x[0], len(x[1])), pedigreeDict.items()))

    # This uses the csvDictWriter function to write the individCountDict to a csv file named IndividualCount.csv
    csvDictWriter(pedigreeDict,
                  writePath, file_name, logger)

    individual_sum = sum(pedigreeDict.values())

    logger.info("Using the file {}, the total number of individuals found within the pedigree is {}".format(
        file_name, individual_sum))

    print("The total number of individuals found within the pedigree is {}".format(
        individual_sum))


###############################################################################


def searchPedigree(inputPath, outputPath, drop_value, reformat, pedigree_size, fileName):
    '''This function will search through the provided pedigree file and output two csv files. One file is a list of the variant index positions and a list of individuals with that variant. Then another file is made with the variant index and then the number of individuals that are parts of pedigrees that carry that variant.'''

    logger = logging.getLogger(outputPath+'/search_pedigree_analysis.log')

    # These next lines read in the list of carriers for each variant as a dataframe and can drop any variant from the dataframe
    try:

        ind_var_carrier_df = pd.read_csv(
            inputPath[0], sep=",", header=None, names=["MEGA_ID", "IID"])  # Read in the csv will all the carriers of a variant

    except FileNotFoundError:

        print('There was no file containing a list of carriers found at {}'.format(
            inputPath[1]))

        logger.error('There was no file containing a list of carriers found at {}'.format(
            inputPath[1]))

        sys.exit(1)

    logger.info(
        'Using the list of individual carriers found at {}.'.format(inputPath[0]))

    #This section is purely for dropping any variants that may be overwhelming in the sample ##############

    if drop_value != None:

        drop_list = ind_var_carrier_df[ind_var_carrier_df["MEGA_ID"].isin(
            drop_value)].index.tolist()

        logger.info(
            'The two probes {} were dropped from the results'.format(drop_value))

        # This is the function used to drop variants
        ind_var_carrier_df = drop_variant(ind_var_carrier_df, drop_list)

    ###########################################################################

    # This next part reads in the pedigree file
    try:
        pedigree_df = pd.read_csv(inputPath[1], sep="\t")

    except FileNotFoundError:

        print('The provided network file was not found at {}'.format(
            inputPath[1]))

        logger.error('The provided network file was not found at {}'.format(
            inputPath[1]))

        sys.exit(1)

    logger.info('Using the network file found at {}'.format(inputPath[1]))

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

    full_network_subset = copy.deepcopy(network_subset)

    network_subset = network_subset.query("FID != IID")

    ##################################################
    # This section will

    if reformat == True:

        file_directory = check_dir(outputPath, "reformated")

        single_var_df = pd.read_csv(
            writePath(inputPath[2], "single_var_list_reformat.csv"), sep=",")

        # This will produce a list of every individual in the pedigree if it is wanted
        if pedigree_size == "full":

            full_single_var_df = copy.deepcopy(single_var_df)

            full_single_var_df = full_single_var_df[full_single_var_df["IID"].isin(
                full_network_subset["IID"].values.tolist())]

            full_network_subset_reformat = pd.merge(
                full_network_subset, full_single_var_df, on="IID")

            write_directory = writePath(
                file_directory, "all_ind_in_ped-reformat.csv")

            logger.info("Writing the list of all individuals matched between the provided carrier list and the provided network file including cases where the FID equals the IID to the directory, {}".format(write_directory))

            full_network_subset_reformat.to_csv(
                write_directory, index=False)

            print("producing a file at {} of IIDs of all carriers of interest within\
                  the fam file...".format(write_directory))

        # This section only produces a list of individuals where the FID does not equal the IID

        print("producing a file with only individuals where IID does not equal FID...")

        single_var_df = single_var_df[single_var_df["IID"].isin(
            network_subset["IID"].values.tolist())]

        network_subset_reformat = pd.merge(
            network_subset, single_var_df, on="IID")

        subset_write_directory = writePath(
            file_directory, "ind_in_ped-reformat.csv")

        logger.info("Writing a list of individual carriers found within the provided network file where the FID does not equal the IID to a csv file at {}".format(
            subset_write_directory))

        network_subset_reformat.to_csv(
            subset_write_directory, index=False)

        print("producing a file at {} with only individuals where IID does not \
            equal FID...".format(subset_write_directory))

    #####################################################
    # This section creates the useful python formated documents.

    # This first function writes a dictionary of all individuals found for each variant in the pedigree to a csv file
    logger.info("Writing the dictionary containing all individual carriers found within the provided network file to a csv file titled 'all_ind_in_pedigree.csv'.")

    csvDictWriter(
        pedigree_iid_dict, outputPath, fileName, logger)

    # This function creates a csv file of the number of IIDs found within the .fam file for each variant
    logger.info(
        "Determining the number of all individuals found per variant. This file is found at pedigree_count.csv")

    pedigreeCount(pedigree_iid_dict, outputPath, "pedigree_count.csv", logger)

    # This function outputs a csv file of the number of individuals found per network and it list the IIDs of the individuals found per network

    network_sizes(network_subset, outputPath,
                  "ind_network_counts.csv", "ind_network_list.csv",
                  pedigree_df)

    # This function just passes the pedigree_iid_dict to the next function to determine if there are multiple individuals in a specific network who carry the same variant.
    multi_ind_in_pedigree(
        pedigree_iid_dict, network_subset, outputPath, logger)

####################################################################


def multi_ind_in_pedigree(carrier_dict, pedigree, output_path, logger):
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
    logger.info(
        "Determining which networks contain multiple individuals that carry the same variant...")

    csvDictWriter(carrier_dict, output_path,
                  "multiple_ind_in_pedigree.csv", logger)

    logger.info("Determining a list of multiple individuals that carry the same variants found per network and the number of individuals found per network. These files are found at ind_network_counts.csv and ind_network_list.csv")

    pedigreeCount(carrier_dict, output_path, "multi_ped_count.csv", logger)
