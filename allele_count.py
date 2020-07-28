# Determining the allele count in families
import pandas as pd
import numpy as np
import sys
import logging
import sidetable

################################################
from check_directory import check_dir
from write_path import writePath

################################################


def allele_counts(input_path, fam_file_path, output_path):
    '''This function determines the allele counts for specific variants in each family'''

    logger = logging.getLogger(output_path+"/allele_count_analysis.log")

    #Reading in the files ##########################

    # Reading in the recoded raw file
    try:
        raw_file = pd.read_csv(input_path[0], sep=" ")

    except FileNotFoundError:

        print("The raw recoded file at {} was not found.".format(
            input_path[0]))

        logger.info("The raw recoded file at {} was not found.".format(
            input_path[0]))

        sys.exit(1)

    logger.info("Using the raw recoded file at {}".format(input_path[0]))

    # reading in the file of list of individuals per network
    try:
        networks_df = pd.read_csv(input_path[1], sep=",")

    except FileNotFoundError:

        print("The list of matched network file at {} was not found.".format(
            input_path[1]))

        logger.info("The file containing lists of individuals per network was not found at {}.".format(
            input_path[1]))

        sys.exit(1)

    logger.info("Using the file containing lists of individuals found in each network. This file is found at {}.".format(
        input_path[1]))

    # Reading in the network .fam file
    try:
        pedigree_df = pd.read_csv(fam_file_path, sep="\t")

    except FileNotFoundError:

        print("The full network pedigree file at {} was not found.".format(
            fam_file_path))

        logger.info("The network file at {} was not found.".format(
            fam_file_path))

        sys.exit(1)

    logger.info("Using the network file found at {}.".format(fam_file_path))

    ################################################
    # Getting the names of all the columns to iterate through
    column_list = list(raw_file.columns[6:].values)

    # Creating a dictionary to keep the networks FID, the variant and te number of alleles in each network
    allele_count_dict = dict()

    ###################################################
    # Iterating through all the rows of the pedigree data
    network_list = networks_df["FID"].values.tolist()

    for network in network_list:

        # Next two lines create a subset of the pedigree for a specific network and then it pulls all the
        # iids from that network into a list
        pedigree_subset_df = pedigree_df[pedigree_df["FID"] == network]

        iid_list = pedigree_subset_df["IID"].values.tolist()

        # This subsets the full raw file for only those value also in the row from the pedigree
        raw_file_subset = raw_file[raw_file["IID"].isin(iid_list)]

        # This iterates through the columns of the dataframe row
        for column in column_list:

            ############################################

            variant_allele_frequency = raw_file_subset.stb.freq([column])

            variant_allele_count = 0

            major_allele_count = 0

            try:
                major_allele_count = variant_allele_frequency[variant_allele_frequency[column]
                                                              == 0.0].iloc[0]['count']*2

            except IndexError:
                pass

            if 1.0 in variant_allele_frequency[column].unique():

                subset = variant_allele_frequency[variant_allele_frequency[column] == 1.0]

                subset_count = subset.iloc[0]['count']

                variant_allele_count = subset_count

                major_allele_count += subset_count

            if 2.0 in variant_allele_frequency[column].unique():

                subset = variant_allele_frequency[variant_allele_frequency[column] == 2.0]

                subset_count = subset.iloc[0]['count']

                variant_allele_count = subset_count*2

            ############################################

            # Checks if the index_list is not empty
            if variant_allele_count != 0:

                if "IID" and "Network" and "Variant ID" in allele_count_dict:

                    allele_count_dict["Network"].append(network)

                    allele_count_dict["Variant ID"].append(column)

                    allele_count_dict["Variant Allele Count"].append(
                        int(variant_allele_count))

                    allele_count_dict["Major Allele Count"].append(
                        int(major_allele_count))

                else:

                    allele_count_dict["Network"] = [network]

                    allele_count_dict["Variant ID"] = [column]

                    allele_count_dict["Variant Allele Count"] = [
                        int(variant_allele_count)]

                    allele_count_dict["Major Allele Count"] = [
                        int(major_allele_count)]

    allele_count_df = pd.DataFrame(
        allele_count_dict, columns=["Network", "Variant ID", "Variant Allele Count", "Major Allele Count"])

    # Determining the allele frequency for the variant allele:
    allele_count_df["Variant Allele Frequency"] = allele_count_df['Variant Allele Count'] / \
        (allele_count_df['Variant Allele Count'] +
         allele_count_df['Major Allele Count'])

    reformat_directory = check_dir(output_path, "allele_counts")

    logger.info("Writing two files to the directory {}. The first file contains the variant with the highest number of alleles per network. This file is called 'allele_count.csv'. The second file just drops any variants where there are multiple variants that have the same max allele count. this file is called 'grouped_allele_counts.csv'.")

    allele_count_df.to_csv(
        writePath(reformat_directory, "allele_count.csv"), index=False)

    # Grouping the resulting output by Allele_counts so that there are no duplicate rows for networks that have the same number of alleles for multiple variants

    grouped_df = allele_count_df.groupby(level=0, group_keys=False).apply(
        lambda x: x.loc[x['Network'] == x["Network"].max()])

    grouped_df = grouped_df.drop_duplicates(subset='Network', keep='first')

    grouped_df = grouped_df.drop(["Variant ID", "Major Allele Count",
                                  "Variant Allele Frequency"], axis=1)

    grouped_df.to_csv(
        writePath(reformat_directory, "grouped_allele_counts.csv"), index=False)

    # Grouping the output by Allele count so that we know the number of networks that carry a certain number of variant alleles

    dist_of_alleles_df = grouped_df.stb.freq(['Variant Allele Count'])

    dist_of_alleles_df = dist_of_alleles_df.drop(
        ["cumulative_count", "cumulative_percent"], axis=1)

    dist_of_alleles_df = dist_of_alleles_df.round(2)

    dist_of_alleles_df.to_csv(
        writePath(reformat_directory, "allele_count_distribution.csv"), index=False)
