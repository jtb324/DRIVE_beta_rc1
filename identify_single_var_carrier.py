# This file contains the functions used to determine how many total individuals contain some variant and how many individuals contain multiple variants

###################################################################################
# importing modules
from os import path
import pandas as pd
import logging

###################################################################################
# importing necessary functions from other files

from write_path import writePath
from check_directory import check_dir
from csv_writer_class import Csv_Writer_Object
from file_exist_checker import Check_File_Exist
from population_filter import Pop_Filter

###################################################################################
# Function to find the total number of variants

# THIS NEEDS TO BE REVAMPED TO WORK WITH DATAFRAME


def totalVariantID(recodeFile, writeLocation, pop_info, pop_code):
    '''this is a function for the inner loop that will search through each position in the row and when it encouters a one or a two it will add that to the idlist and then return so that the outer loop in the main script moves on to the next row.'''

    logger = logging.getLogger(writeLocation+'/single_variant_analysis.log')

    logger.info('Determining just a list of all individuals who have carriers')

    logger.info('Using raw recode file found at {}'.format(recodeFile[0]))

    recode_file = recodeFile[0]
    # Check to subset the recodeFile
    if pop_code:
        dataset_filter = Pop_Filter(pop_info, recodeFile[0])

        pop_info_df, recode_df = dataset_filter.load_files()

        pop_info_subset_df = dataset_filter.get_pop_info_subset(
            pop_info_df, pop_code)

        recode_file_df = dataset_filter.filter_recode_df(
            pop_info_subset_df, recode_df)

    if path.exists(recodeFile[0]):

        with open(recodeFile[0]) as geno_file:

            headerLine = next(geno_file)  # This skips the 1st row

            # This next two lines create lists for the total variants and the multivariants ids
            totalVariantList = []

            for row in geno_file:  # This iterates through each row in the file

                row = row.split()  # This will split the row by white space

                genoRow = row[6:]

                if '1' in genoRow or '2' in genoRow:

                    totalVariantList.append(row[1])

            print("The total number of individual carrier of at least one desired variant is: {}".format(
                len(totalVariantList)))

            logger.info("The total number of individual carrier of at least one desired variant is: {}".format(
                len(totalVariantList)))

            writeDirectory = writePath(writeLocation, "totalVariantIDList.txt")

            MyFile = open(
                writeDirectory, 'w')

            for element in totalVariantList:
                MyFile.write(element)
                MyFile.write('\n')
            MyFile.close()

    else:

        print("The raw recoded file at {} was not found.".format(
            recodeFile[0]))
        logger.error(
            "The raw recoded file at {} was not found.".format(recodeFile[0]))

############################################################################################
# This function determines all the individuals who have a specific variant


def singleVariantAnalysis(recodeFile, write_path, reformat, fileName, pop_info, pop_code):
    '''This function returns a csv containing a list of individuals who carry each variants. It takes a recoded variant file, a path to write the output to, and a file name'''

    # totalVariantID(recodeFile, write_path)

    logger = logging.getLogger(write_path+'/single_variant_analysis.log')

    file_checker = Check_File_Exist(recodeFile[0], logger)

    raw_file = file_checker.check_file_exist(separator=" ")

    logger.info('Using raw recode file found at {}'.format(recodeFile))

    # subsetting the raw_file for a specific population if the population code, pop_code, is provided

    if pop_code:

        dataset_filter = Pop_Filter(pop_info, raw_file)

        pop_info_df, recode_df = dataset_filter.load_files()

        pop_info_subset_df = dataset_filter.get_pop_info_subset(
            pop_info_df, pop_code)

        raw_file = dataset_filter.filter_recode_df(
            pop_info_subset_df, recode_df)

    column_list = list(raw_file.columns[6:].values)

    var_dict = dict()

    var_dict_reformat = dict()

    iid_list_reformat = []

    variant_list_reformat = []

    for column in column_list:

        iid_list = []

        index_list = raw_file.index[raw_file[column].isin([1, 2])].tolist()

        for i in index_list:
            iid_list.append(raw_file.loc[i, "IID"])

        if column in var_dict:  # This checks to see if the indexTuple is already a key in the multiVarDict

            # If true then it just appends the IID to the value of the multiVarDict
            var_dict[column].append(iid_list)

        else:

            # If false then it creates a new multiVarDict input with that key and value
            var_dict[column] = iid_list

        if reformat == True:

            for i in index_list:

                if "IID" and "Variant ID" in var_dict_reformat:

                    var_dict_reformat["IID"].append(raw_file.loc[i, "IID"])

                    var_dict_reformat["Variant ID"].append(column)

                else:

                    var_dict_reformat["IID"] = [raw_file.loc[i, "IID"]]
                    var_dict_reformat["Variant ID"] = [column]

    if var_dict_reformat:

        var_reformat_df = pd.DataFrame(
            var_dict_reformat, columns=["IID", "Variant ID"])

        reformat_directory = check_dir(write_path, "reformated")

        var_reformat_df.to_csv(
            writePath(reformat_directory, "single_var_list_reformat.csv"), index=False)

        logger.info('The list of IIDs for each probe id were written to this filepath, {}'.format(
            writePath(reformat_directory, "single_var_list_reformat.csv")
        ))

    elif reformat and not bool(var_dict_reformat):

        print("There were no individuals found so there dictionary was not written to a csv file.")

    csv_writer = Csv_Writer_Object(var_dict, write_path, fileName, logger)

    csv_writer.log_file_path()

    csv_writer.write_to_csv()
