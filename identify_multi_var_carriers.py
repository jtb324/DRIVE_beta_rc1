##########################################################################
# importing modules
import sys
import pandas as pd
import logging

##########################################################################
# importing necessary functions from other files

from write_path import writePath
from check_directory import check_dir
from csv_writer_class import Csv_Writer_Object
from file_exist_checker import Check_File_Exist

##########################################################################

# Function that groups individuals by which variants they carry


def multiVariantAnalysis(recodeFile, write_path, reformat, fileName):
    '''This function preforms the main multiple variant analysis and will make two dictionaries. One multiVarDict contains key that are the index of each variant from the original PLINK recode file (starts at the seventh position because the first 6 values are not important info in this function) and then the values are a list of individuals who contain those variants. The second multiVarDict contains the same keys, but the values are the number of individuals which carry those variants'''

    logger = logging.getLogger(write_path+'/multi_variant_analysis.log')

    #Reading in the file ########################################

    file_checker = Check_File_Exist(recodeFile[0], logger)

    raw_file = file_checker.check_file_exist()

    logger.info('Using the raw recode file at {}'.format(recodeFile))

    #############################################################

    column_list = list(raw_file.columns[6:].values)

    multi_var_carriers = dict()

    multi_var_carriers_reformat = dict()

    for ind in raw_file.index:

        index_1 = raw_file.loc[ind, column_list][raw_file.loc[ind,
                                                              column_list] == 1].index.tolist()

        index_2 = raw_file.loc[ind, column_list][raw_file.loc[ind,
                                                              column_list] == 2]. index.tolist()

        index_list = index_1 + index_2

        index_tuple = tuple(index_list)

        if len(index_tuple) > 1:

            if index_tuple in multi_var_carriers:

                multi_var_carriers[index_tuple].append(
                    raw_file.loc[ind, "IID"])

            else:

                multi_var_carriers[index_tuple] = [raw_file.loc[ind, "IID"]]

            if reformat == True:

                for index in index_tuple:

                    if "IID" and "Variant List" in multi_var_carriers_reformat:

                        multi_var_carriers_reformat["IID"].append(
                            raw_file.loc[ind, "IID"])
                        multi_var_carriers_reformat["Variant List"].append(
                            index)

                    else:

                        multi_var_carriers_reformat["IID"] = [
                            raw_file.loc[ind, "IID"]]
                        multi_var_carriers_reformat["Variant List"] = [
                            index]

    # This converts the reformated dictionary to a dataframe
    if multi_var_carriers_reformat:

        reformat_df = pd.DataFrame(multi_var_carriers_reformat, columns=[
            "IID", "Variant List"])

        reformat_directory = check_dir(write_path, "reformated")

        # This writes the reformated dataframe as a csv file
        reformat_df.to_csv(
            writePath(reformat_directory, "multi_var_reformated.csv"), index=False)

        logging.info('reformated file written to {}'.format(
            writePath(reformat_directory, "multi_var_reformated.csv")))

    elif reformat and not bool(multi_var_carriers_reformat):

        print("There were no individuals found within the reformated csv so the dictionary was not written to a csv file.")

    # This section writes the original formated multivar_carrier dictionary to a csv
    logger.info(
        'Writing the number of individuals who carry multiple variants to a csv file...')

    csv_writer = Csv_Writer_Object(
        multi_var_carriers, write_path, fileName, logger)

    csv_writer.log_file_path()

    csv_writer.write_to_csv()

    # This passes the multiVarDict to the individualCount function to determine how many individuals have each combination of variants
    individualCount(multi_var_carriers, write_path, logger)

############################################################################################
# Function that counts how many individuals carry a set of variants


def individualCount(multiVarDict, write_path, logger):
    '''This function will create a new multiVarDict where the keys are the index of each variant and the values are the number of individuals containing those variants'''

    # This uses a map function. The .items makes of tuple of key:value pairs and then
    # the lambda function takes the items as a input and updates the original dictionary by
    # by assigning the length of the second element of the tuple to the correct key
    multiVarDict = dict(map(lambda x: (x[0], len(x[1])), multiVarDict.items()))

    # This uses the csvDictWriter function to write the individCountDict to a csv file named IndividualCount.csv

    csv_writer = Csv_Writer_Object(multiVarDict,
                                   write_path, "IndividualCount.csv", logger)

    csv_writer.log_file_path()

    csv_writer.write_to_csv()
