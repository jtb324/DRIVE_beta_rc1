# This file contains the functions used to determine how many total individuals contain some variant and how many individuals contain multiple variants

###################################################################################
# importing modules
import pandas as pd
import logging
import glob
import os

###################################################################################
# importing necessary functions from other files

import file_creator_scripts
import population_filter_scripts

###################################################################################
# Function to find the total number of variants

# THIS NEEDS TO BE REVAMPED TO WORK WITH DATAFRAME


def totalVariantIDList(iid_list: set, writeLocation: str, chromo_name: str):
    '''this is a function for the inner loop that will search through each position in the row and when it encouters a one or a two it will add that to the idlist and then return so that the outer loop in the main script moves on to the next row.'''

    logger = logging.getLogger(writeLocation+'/single_variant_analysis.log')

    print("The total number of individual carrier of at least one desired variant is: {}".format(
        len(iid_list)))

    logger.info("The total number of individual carrier of at least one desired variant is: {}".format(
        len(iid_list)))

    if len(chromo_name) == 34:

        file_name_head = chromo_name[21:30]

    elif len(chromo_name) == 35:

        file_name_head = chromo_name[21:31]

    file_name = "".join([file_name_head, ".total_variant_ID_list.txt"])

    writeDirectory = file_creator_scripts.writePath(writeLocation, file_name)

    MyFile = open(
        writeDirectory, 'w')

    for element in iid_list:
        MyFile.write(element)
        MyFile.write('\n')
    MyFile.close()

###########################################################################################


def find_all_files(input_file_path: str):
    '''This function will find a function of all files within a specified directory'''
    cur_dir = os.getcwd()

    os.chdir(input_file_path)
    print(f"current directory: {os.getcwd()}")
    recode_file_list = []

    for file in glob.glob("*.raw"):

        full_file_path = "".join([input_file_path, file])

        recode_file_list.append((full_file_path, file))

    # print(f"{len(recode_file_list)} recoded raw file found with the {os.getcwd()}")

    os.chdir(cur_dir)

    return recode_file_list


###########################################################################################
# This function determines all the individuals who have a specific variant


def singleVariantAnalysis(recodeFile: list, write_path: str, reformat: bool, fileName: str, pop_info: str, pop_code: str):
    '''This function returns a csv containing a list of individuals who carry each variants. It takes a recoded variant file, a path to write the output to, and a file name'''

    recode_file_list = find_all_files(recodeFile[0])

    for file_tuple in recode_file_list:

        if len(file_tuple[1]) == 34:

            file_prefix = file_tuple[1][21:30]

            output_fileName = "".join([file_prefix, ".", fileName])
        elif len(file_tuple[1]) == 35:

            file_prefix = file_tuple[1][21:31]

            output_fileName = "".join([file_prefix, ".", fileName])

        recodeFile = file_tuple[0]

        logger = logging.getLogger(write_path+'/single_variant_analysis.log')

        file_checker = file_creator_scripts.Check_File_Exist(recodeFile, logger)

        raw_file = file_checker.check_file_exist(separator=" ")

        logger.info('Using raw recode file found at {}'.format(recodeFile))

        # subsetting the raw_file for a specific population if the population code, pop_code, is provided
        print(pop_code)
        if pop_code:

            print(f"this is the pop code: {pop_info}")
            print(raw_file)
            dataset_filter = population_filter_scripts.Pop_Filter(pop_info, raw_file)

            pop_info_df, recode_df = dataset_filter.load_files()

            pop_info_subset_df = dataset_filter.get_pop_info_subset(
                pop_info_df, pop_code)

            raw_file = dataset_filter.filter_recode_df(
                pop_info_subset_df, recode_df)

        print(raw_file)
        column_list = list(raw_file.columns[6:].values)

        var_dict = dict()

        var_dict_reformat = dict()

        total_id_set = set()

        for column in column_list:

            iid_list = []

            index_list = raw_file.index[raw_file[column].isin([1.0, 2.0])].tolist()

            for i in index_list:
                iid_list.append(raw_file.loc[i, "IID"])
                total_id_set.add(raw_file.loc[i, "IID"])

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

            reformat_directory = file_creator_scripts.check_dir(write_path, "reformated")

            var_reformat_df.to_csv("".join([
                reformat_directory, "/", file_prefix, ".single_var_list_reformat.csv"]), index=False)

            logger.info('The list of IIDs for each probe id were written to this filepath, {}'.format(
                file_creator_scripts.writePath(reformat_directory, "single_var_list_reformat.csv")
            ))

        elif reformat and not bool(var_dict_reformat):

            print(
                "There were no individuals found so there dictionary was not written to a csv file.")

        totalVariantIDList(total_id_set, write_path,
                           file_tuple[1])

        csv_writer = file_creator_scripts.Csv_Writer_Object(
            var_dict, write_path, output_fileName, logger)

        csv_writer.log_file_path()

        csv_writer.write_to_csv()
