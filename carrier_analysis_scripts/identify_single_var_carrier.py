# This file contains the functions used to determine how many total individuals contain some variant and how many individuals contain multiple variants

###################################################################################
# importing modules
import pandas as pd
import logging
import glob
import os
import re

###################################################################################
# importing necessary functions from other files

import file_creator_scripts
import population_filter_scripts

###################################################################################
# Function to find the total number of variants

# THIS NEEDS TO BE REVAMPED TO WORK WITH DATAFRAME


def totalVariantIDList(iid_list: set, writeLocation: str, file_name_head: str):
    """this is a function for the inner loop that will search through each position in the row and when it encouters a one or a two it will add that to the idlist and then return so that the outer loop in the main script moves on to the next row."""

    file_name = "".join([file_name_head, ".total_variant_ID_list.txt"])

    writeDirectory = file_creator_scripts.writePath(writeLocation, file_name)

    MyFile = open(writeDirectory, "w")

    for element in iid_list:
        MyFile.write(element)
        MyFile.write("\n")
    MyFile.close()


###########################################################################################


def find_all_files(input_file_path: str):
    """This function will find a function of all files within a specified directory"""
    print("Finding all the files")
    cur_dir = os.getcwd()

    os.chdir(input_file_path)

    recode_file_list = []

    for file in glob.glob("*.raw"):

        full_file_path = "".join([input_file_path, file])

        recode_file_list.append((full_file_path, file))

    # print(f"{len(recode_file_list)} recoded raw file found with the {os.getcwd()}")

    os.chdir(cur_dir)

    return recode_file_list


###########################################################################################
# This function determines all the individuals who have a specific variant


def singleVariantAnalysis(
    recodeFile: str,
    write_path: str,
    pop_info: str,
    pop_code: str,
):
    """This function returns a csv containing a list of individuals who carry each variants. It takes a recoded variant file, a path to write the output to, and a file name"""
    print("finding all the carriers")

    try:

        # making a directory to put the files from this step in

        os.mkdir("".join([write_path, "carrier_analysis_output/"]))

    except FileExistsError:
        pass

    output_path: str = "".join([write_path, "carrier_analysis_output/"])
    recode_file_list = find_all_files(recodeFile)
    for file_tuple in recode_file_list:
        match = re.search(r".chr\d\d_", file_tuple[1])

        chr_num: str = match.group(0)

        chr_num: str = chr_num.strip(".")

        file_prefix: str = chr_num.strip("_")

        # if len(file_tuple[1]) == 34:

        #     file_prefix = file_tuple[1][21:30]

        #     output_fileName = "".join(
        #         [file_prefix, ".", "single_variant_carrier.csv"])

        # elif len(file_tuple[1]) == 35:

        #     file_prefix = file_tuple[1][21:31]

        output_fileName = "".join(
            [file_prefix, ".", "single_variant_carrier.csv"])

        recodeFile = file_tuple[0]

        file_checker = file_creator_scripts.Check_File_Exist(recodeFile)

        raw_file = file_checker.check_file_exist(separator=" ")

        # subsetting the raw_file for a specific population if the population code, pop_code, is provided

        if pop_code:

            dataset_filter = population_filter_scripts.Pop_Filter(
                pop_info, raw_file)

            pop_info_df, recode_df = dataset_filter.load_files()

            pop_info_subset_df = dataset_filter.get_pop_info_subset(
                pop_info_df, pop_code)

            raw_file = dataset_filter.filter_recode_df(pop_info_subset_df,
                                                       recode_df)

        column_list = list(raw_file.columns[6:].values)

        var_dict = dict()

        var_dict_reformat = dict()

        total_id_set = set()

        for column in column_list:

            iid_list = []

            index_list = raw_file.index[raw_file[column].isin([1.0,
                                                               2.0])].tolist()

            for i in index_list:
                iid_list.append(raw_file.loc[i, "IID"])
                total_id_set.add(raw_file.loc[i, "IID"])

            if (
                    column in var_dict
            ):  # This checks to see if the indexTuple is already a key in the multiVarDict

                # If true then it just appends the IID to the value of the multiVarDict
                var_dict[column].append(iid_list)

            else:

                # If false then it creates a new multiVarDict input with that key and value
                var_dict[column] = iid_list

            for i in index_list:

                if "IID" and "Variant ID" in var_dict_reformat:

                    var_dict_reformat["IID"].append(raw_file.loc[i, "IID"])

                    var_dict_reformat["Variant ID"].append(column)

                else:

                    var_dict_reformat["IID"] = [raw_file.loc[i, "IID"]]
                    var_dict_reformat["Variant ID"] = [column]

        if var_dict_reformat:

            var_reformat_df = pd.DataFrame(var_dict_reformat,
                                           columns=["IID", "Variant ID"])

            reformat_directory = file_creator_scripts.check_dir(
                output_path, "reformated")

            var_reformat_df.to_csv(
                "".join([
                    reformat_directory,
                    "/",
                    file_prefix,
                    ".single_var_list_reformat.csv",
                ]),
                index=False,
            )

        elif not bool(var_dict_reformat):

            print(
                "There were no individuals found so there dictionary was not written to a csv file."
            )

        totalVariantIDList(total_id_set, output_path, file_prefix)

        csv_writer = file_creator_scripts.Csv_Writer_Object(
            var_dict, output_path, output_fileName)

        csv_writer.write_to_csv()
