# This file contains the functions used to determine how many toatl individuals contain some variant and how many individuals contain multiple variants

###################################################################################
import os
import csv
import pandas as pd
import numpy as np
from write_path import writePath
from check_directory import check_dir
###################################################################################
# Function to find the total number of variants


def totalVariantID(recodeFile, writeLocation):
    '''this is a function for the inner loop that will search through each position in the row and when it encouters a one or a two it will add that to the idlist and then return so that the outer loop in the main script moves on to the next row.'''

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

        writeDirectory = writePath(writeLocation, "totalVariantIDList.txt")

        MyFile = open(
            writeDirectory, 'w')

        for element in totalVariantList:
            MyFile.write(element)
            MyFile.write('\n')
        MyFile.close()

############################################################################################
# This function determines all the individuals who hace a specific variant


def singleVariantAnalysis(recodeFile, write_path, reformat, fileName):
    '''This function returns a csv containing a list of individuals who carry each variants. It takes a recoded variant file, a path to write the output to, and a file name'''

    raw_file = pd.read_csv(recodeFile[0], sep=" ")

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

    var_reformat_df = pd.DataFrame(
        var_dict_reformat, columns=["IID", "Variant ID"])

    reformat_directory = check_dir(write_path, "reformated")

    var_reformat_df.to_csv(
        writePath(reformat_directory, "single_var_list_reformat.csv"), index=False)

    csvDictWriter(
        var_dict, write_path, fileName)


############################################################################################
# Function that counts how many individuals carry a set of variants


def individualCount(multiVarDict, writePath):
    '''This function will create a new multiVarDict where the keys are the index of each variant and the values are the number of individuals containing those variants'''

    # This uses a map function. The .items makes of tuple of key:value pairs and then
    # the lambda function takes the items as a input and updates the original dictionary by
    # by assigning the length of the second element of the tuple to the correct key
    multiVarDict = dict(map(lambda x: (x[0], len(x[1])), multiVarDict.items()))

    # This uses the csvDictWriter function to write the individCountDict to a csv file named IndividualCount.csv
    csvDictWriter(multiVarDict,
                  writePath, "IndividualCount.csv")

################################################################################################
# Function that groups individuals by which variants they carry


def multiVariantAnalysis(recodeFile, write_path, reformat, fileName):
    '''This function preforms the main multiple variant analysis and will make two dictionaries. One multiVarDict contains key that are the index of each variant from the original PLINK recode file (starts at the seventh position because the first 6 values are not important info in this function) and then the values are a list of individuals who contain those variants. The second multiVarDict contains the same keys, but the values are the number of individuals which carry those variants'''

    raw_file = pd.read_csv(recodeFile[0], sep=" ")

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

                if "IID" and "Variant List" in multi_var_carriers_reformat:

                    multi_var_carriers_reformat["IID"].append(
                        raw_file.loc[ind, "IID"])
                    multi_var_carriers_reformat["Variant List"].append(
                        index_list)

                else:

                    multi_var_carriers_reformat["IID"] = [
                        raw_file.loc[ind, "IID"]]
                    multi_var_carriers_reformat["Variant List"] = [index_list]

    # This converts the reformated dictionary to a dataframe
    reformat_df = pd.DataFrame(multi_var_carriers_reformat, columns=[
                               "IID", "Variant List"])

    reformat_directory = check_dir(write_path, "reformated")

    # This writes the reformated dataframe as a csv file
    reformat_df.to_csv(
        writePath(reformat_directory, "multi_var_reformated.csv"), index=False)

    # This writes the original formated multivar_carrier dictionary to a csv
    csvDictWriter(
        multi_var_carriers, write_path, fileName)

    # This passes the multiVarDict to the individualCount function to determine how many individuals have each combination of variants
    individualCount(multi_var_carriers, write_path)


######################################################################################################
# Function that writes multiVarDict to a csv file

# This function writes the passed multiVarDict to a file whose title is given by the name argument
def csvDictWriter(multiVarDict, directoryName, fileName):

    # This line opens the csv file with write permissions
    with open(writePath(directoryName, fileName), 'w') as csvfile:

        # this line creates a writer object that the multiVarDict will be written to
        writer = csv.writer(csvfile)

        for key, value in multiVarDict.items():  # This line iterates through the keys and values in the multiVarDict

            # This line writes the keys and values to a row
            writer.writerow([key, value])

    csvfile.close()  # This closes the csv file
