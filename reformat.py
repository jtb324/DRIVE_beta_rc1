# This script is designed to take the output and provide something more like davids ideal format
import pandas as pd
import glob
import os
import re
import argparse
import numpy as np

# need to get list of carriers
# First function will be designed to get a glob of files


def get_files(directory: str, file_id: str) -> list:

    cur_dir: str = os.getcwd()
    os.chdir(directory)

    file_list: list = []

    for file in glob.glob(file_id):

        full_file_path: str = "".join([directory, file])

        file_list.append(full_file_path)

    os.chdir(cur_dir)

    return file_list


def form_genotype_df(map_file_path: str, ped_file_path: str) -> pd.DataFrame:
    '''This function will form a dictionary of dictionaries where the outer key is the iid, 
    the inner key is the variant and the inner value is the genotype string'''

    # form three list for the iid, the variant, and the genotype string
    iid_list: list = []
    variant_list: list = []
    geno_list: list = []

    # opening the ped file
    with open(ped_file_path, "r") as ped_file:

        # iterating through row
        for row in ped_file:

            # split row
            row: list = row.split(" ")

            # getting only the genotype portion of the row

            geno_row: list = row[6:]

            # getting the iid
            iid: str = row[1]

            # using enumerate to iterate through each row of the file and get the index of the row
            with open(map_file_path, "r") as map_file:

                for index, map_row in enumerate(map_file):

                    map_row: list = map_row.split("\t")

                    # getting the variant id
                    variant_id: str = map_row[1]

                    # getting the genotype index positions
                    start_indx: int = index * 2

                    second_indx: int = start_indx + 1

                    # getting the genotypes from the geno_row
                    allele_1: str = geno_row[start_indx]

                    allele_2: str = geno_row[second_indx]

                    # appending to all the list
                    iid_list.append(iid)

                    variant_list.append(variant_id)

                    geno_list.append("".join([allele_1, allele_2]).strip("\n"))

    # inserting values into the genotyped dict
    genotype_dict = {
        "IID": iid_list,
        "Variant": variant_list,
        "Genotype": geno_list
    }

    genotype_df = pd.DataFrame.from_dict(genotype_dict)

    return genotype_df


def get_chr_num(file: str, pattern: str, alt_pattern: str) -> str:
    '''This function returns the chromosome number'''

    match = re.search(pattern, file)

    if match:

        chr_num: str = match.group(0)

        # removing the _ in the file name
        chr_num = chr_num[:len(chr_num)-1]

    else:

        match = re.search(alt_pattern, file)

        chr_num: str = match.group(0)

        # removing the _ in the file name
        chr_num = chr_num[:len(chr_num)-1]

    return chr_num


def search_allpair_file(allpair_file: str, carrier_list: list) -> list:
    '''This function will return a list of carriers that also share segments'''

    # loading only the columns that contain pairs from the allpair file
    allpair_df: pd.DataFrame = pd.read_csv(
        allpair_file, sep="\t", usecols=["pair_1", "pair_2"])

    # getting all the grids for pair one
    pair_1_list: list = allpair_df.pair_1.values.tolist()

    # getting all the grids for pair two into a list
    pair_2_list: list = allpair_df.pair_2.values.tolist()

    # combining the two pair list into a set so that repeating values get dropped
    pair_iid_set: set = set(pair_1_list + pair_2_list)

    # getting a list of all the carriers from carrier_list that are also in the set
    sharing_carrier_list: list = [
        grid for grid in carrier_list if grid in list(pair_iid_set)]

    # returning the list of individuals that are carriers for this variant and also share segments
    return sharing_carrier_list


def subset_genotype(geno_df: pd.DataFrame, carrier_list: list) -> pd.DataFrame:
    '''This function will suvbset the genotype df for only those who are carriers.
    This does not necessarily include only those who are confirmed carriers'''

    print(len(geno_df))
    # subsetting the dataframe to only the list of identified carriers
    geno_df_subset: pd.DataFrame = geno_df[geno_df["IID"].isin(carrier_list)]

    print(len(geno_df_subset))

    return geno_df_subset


def add_column(df_subset: pd.DataFrame, confirmed_carrier_list: list) -> pd.DataFrame:
    '''This function will add a column to the df_subset that contains either a one 
    to indicate that the grid is a carrier confirmed by shared segment or a 0 if they are not'''

    # This line just adds a column that indicates the cconfirmed status
    df_subset["confirmed_status"] = np.where(
        df_subset["IID"].isin(confirmed_carrier_list), 1, 0)

    print(df_subset)
    return df_subset


def merge_df(ideal_df: pd.DataFrame, modified_df: pd.DataFrame) -> pd.DataFrame:
    '''This function will return merge two dataframes and will return the merged dataframe'''

    # Combining the two dataframes into one without droping any rows
    final_df: pd.DataFrame = pd.concat([ideal_df, modified_df])

    return final_df


def run(args):
    "function to run"
    # Getting list of the carrier files, the map files, the ped files and the allpair_files
    carrier_files: list = get_files(args.directory, "*single_variant_list.csv")

    map_files: list = get_files(args.plink_dir, "*.map")

    ped_files: list = get_files(args.plink_dir, "*.ped")

    allpair_files: list = get_files(args.directory, "*allpair.new.txt")

    # Also need to create an empty df to be merged with others later to later
    ideal_format_df: pd.DataFrame = pd.DataFrame()

    # Iterating through the ped files
    for file in ped_files:

        # This gets the chromosome number for the files
        chr_num: str = get_chr_num(file, r"chr\d_", r"chr\d\d_")

        map_file: str = [
            map_file for map_file in map_files if chr_num in map_file][0]

        print(map_file)

        genotype_df: pd.DataFrame = form_genotype_df(map_file, file)

        # getting the correct carrier file based off of the chromosome
        car_file: str = [carrier_file for carrier_file in carrier_files if "".join(
            [chr_num, "_"]) in carrier_file][0]

        print(car_file)

        # load the car_file into a dataframe
        with open(car_file, "r") as carrier_file:

            for row in carrier_file:

                # getting the variant id
                variant: str = row.split(" ")[0]

                # getting the list of carriers
                carrier_str = row.split(" ")[1]

                # stripping characters away to get a list of grids that carry the variant
                carrier_list: list = carrier_str.strip("[]").replace(
                    "'", "").replace(",", "").split(" ")

                # using list comprehension to get the allpair file for a specific chromosome and variant
                allpair_file: str = [file for file in allpair_files if "".join(
                    [chr_num, "."]) in file and variant in file][0]

                # getting a list of carriers that share segments and therefore are "confirmed"
                confirmed_carrier_list: list = search_allpair_file(
                    allpair_file, carrier_list)

                # TODO:Need to create a function that will subset genotype_df for the whole list of carriers and then it will add a 1 to
                # the iids that share a segment and a 0 to those that don't. Then need to combine this dataframe with others iteratively

                # This will subset the dataframe for the carriers
                subset_df: pd.DataFrame = subset_genotype(
                    genotype_df, carrier_list)

                # Adding a column for the carrier status to the file
                modified_geno_df: pd.DataFrame = add_column(
                    subset_df, confirmed_carrier_list)
                print("this is the modified df")

                print(modified_geno_df)

                # Combining the dataframes
                ideal_format_df = merge_df(ideal_format_df, modified_geno_df)

                print(ideal_format_df)

    # writing the ideal_format_df to a csv file
    output_path = "".join([args.output, "confirmed_carriers.txt"])
    ideal_format_df.to_csv(output_path, sep="\t")


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument("--input", help="This argument takes the initial input file which contains variants for multiple chromosomes and splits it into multiple files, one for each chromosome",
                        dest="", type="", nargs="", required="Bool")

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
