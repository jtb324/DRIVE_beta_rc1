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


def subset_genotype(geno_df: pd.DataFrame, carrier_list: list, variant_id: str) -> pd.DataFrame:
    '''This function will suvbset the genotype df for only those who are carriers.
    This does not necessarily include only those who are confirmed carriers'''

    # subsetting the dataframe to only the list of identified carriers
    geno_df_subset: pd.DataFrame = geno_df[geno_df["IID"].isin(
        carrier_list)]

    geno_df_subset = geno_df_subset[geno_df_subset.Variant == variant_id]

    return geno_df_subset


def add_column(df_subset: pd.DataFrame, confirmed_carrier_list: list) -> pd.DataFrame:
    '''This function will add a column to the df_subset that contains either a one
    to indicate that the grid is a carrier confirmed by shared segment or a 0 if they are not'''

    # This line just adds a column that indicates the cconfirmed status
    df_subset["confirmed_status"] = np.where(
        df_subset["IID"].isin(confirmed_carrier_list), 1, 0)

    return df_subset


def write_to_file(df: pd.DataFrame, output_path):
    '''This function will write each row of the df to a txt file. This has to be done to avoid memory issues'''

    # writing the dataframe to a file
    df.to_csv(output_path, header=None, sep="\t", index=None, mode="a")


def check_no_carrier(no_carrier_file: str, variant_id: str) -> int:
    '''This function will check if the variant is in a file called no_carriers_in_file.txt'''

    # if this file is not present than the function needs to return 0
    if os.path.isfile(no_carrier_file):

        return 0

    else:
        # loading the file into a dataframe
        no_carrier_df: pd.DataFrame = pd.read_csv(
            no_carrier_file, sep="\t", names=["variant", "chr"])

        # getting a list of all variants that had no carriers
        variant_list: list = no_carrier_df.variant.values.tolist()

        # figuring out if the current id is not in the no_carrier list
        if variant_id in variant_list:

            # returns 1 if it is
            return 1

        else:

            print(f"The variant {variant_id} is not a valid variants")

            # returns 0 if the variant is not found
            return 0


def run(args):
    "function to run"

    # defining an output path and file name for the output file
    output_path: str = "".join([args.output, "confirmed_carriers.txt"])

    # checking if the file exist
    if os.path.isfile(output_path):

        # removing the file if it exist
        os.remove(output_path)

    # writing the header line to the output text file
    with open(output_path, "w") as output_file:

        output_file.write(
            f"{'IID'}\t{'variant_id'}\t{'genotype'}\t{'confirmed_status'}\n")

    # Getting list of the carrier files, the map files, the ped files and the allpair_files
    carrier_files: list = get_files(args.directory, "*single_variant_list.csv")

    map_files: list = get_files(args.plink_dir, "*.map")

    ped_files: list = get_files(args.plink_dir, "*.ped")

    allpair_files: list = get_files(args.allpair_dir, "*allpair.txt")

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
                variant: str = row.split(",", maxsplit=1)[0]

                # getting the list of carriers
                carrier_str: str = row.split(",", maxsplit=1)[1]

                # stripping characters away to get a list of grids that carry the variant
                carrier_str = carrier_str.replace('"', '')
                carrier_str = carrier_str.strip("[]")

                carrier_str = carrier_str.replace(
                    "'", "")
                carrier_str = carrier_str.replace(']', '')
                # carrier_str = carrier_str.strip("")

                carrier_str = carrier_str.replace(
                    ",", "")

                carrier_str = carrier_str.replace(
                    "\n", "")

                carrier_list: list = carrier_str.split(" ")

                # using list comprehension to get the allpair file for a specific chromosome and variant
                print(f"{chr_num},{variant}")

                # There is an error if it does not find an allpair_file. Some of these files don't exist because there are no carriers
                try:
                    allpair_file: str = [file for file in allpair_files if "".join(
                        [chr_num, "."]) in file and variant in file][0]

                except IndexError:

                    # returns either a 1 or 0 if the variant is in the no_carrier_file.txt list
                    carrier_int: int = check_no_carrier(
                        args.no_carrier_file, variant)

                    if carrier_int == 1:

                        # print statement that explains that ther is no variant
                        print(f"The variant, {variant}, has no carriers")
                        # skips to the next iteration of the for loop
                        continue

                    elif carrier_int == 0:

                        print(f"The variant, {variant}, failed")

                        # checking if the file exist
                        if os.path.isfile("".join([args.output, "failed_variants.txt"])):

                            # removing the file if it exist
                            os.remove(
                                "".join([args.output, "failed_variants.txt"]))

                        # writing the variant that failed to a file
                        with open("".join([args.output, "failed_variants.txt"]), "a+") as file:

                            file.write(f"{variant}\n")

                        continue

                print(allpair_file)
                # getting a list of carriers that share segments and therefore are "confirmed"

                confirmed_carrier_list: list = search_allpair_file(
                    allpair_file, carrier_list)

                # This will subset the dataframe for the carriers
                subset_df: pd.DataFrame = subset_genotype(
                    genotype_df, carrier_list, variant[:-2])

                # Adding a column for the carrier status to the file
                modified_geno_df: pd.DataFrame = add_column(
                    subset_df, confirmed_carrier_list)

                print(modified_geno_df)

                # Combining the dataframes
                write_to_file(modified_geno_df, output_path)


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument("-d", help="This argument list the path for the carrier csv files. These files should end with single_variant_list.csv",
                        dest="directory", type=str, required=True)

    parser.add_argument("-p", help="This argument list the path for the directory containing the PLINK output files. These files will be the map and ped files",
                        dest="plink_dir", type=str, required=True)

    parser.add_argument("-a", help="This argument list the path for the allpair.txt files.",
                        dest="allpair_dir", type=str, required=True)

    parser.add_argument("-o", help="This argument list the path for the output file",
                        dest="output", type=str, required=True)

    parser.add_argument("-n", help="This argument list the path for the no_carrier_in_file.txt",
                        dest="no_carrier_file", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
