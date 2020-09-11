from itertools import count
import pandas as pd
import glob
import os
import re
import argparse
from population_filter import Pop_Filter
# need to gather all of the single var list
cur_dir = os.getcwd()


def get_var_files(directory: str, file_id: str) -> list:

    cur_dir = os.getcwd()

    os.chdir(directory)

    file_list = []

    for file in glob.glob(file_id):

        full_file_path = "".join([directory, file])

        file_list.append(full_file_path)

    os.chdir(cur_dir)

    return file_list


def get_chr_num(carrier_file: str) -> str:

    match = re.search(r'chr\d_', carrier_file)

    if match:
        chr_num = match.group(0)

        # removing the _ in the file name
        chr_num = chr_num[:len(chr_num)-1]

    else:
        match = re.search(r'chr\d\d_', carrier_file)

        chr_num = match.group(0)

        # removing the _ in the file name
        chr_num = chr_num[:len(chr_num)-1]

    return chr_num


def get_allele_frq(carrier_file_list: list, raw_file_list: list, pop_info_filepath: str, pop_code: str, output_path: str):
    with open("".join([output_path, "allele_frequencies.txt"]), "a+") as myFile:
        myFile.write("chr\tvariant_id\tallele_freq\n")
        for file in carrier_file_list:

            chr_num = get_chr_num(file)
            chr_num = "".join([chr_num, "_"])

            car_raw_file_tuple = [(carrier_file, raw_file)
                                  for carrier_file in carrier_file_list for raw_file in raw_file_list if chr_num in carrier_file and chr_num in raw_file]

            carrier_file = car_raw_file_tuple[0][0]

            raw_file = car_raw_file_tuple[0][1]

            # getting all the variants in a list
            carrier_file_df = pd.read_csv(carrier_file, sep=",", header=None, names=[
                "variant_id", "iid_list"])

            variant_list = carrier_file_df.variant_id.values.tolist()

            # loading in the raw file into a dataframe and filtering it just for the desired population code using the Pop_Filter code
            pop_filter = Pop_Filter(pop_info_filepath, raw_file)

            pop_info_df, recode_df = pop_filter.load_files()

            pop_info_df = pop_filter.get_pop_info_subset(pop_info_df, pop_code)

            raw_file_df = pop_filter.filter_recode_df(pop_info_df, recode_df)

            print(raw_file_df.columns.tolist())
            for variant in variant_list:
                print(variant)
                # get all rows for that variant in the raw-file_df

                total_allele_count = len(raw_file_df)*2
                print(total_allele_count)

                carry_allele_df = raw_file_df[raw_file_df[variant].isin(
                    [1, 2])][variant]

                minor_allele_count = carry_allele_df.count()

                allele_frq = minor_allele_count/total_allele_count
                print(allele_frq)
                print(chr_num[:len(chr_num)-1])
                myFile.write(
                    f"{chr_num[:len(chr_num)-1]}\t{variant}\t{allele_frq}\n")

        myFile.close()


def run(args):
    "function to run"
    carrier_file_list = get_var_files(
        args.car_dir, "*.single_variant_list.csv")

    raw_file_list = get_var_files(args.raw_dir, "*.raw")

    get_allele_frq(carrier_file_list, raw_file_list,
                   args.pop_file, args.pop_code, args.output)


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument("-c", help="This argument takes the directory for all the .single_var_list.csv",
                        dest="car_dir", type=str, required=True)

    parser.add_argument("-r", help="This argument takes the directory for all the raw files",
                        dest="raw_dir", type=str, required=True)

    parser.add_argument("-i", help="This flag gives the path to a populationn info file that can filter down the raw files for a certain population",
                        dest="pop_file", type=str, required=True)

    parser.add_argument("-p", help="This flag gives the desired population code",
                        dest="pop_code", type=str, required=True)

    parser.add_argument("-o", help="This flag gives the output file path",
                        dest="output", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
