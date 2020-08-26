#!/usr/bin/python

# THese are modules used
import re
import argparse
import pandas as pd
import sys
import os
import shutil

from IBDinput_byBP_classes import Shared_Segment_Convert, Pre_Shared_Segment_Converter

####################################################################################################
# Initializing the parser


def run(args):
    print("running")
    ###########################################################
    # This first section will be used to get the shared segment files for each chromosome

    preformater = Pre_Shared_Segment_Converter(
        args.input, args.pheno, args.format, args.output)

    segment_file_list = preformater.gather_segment_files(args.suffix)

    # also need to get all of the possible variant files per chromosome
    chr_var_file_list = preformater.gather_chromosome_files()
    # Need to make sure that the proper segment file is passed with the proper chromosome file

    for chromo_file in chr_var_file_list:
        # determining what chromosome is being evaluated
        print(chromo_file)
        match = re.search(r'chr\d*', chromo_file)

        print(match)

        chr_num = match.group(0)

        # iterating through the segment_file_list to find the shared segment file for the right chromosome
        for segment_file in segment_file_list:

            print(segment_file)

            # checking if the chr_num matches
            if chr_num in segment_file:
                print("match found")

                var_info_file_path, variant_directory = preformater.create_variant_lists(
                    chromo_file, args.var_file)

                iid_file_list = preformater.get_iid_files(variant_directory)
        # reading the variant baseposition txt file into a dataframe
                print(var_info_file_path)
                variant_bp_df = pd.read_csv(var_info_file_path, header=None, names=[
                    "output_file_name", "variant_bp", "variant_id"], sep="\t")

                print(variant_bp_df)

                variant_bp_list = variant_bp_df.variant_bp.values.tolist()
                print(variant_bp_list)

                variant_id_list = variant_bp_df.variant_id.values.tolist()
                print(variant_id_list)

                for var_info_tuple in zip(variant_bp_list, variant_id_list):
                    print(var_info_tuple)
                    variant_position = int(var_info_tuple[0])

                    variant_id = str(var_info_tuple[1])

                    for iid_file in iid_file_list:

                        if variant_id in iid_file:

                            pheno_file = iid_file

                            ibd_file_converter = Shared_Segment_Convert(
                                segment_file, pheno_file, args.output, args.format, args.min, args.thread, variant_position, variant_id)

                            parameter_dict = ibd_file_converter.generate_parameters()

                            uniqID, dupID = ibd_file_converter.build_id_pairs()

                            IBDdata, IBDindex = ibd_file_converter.create_ibd_arrays()

                            ibd_file_converter.run_parallel(
                                IBDdata, IBDindex, parameter_dict, uniqID)

        # delete the variant_list directory for the next iteration
        # Also delete the variant_info.txt file because this is only good for the selected values


def main():
    parser = argparse.ArgumentParser(
        description="This script converts the output of the  IBD detection programs to a human readable form")

    parser.add_argument("-i", '--input', help="This argument just list the path to the IBD detection software output",
                        dest="input", type=str, required=True)

    parser.add_argument("-p", "--pheno", help="This argument just list the directory to the single_variant_list.csv files",
                        dest="pheno", type=str, required=True)

    parser.add_argument("-o", "--output", help="This argument list the output directory",
                        dest="output", type=str, required=True)

    parser.add_argument("-f", "--format", help="This argument specifies the IBD program used",
                        dest="format",  type=str, required=True)

    parser.add_argument("-m", "--min", help="This argument specifies the minimum cM threshold",
                        dest="min", type=str, required=True)

    parser.add_argument("-t", "--thread", help="This argument specifies the threadcount to be used in parallel",
                        dest="thread", type=str, required=True)

    parser.add_argument("-s", "--suffix", help="This argument specifies the file suffix used because files may be named different for different IBD programs used",
                        dest="suffix", type=str, required=True)

    parser.add_argument("-v", "--var", help="This argument provides the path to the original variant file. It is used to get the base position of the variants",
                        dest="var_file", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
