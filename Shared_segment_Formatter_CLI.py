#!/usr/bin/python

# THese are modules used
import re
import argparse
import pandas as pd
import sys
import os.path
import shutil
import multiprocessing as mp

from IBDinput_byBP_classes import Shared_Segment_Convert, Pre_Shared_Segment_Converter

####################################################################################################

# Initializing the parser


def run(args):
    # Instantiating a class object to be raised to break the inner loop

    print("running")
    ###########################################################
    # This first section will be used to get the shared segment files for each chromosome

    if os.path.exists("".join([args.output, "no_carriers_in_file.txt"])):
        print("removing the no_carriers_in_network.txt file...")
        os.remove("".join([args.output, "no_carriers_in_file.txt"]))

    preformater = Pre_Shared_Segment_Converter(
        args.input, args.pheno, args.format, args.output, args.map_file)

    segment_file_list = preformater.gather_segment_files(args.suffix)
    # also need to get all of the possible variant files per chromosome
    chr_var_file_list = preformater.gather_chromosome_files()

    # getting the list of map files
    map_file_list = preformater.get_map_files()

    # Need to make sure that the proper segment file is passed with the proper chromosome file
    for chromo_file in chr_var_file_list:

        print("In outer loop")

        match = re.search(r'chr\d_', chromo_file)

        if match:

            chr_num = match.group(0)

            # removing the _ in the file name
            chr_num = chr_num[:len(chr_num)-1]

            # adding a .
            chr_num = "".join([chr_num, "."])

        else:

            match = re.search(r'chr\d\d_', chromo_file)

            chr_num = match.group(0)

            # removing the _ in the file name
            chr_num = chr_num[:len(chr_num)-1]

            # adding a .
            chr_num = "".join([chr_num, "."])

        # dropping the 0 to try and match the naming in the segment files
        segment_chr_pattern = "".join(chr_num.split("0"))

        # Creating a tuple that gets the proper segment_file and the proper map file that corresponds to that chr

        segment_map_tuple = [(segment_file, map_file)
                             for segment_file in segment_file_list for map_file in map_file_list if chr_num in segment_file and match.group(0) in map_file]

        # If the shared segment file chr is not "chr0#" it will try the format "chr#"
        if not segment_map_tuple:

            segment_map_tuple = [(segment_file, map_file)
                                 for segment_file in segment_file_list for map_file in map_file_list if segment_chr_pattern in segment_file and match.group(0) in map_file]

        # iterating through the segment_file_list to find the shared segment file for the right chromosome

        # This checks to see if the tuple is empty or not
        if segment_map_tuple:

            # This gets the values out of the tuple
            # The first element is the segment file
            segment_file = segment_map_tuple[0][0]

            # The second element is the map file
            map_file = segment_map_tuple[0][1]

            var_info_df, variant_directory = preformater.create_variant_lists(
                chromo_file, args.var_file, map_file)

            # This checks if the var_info_file_path and the variant_directory are empty string because
            # this would mean that the the chromo_file only has one variant and it has no carriers
            if variant_directory == "":
                print(
                    f"There were no carriers in the file {chromo_file}")
                continue

            iid_file_list = preformater.get_iid_files(
                variant_directory)

            variant_bp_list = var_info_df.site.values.tolist()

            variant_id_list = var_info_df.variant_id.values.tolist()

            for var_info_tuple in zip(variant_bp_list, variant_id_list):

                variant_position = int(var_info_tuple[0])

                variant_id = str(var_info_tuple[1])

                iid_file = [
                    iid_file for iid_file in iid_file_list if variant_id in iid_file][0]

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

    parser.add_argument("-mf", "--map_file", help="This argument provides the path to the different map files",
                        dest="map_file", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
