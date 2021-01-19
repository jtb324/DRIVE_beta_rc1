#!/usr/bin/env python

import logging
import argparse
import os.path
from os import path
import sys
from datetime import datetime

import carrier_analysis_scripts
import create_network_scripts
import file_creator_scripts
import allele_frequency_analysis_scripts
import pre_shared_segments_analysis_scripts
import plink_initial_format_scripts
import haplotype_segments_analysis
import full_analysis


def run(args):
    # Asking for user input for constants that will be used throughout the program
    ANALYSIS_TYPE: str = input("Please input an analysis type: ").strip(" ")

    if ANALYSIS_TYPE not in ["gene", "maf", ""]:
        print(
            "Please choose one of the two allowed analysis types ('gene'/'maf', or "
            ")")
        sys.exit(1)

    MIN_CM: int = int(
        input(
            "Please input a value for the minimum centimorgan threshold: (Default is 3 cM) "
        ))

    if not MIN_CM:
        MIN_CM = 3

    THREADS: int = int(
        input(
            "Please enter the number of threads you wish to use during this process. The default value is 3. (Bear in mind that this number will be used for all parallelized steps): "
        ))

    if not THREADS:
        THREADS = 3

    if ANALYSIS_TYPE == "maf":

        MAF_FILTER: int = int(
            input(
                "Please input a minor allele frequency threshold to be used. (The default is 0.05): "
            ))

    ILASH_PATH: str = "/data100t1/share/BioVU/shapeit4/Eur_70k/iLash/min100gmap/"

    HAPIBD_PATH: str = "/data100t1/share/BioVU/shapeit4/Eur_70k/hapibd/"

    IBD_PATHS_LIST: list = [ILASH_PATH, HAPIBD_PATH]
    # TO
    if ANALYSIS_TYPE.lower() == "gene":

        # getting the output directory to be to the variants_of_interest
        # subdirectory

        plink_file_path: str = plink_initial_format_scripts.split_input_and_run_plink(
            args.var_file,
            args.output,
            args.recode_options,
            args.binary_file,
            "".join([args.output, "plink_output_files/"]),
        )

        # TODO: need to return a value for where the plink files should be

    elif ANALYSIS_TYPE.lower() == "maf":
        print(
            f"extracting all variants below a minor allele frequency of {args.maf_filter}"
        )
        #TODO: adjust this so that it accept a range of variants not a numeric range
        if args.range:
            START: int = int(args.range[0])
            END: int = int(args.range[1])

            plink_runner = plink_initial_format_scripts.PLINK_Runner(
                args.recode_options,
                args.output,
                args.binary_file,
                maf_filter=MAF_FILTER,
                start=START,
                end=END,
            )
        else:
            plink_runner = plink_initial_format_scripts.PLINK_Runner(
                args.recode_options,
                args.output,
                args.binary_file,
                maf_filter=MAF_FILTER,
            )

        plink_file_path: str = plink_runner.run_PLINK_maf_filter()
        # TODO: need to return a string listing the location of the plink
        # output files
    else:
        print("no analysis method was passed to the program.")
        print(
            "the program will now assume that the user has already gathered raw files for analysis."
        )
        # Setting the plink_file_path if the top two steps are skipped
        plink_file_path = "".join([args.output, "plink_output_files/"])
    # TODO: need to adjust the section so that it takes the proper options. At this moment this feature is not being used so it is commented out, but it is worth keeping in the program
    # if args.analysis == "multi":

    #     print("generating list of individuals carrying multiple variants....")

    #     carrier_analysis_scripts.multiVariantAnalysis(
    #         args.input, args.output, args.compatible_format, "multi_variant_list.csv"
    #     )

    #     logging.info(
    #         "Finished creating a list of individuals carrying multiple variants"
    #     )

    print("generating list of individuals at each probe id...")

    # The args.input should be a directory indicating where the raw files are located
    carrier_analysis_scripts.singleVariantAnalysis(
        plink_file_path,
        args.output,
        args.pop_info,
        args.pop_code,
    )

    # The above function outputs files to a subdirectory called "carrier_analysis_output"

    print(
        "determining the minor allele frequency within the provided binary file..."
    )
    allele_frequency_analysis_scripts.determine_maf(
        "".join([args.output, "carrier_analysis_output/"]),
        "".join([args.output, "plink_output_files/"]), args.pop_info,
        args.pop_code, args.output)

    print("converting the ibd output to a human readable version...")

    IBD_search_output_files: str = "".join(
        [args.output, "formated_ibd_output/"])

    if not path.exists(IBD_search_output_files):

        os.mkdir(IBD_search_output_files)

    for program in args.ibd_programs:

        suffix_dict: dict = {
            "ilash": ".match.gz",
            "hapibd": ".ibd.gz",
        }

        file_suffix: str = suffix_dict[program]

        # getting the correct ibd_file_path
        ibd_file: str = [
            file for file in IBD_PATHS_LIST if program in file.lower()
        ][0]

        print(ibd_file)

        pre_shared_segments_analysis_scripts.convert_ibd(
            ibd_file, "".join([args.output, "carrier_analysis_output/"]),
            program, IBD_search_output_files,
            "".join([args.output,
                     "plink_output_files/"]), file_suffix, MIN_CM, THREADS)
    print("combining segment output...")

    pre_shared_segments_analysis_scripts.combine_output(
        IBD_search_output_files, args.ibd_programs, IBD_search_output_files,
        "".join([args.output, "carrier_analysis_output/reformated/"]))

    pre_shared_segments_analysis_scripts.reformat_files(
        "".join([args.output, "carrier_analysis_output/"]),
        "".join([args.output, "plink_output_files/"]),
        IBD_search_output_files,
        IBD_search_output_files,
        "".join([IBD_search_output_files, "no_carriers_in_file.txt"]),
    )

    print(
        "generating pdf files of networks of individuals who share segments..."
    )

    # This dictionaru keeps track of how many carriers are actually in the network. It needs to be a global variable so that it is just extended for each variant instead of recreated
    carrier_in_network_dict = dict()

    if not path.exists("".join([args.output, "networks/"])):
        os.mkdir("".join([args.output, "networks/"]))

    output = "".join(["".join([args.output, "networks/"]), "network_imgs"])

    if not path.exists(output):
        os.mkdir(output)

    carrier_in_network_dict = create_network_scripts.create_networks(
        IBD_search_output_files,
        "".join([args.output,
                 "carrier_analysis_output/"]), carrier_in_network_dict, output)

    # Writing the dictionary to a csv file
    csv_writer = file_creator_scripts.Csv_Writer_Object(
        carrier_in_network_dict, "".join([args.output, "networks/"]),
        "carriers_in_networks.csv")
    
    csv_writer.write_to_csv()

    print(
        "getting information about the haplotypes for the confirmed carriers")

    if not path.exists("".join([args.output, "haplotype_analysis/"])):

        os.mkdir("".join([args.output, "haplotype_analysis/"]))

    haplotype_segments_analysis.get_segment_lengths(
        "".join([IBD_search_output_files, "confirmed_carriers.txt"]),
        "".join([args.output, "haplotype_analysis/"]),
        ILASH_PATH,
        HAPIBD_PATH,
        THREADS,
        "".join([args.output, "plink_output_files/"]),
        "".join([args.output, "carrier_analysis_output/", "reformated/"]),
        "".join([args.output, "networks/", "network_imgs/network_groups.csv"]),
        IBD_search_output_files,
    )

    print(
        "generating a file contain the genotypes for all IIDs in the provided file"
    )

    full_analysis.get_all_genotypes(
        "".join([args.output, "plink_output_files/"]),
        "".join([IBD_search_output_files, "confirmed_carriers.txt"]),
        args.pop_info, "".join([args.output,
                                "haplotype_analysis/"]), args.pop_code)

    # TODO": Fix the timing issue so that it gives the correct time
    finishing_time = datetime.utcnow()
    print(
        f"The program successfully finished at {finishing_time.strftime('%H:%M:%S')}"
    )


def main():
    parser = argparse.ArgumentParser(
        description="This identifies individuals who have a specific variant in a raw file from PLINK"
    )

    parser.add_argument(
        "--bfile",
        "-b",
        help="This argument will list the directory to the bim file which give them"
        "genotype information that PLINK uses",
        dest="binary_file",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--recode_options",
        help="This argument list the recode option used to run plink",
        dest="recode_options",
        nargs="+",
        type=str,
        required=True,


    parser.add_argument(
        "--maf_filter",
        help="This argument list the maximum minor allele frequency to be used (ex: 0.01)",
        dest="maf_filter",
        type=str,
        required=False,
    )


    parser.add_argument(
        "--output",
        "-o",
        help="This is the directory that text files containing the ids of the individuals who have desired variants will be written to.",
        dest="output",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--analysis",
        help="This tag indicates that the multiVariantAnalysis function will be called "
        "to analyze how many individuals carry multiple variants. Two csv files are "
        "made which contain the indices of the variants and a list of the "
        "individuals that contain those variants. This accepts single, total, multi"
        ", matchPED, allele_counts, draw_networks",
        dest="analysis",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--ibd_programs",
        help="This argument list which ibd programs were used",
        dest="ibd_programs",
        type=str,
        nargs="+",
        required=False,
    )

    parser.add_argument(
        "--variant_file",
        help=
        "This argument provides a path to a file that list all individuals that "
        "carry a specific variant",
        dest="var_file",
        type=str,
        default=False,
    )

    parser.add_argument(
        "--pop_info",
        help="This argument provides the file path to a file containing the population "
        "distribution of a dataset for each grid. This file should be a text file "
        "and at least contain two columns, where one column is 'Pop', the "
        "population code for each grid based on 1000 genomes, and then the second "
        "column is 'grid', which list the IIDs for each each grid.",
        dest="pop_info",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--pop_code",
        help="This argument can be a single population code or a list of population "
        "codes that someone is looking for. Population codes have to match the "
        "1000 genomes.",
        dest="pop_code",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--range",
        help=
        "This argument will list the  start and end bp of the range you wish to look at. The argument should be formated like '--range START END'.",
        dest="range",
        nargs="+",
        type=list,
        required=False,
    )

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
