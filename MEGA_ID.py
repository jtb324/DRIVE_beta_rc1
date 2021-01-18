#!/usr/bin/env python

import logging
import argparse
import os.path
from os import path
import sys

import carrier_analysis_scripts
import create_network_scripts
import file_creator_scripts
import allele_frequency_analysis_scripts
import pre_shared_segments_analysis_scripts
import plink_initial_format_scripts
import haplotype_segments_analysis
import full_analysis


def run(args):

    MIN_CM: int = int(
        input("Please input a value for the minimum centimorgan threshold: "))
    THREADS: int = int(input(
        "Please enter the number of threads you wish to use during this process. (Bear in mind that this number will be used for all parallelized steps):"))

    # TODO: rewrite this so that the program runs from start to finish
    if args.analysis.lower() == "gene":

        # getting the output directory to be to the variants_of_interest
        # subdirectory

        plink_initial_format_scripts.split_input_and_run_plink(
            args.var_file,
            args.output,
            args.recode_options,
            args.binary_file,
            args.plink_dir,
        )

        # TODO: need to return a value for where the plink files should be

    elif args.analysis.lower() == "maf":
        print(
            f"extracting all variants below a minor allele frequency of {args.maf_filter}"
        )

        plink_runner = plink_initial_format_scripts.PLINK_Runner(
            args.binary_file,
            args.recode_options,
            args.output,
            args.maf_filter,
        )

        plink_file_path: str = plink_runner.run_PLINK_maf_filter()
        # TODO: need to return a string listing the location of the plink
        # output files
    else:
        print("no analysis method was passed to the program.")
        print(
            "the program will now assume that the user has already gathered raw files for analysis."
        )

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
        "".join([args.output, "variants_of_interest/"]),
        args.output,
        args.pop_info,
        args.pop_code,
    )

    # The above function outputs files to a subdirectory called "carrier_analysis_output"

    print(
        "determining the minor allele frequency within the provided binary file..."
    )
    allele_frequency_analysis_scripts.determine_maf(
        "".join([args.output, "carrier_analysis_output/"]
                ), "".join([args.output, "variants_of_interest/"]),
        args.pop_info, args.pop_code, args.output
    )

    print("converting the ibd output to a human readable version...")

    IBD_search_output_files: str = "".join(
        [args.output, "formated_ibd_output/"])

    for program in args.ibd_programs:

        suffix_dict: dict = {
            "ilash": ".match.gz",
            "hapibd": ".ibd.gz",
        }

        file_suffix: str = suffix_dict[program]

        pre_shared_segments_analysis_scripts.convert_ibd(
            args.ibd_files,
            "".join([args.output, "carrier_analysis_output/"],
            program,
            IBD_search_output_files,
            "".join([args.output, "variants_of_interest/"]),
            file_suffix,
            MIN_CM,
            THREADS,
        )

    print("combining segment output...")

    pre_shared_segments_analysis_scripts.combine_output(
        IBD_search_output_files, args.ibd_programs, IBD_search_output_files,
        "".join([args.output, "variants_of_interest/"]))

    pre_shared_segments_analysis_scripts.reformat_files(
        "".join([args.output, "carrier_analysis_output/"]),
        "".join([args.output, "variants_of_interest/"]),
        IBD_search_output_files
        IBD_search_output_files,
        "".join([IBD_search_output_files, "no_carriers_in_file.txt"]),
    )

    elif args.analysis == "draw_networks":

        print("generating pdf files of networks of individuals who share segments...")

        # This dictionaru keeps track of how many carriers are actually in the network. It needs to be a global variable so that it is just extended for each variant instead of recreated
        carrier_in_network_dict=dict()


        if not path.exists("".join([args.output, "networks/"])):
            os.mkdir("".join([args.output, "networks/"]))

        output="".join(["".join([args.output, "networks/"]), "network_imgs"])

        if not path.exists(output):
            os.mkdir(output)

        carrier_in_network_dict=create_network_scripts.create_networks(
            IBD_search_output_files,  "".join(
                [args.output, "variants_of_interest/"]),
            carrier_in_network_dict,  "".join([args.output, "networks/"])
        )

        # Writing the dictionary to a csv file
        csv_writer=file_creator_scripts.Csv_Writer_Object(
            carrier_in_network_dict,  "".join(
                [args.output, "networks/"]), "carriers_in_networks.csv"
        )

        csv_writer.write_to_csv()

    print("getting information about the haplotypes for the confirmed carriers")

    if not path.exists("".join([args.output, "haplotype_analysis/"]))

        os.mkdir("".join([args.output, "haplotype_analysis/"]))

    haplotype_segments_analysis.get_segment_lengths(
        "".join([args.output, "network/", "confirmed_carriers.txt"]),
        "".join([args.output, "haplotype_analysis/"])
        args.ilash_dir,
        args.hapibd_dir,
        THREADS,
        IBD_search_output_files,
         "".join([args.output, "variants_of_interest/", "reformated/"]),
        "".join([args.output, "networks/", "network_imgs/network_groups.csv"]),
        IBD_search_output_files,
    )

    print(
        "generating a file contain the genotypes for all IIDs in the provided file"
    )

    full_analysis.get_all_genotypes(
        IBD_search_output_files, "".join([args.output, "network/", "confirmed_carriers.txt"]
                                         ), args.pop_info, "".join([args.output, "haplotype_analysis/"]), args.pop_code
    )


def main():
    parser=argparse.ArgumentParser(
        description="This identifies individuals who have a specific variant in a raw file from PLINK"
    )

    parser.add_argument(
        "--vfile",
        help="This is the pathway for the input file which contains a list of all the snp ids using the mega array ids and then the chromsome numbers. This file should have at least two columns titled 'SNP' and 'Chr'(The capitalization matters.).",
        dest="variant_file",
        type=str,
    )

    parser.add_argument(
        "--rfile",
        "-r",
        help="This is the pathway for the PLINK recoded input file. ",
        dest="raw_file",
        type=str,
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
        "--cfile",
        "-c",
        help="This argument supplies the path for the carrier files that end in '.single_var_list.csv",
        dest="carrier_file",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--recode_options",
        help="TThis argument list the recode option used to run plink",
        dest="recode_options",
        nargs="+",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--ilash_files",
        help="This argument list the directory for all ilash files ending in .match.gz",
        dest="ilash_dir",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--hapibd_files",
        help="This argument list the directory for all the hapibd files ending in .ibd.gz",
        dest="hapibd_dir",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--nfiles",
        help="This argument list the path for the 'network_groups.csv file which list what network pairs are a part of'",
        dest="network_files",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--rfile",
        "-r",
        help="This argument supplies the path for the raw output files from PLINK. These files will be found in the directory ./variants_of_interest/ if the split_file was run",
        dest="raw_file",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--converted_ibd",
        help="This argument list the directory for all of the ibd files that were converted into a human readbale path. THese files should end in small.txt.gz",
        dest="segment_files",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--reformated_carriers",
        help="This argument list the directory for all of the reformated carrier files list. These files should be in a subdirected called reformat",
        dest="reformat_files",
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
        "--thread",
        "-t",
        help="This argument list how many workers the program can use for multiprocessing",
        dest="threads",
        type=str,
        required=False,
        default=1,
    )

    parser.add_argument(
        "--ifiles",
        "-ibd",
        help="This argument list the directory that contains the output from the ibd files",
        dest="ibd_files",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--plink_dir",
        "-pd",
        help="This argument list the directory that contains the plink files for the run. These should be map and ped files.",
        dest="plink_dir",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--no_carrier_file",
        "-nc",
        help="This argument list the file path to the no_carriers.txt file",
        dest="no_carriers_file",
        type=str,
        required=False,
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
        "--file_suffix",
        "-fs",
        help="This argument list the suffix for the different ibd files because it could change",
        dest="ibd_file_suffix",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--var_input",
        "-vi",
        help="This argument list the path to the input file that list all the variants. This should be a csv or xlsx file",
        dest="var_input_file",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--maf_filter",
        "-maf",
        help="This argument list the max minor allele frequency used to extract variants from teh provided binary file",
        dest="maf_filter",
        type=str,
        required=False,
        default=0.01,
    )

    parser.add_argument(
        "--drop_var",
        help="This functionality is used to drop variants from a file if needed to for some reason. This is passed into the searchPedigree function incase maybe a certain variant is too common and can be removed",
        dest="drop_var",
        type=str,
        nargs="+",
    )

    parser.add_argument(
        "--variant_file",
        help="This argument provides a path to a file that list all individuals that "
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

    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
