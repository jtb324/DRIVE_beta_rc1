#!/usr/bin/env python

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
import utility_scripts


def run(args):
    # creating the README for the main parent directory
    readme = utility_scripts.Readme("_README.md", args.output)

    readme.rm_previous_file()

    readme.write_header(args.output)

    readme.create_date_info()

    readme.add_line(utility_scripts.main_parameter_text)
    readme.add_line(utility_scripts.main_directory_header)
    readme.add_line(utility_scripts.main_directory_text)

    # Next few lines give settings for the logger

    # Creating a logfile

    file_name: str = "".join([args.output, "run.log"])

    if path.exists(file_name):
        os.remove(file_name)

    logfile = utility_scripts.LogFile("run.log", args.output)
    logfile.write_header()
    logfile.create_date_info()
    logfile.add_newline("INFO", "Starting the run...\n")

    # Logging some of the info about the users input
    logfile.add_newline("INFO", f"Binary File: {args.binary_file}\n")
    logfile.add_newline("INFO", f"Recode options: {args.recode_options}\n")
    logfile.add_newline("INFO", f"Writing the output to: {args.output}\n")
    logfile.add_newline(
        "INFO", f"The ibd programs being used are: {args.ibd_programs}\n")
    logfile.add_newline(
        "INFO",
        f"The population information file describing the demographics is found at: {args.pop_info}\n"
    )
    logfile.add_newline(
        "INFO", f"The population code being used is: {args.pop_code}\n")

    # Asking for user input for constants that will be used throughout the program
    ANALYSIS_TYPE: str = utility_scripts.ask_for_analysis_type(logfile)

    MIN_CM: str = utility_scripts.ask_for_min_cm(logfile)

    THREADS: str = utility_scripts.ask_for_thread_count(logfile)

    MAF_FILTER: str = utility_scripts.ask_for_maf_filter(logfile)

    ILASH_PATH: str = "/data100t1/share/BioVU/shapeit4/Eur_70k/iLash/min100gmap/"

    HAPIBD_PATH: str = "/data100t1/share/BioVU/shapeit4/Eur_70k/hapibd/"

    logfile.add_newline(
        "INFO", f"The specified path to the iLASH files are {ILASH_PATH}\n")

    logfile.add_newline(
        "INFO",
        f"The specified path to the hapibd files are {HAPIBD_PATH}\n\n")

    IBD_PATHS_LIST: list = [ILASH_PATH, HAPIBD_PATH]
    # TO
    # Next line checks to see if a plink output file already exist and delets it if it does
    if os.path.exists("".join(
        [args.output, "plink_output_files/", "plink_log.log"])):

        os.remove("".join(
            [args.output, "plink_output_files/", "plink_log.log"]))

    if ANALYSIS_TYPE.lower() == "gene":

        # getting the output directory to be to the variants_of_interest
        # subdirectory
        logfile.add_newline(
            "INFO",
            f"generating a list of snps from the file {args.var_file}\n")

        plink_file_path: str = plink_initial_format_scripts.split_input_and_run_plink(
            args.var_file, args.output, args.recode_options, args.binary_file,
            "".join([args.output, "plink_output_files/"]), MAF_FILTER)

        missing_var_count: int = plink_initial_format_scripts.check_for_missing_var(
            plink_file_path, args.var_file)

        logfile.add_newline(
            "INFO",
            f"There were {missing_var_count} missing variants between the provided input file at {args.var_file} and the plink output in the directory {plink_file_path}\n"
        )
    # This next section runs the first level of the program using the maf analysis type
    # this means that the
    elif ANALYSIS_TYPE.lower() == "maf":

        CHR: str = input(
            "Please input the chromosome that the variant of interest is on. Please use a leading 0 for single digit numbers: "
        )

        logfile.add_newline(
            "INFO",
            f"Setting the chromosome of interest to be chromosome {CHR}\n")
        print(
            f"extracting all variants below a minor allele frequency of {MAF_FILTER}"
        )

        if args.range:
            START_RS: str = args.range[0]
            END_RS: str = args.range[1]

            logfile.add_newline(
                "INFO",
                f"beginning analysis for the range starting with the variant {START_RS} and ending at {END_RS}\n"
            )

            plink_runner = plink_initial_format_scripts.PLINK_Runner(
                args.recode_options,
                args.output,
                args.binary_file,
                maf_filter=MAF_FILTER,
                start_rs=START_RS,
                end_rs=END_RS,
                chr_num=CHR,
            )
        else:
            plink_runner = plink_initial_format_scripts.PLINK_Runner(
                args.recode_options,
                args.output,
                args.binary_file,
                maf_filter=MAF_FILTER,
            )

        plink_file_path: str = plink_runner.run_PLINK_maf_filter()

        logfile.add_newline(
            "INFO", f"PLINK output files written to: {plink_file_path}\n")

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

    #     logfile.add_newline(
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

    logfile.add_newline(
        "INFO",
        f"Writing the results of the carrier analysis called: {''.join([args.output, 'carrier_analysis_output/'])}\n"
    )
    # The above function outputs files to a subdirectory called
    # "carrier_analysis_output"

    print(
        "determining the minor allele frequency within the provided binary file..."
    )
    allele_frequency_analysis_scripts.determine_maf(
        "".join([args.output, "carrier_analysis_output/"]),
        "".join([args.output, "plink_output_files/"]), args.pop_info,
        args.pop_code, args.output)

    THRESHOLD: float = 0.10

    print(
        f"checking the variant minor allele frequencies against an arbitrary threshold of {THRESHOLD}"
    )
    # This next function will check to see if variants are above a specified threshold
    # If they are then the user has an option to end the program and remove the variants
    # or just continue on with the program.
    # The function will return a tuple where the first value is a list containing variants
    # that are above a specified threshold and the second value is either 0 or 1
    # where 0 quits the program and 1 continues

    variants_above_threshold_tuple: tuple = allele_frequency_analysis_scripts.check_mafs(
        "".join([
            args.output, "carrier_analysis_output/", "allele_frequencies.txt"
        ]), THRESHOLD)
    variants_above_threshold: list = variants_above_threshold_tuple[0]

    program_end_code: int = variants_above_threshold_tuple[1]

    if program_end_code == 0:

        logfile.add_newline(
            "WARNING",
            f"The variants {', '.join(variants_above_threshold)} were above the arbitrary threshold of {THRESHOLD}\n"
        )

        logfile.add_newline(
            "INFO",
            "ending program so that user can remove the above variants\n")
        sys.exit(1)

    elif program_end_code == 1 and variants_above_threshold:

        logfile.add_newline(
            "WARNING",
            f"The variants {', '.join(variants_above_threshold)} were above the arbitrary threshold of {THRESHOLD}\n"
        )

    print("converting the ibd output to a human readable version...")

    IBD_search_output_files: str = "".join(
        [args.output, "formatted_ibd_output/"])

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

        pre_shared_segments_analysis_scripts.convert_ibd(
            ibd_file, "".join([args.output, "carrier_analysis_output/"]),
            program, IBD_search_output_files,
            "".join([args.output,
                     "plink_output_files/"]), file_suffix, MIN_CM, THREADS)
    print("combining segment output...")

    pre_shared_segments_analysis_scripts.combine_output(
        "".join([IBD_search_output_files, "reformatted_ibd_output/"]),
        args.ibd_programs, IBD_search_output_files,
        "".join([args.output, "carrier_analysis_output/reformatted/"]))

    pre_shared_segments_analysis_scripts.reformat_files(
        "".join([args.output, "carrier_analysis_output/"]),
        "".join([args.output, "plink_output_files/"]),
        IBD_search_output_files,
        IBD_search_output_files,
        "".join([IBD_search_output_files, "no_carriers_in_file.txt"]),
    )
    logfile.add_newline(
        "INFO",
        f"Writing the results from the IBD conversion files to: {IBD_search_output_files}\n"
    )

    print(
        "generating pdf files of networks of individuals who share segments..."
    )

    # This dictionaru keeps track of how many carriers are actually in the network. It needs to be a global variable so that it is just extended for each variant instead of recreated

    if not path.exists("".join([args.output, "networks/"])):
        os.mkdir("".join([args.output, "networks/"]))

    # output = "".join(["".join([args.output, "networks/"]), "network_imgs"])

    network_dir: str = "".join([args.output, "networks/"])

    create_network_scripts.create_networks(
        IBD_search_output_files,
        "".join([args.output,
                 "carrier_analysis_output/"]), network_dir, network_dir)

    logfile.add_newline(
        "INFO",
        f"Writing the results of the network analysis to: {''.join([args.output, 'networks/'])}\n"
    )

    print(
        "getting information about the haplotypes for the confirmed carriers")

    if not path.exists("".join([args.output, "haplotype_analysis/"])):

        os.mkdir("".join([args.output, "haplotype_analysis/"]))

    haplotype_info_path: str = haplotype_segments_analysis.get_segment_lengths(
        "".join([IBD_search_output_files, "confirmed_carriers.txt"]),
        "".join([args.output, "haplotype_analysis/"]),
        ILASH_PATH,
        HAPIBD_PATH,
        THREADS,
        "".join([args.output, "plink_output_files/"]),
        "".join([args.output, "carrier_analysis_output/", "reformatted/"]),
        "".join([args.output, "networks/", "network_groups.csv"]),
        IBD_search_output_files,
    )

    logfile.add_newline(
        "INFO",
        f"Writing the output of the haplotype analysis to: {''.join([args.output, 'haplotype_analysis/'])}\n"
    )
    print(
        "generating a file contain the genotypes for all IIDs in the provided file"
    )

    full_analysis.get_all_genotypes(
        "".join([args.output, "plink_output_files/"]),
        "".join([IBD_search_output_files, "confirmed_carriers.txt"]),
        args.pop_info, "".join([args.output,
                                "haplotype_analysis/"]), args.pop_code)

    logfile.add_newline(
        "INFO",
        f"Writing the result of getting all the genotypes for all IIDs in the provided file to: {''.join([args.output, 'haplotype_analysis/'])}\n"
    )

    # Add the function that gets the haplotype string to this
    logfile.add_newline(
        "INFO",
        f"Writing the most probable haplotypes for each network to the file 'network_haplotypes.txt' at {''.join([args.output, 'haplotype_analysis/'])}\n"
    )
    print("Finding the most probable haplotypes")

    haplotype_segments_analysis.gather_haplotypes(
        haplotype_info_path, "".join([args.output, "haplotype_analysis/"]),
        args.binary_file, args.pop_info, args.pop_code)

    logfile.add_newline("INFO", 'Analysis finished...\n')

    # TODO": Fix the timing issue so that it gives the correct time
    finishing_time = datetime.utcnow()
    print(
        f"The program successfully finished at {finishing_time.strftime('%H:%M:%S')}"
    )


def main():
    parser = argparse.ArgumentParser(
        description=
        "This identifies individuals who have a specific variant in a raw file from PLINK"
    )

    parser.add_argument(
        "--bfile",
        "-b",
        help=
        "This argument will list the directory to the bim file which give them"
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
    )

    parser.add_argument(
        "--output",
        "-o",
        help=
        "This is the directory that text files containing the ids of the individuals who have desired variants will be written to.",
        dest="output",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--analysis",
        help=
        "This flag will pass the analysis type if the user wants to first run plink for the variants. The flag expects the argument to either be 'gene' or 'maf'. If the maf option is chosen than the user also needs to specify a start and end bp range and a chromosome as well as a minor allele frequency threshold",
        dest="analysis",
        type=str,
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
        help=
        "This argument provides the file path to a file containing the population "
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
        help=
        "This argument can be a single population code or a list of population "
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
        type=str,
        required=False,
    )

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
