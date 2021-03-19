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
    # TODO: refactor these readme section
    readme = utility_scripts.Readme("_README.md", args.output)

    readme.rm_previous_file()

    readme.write_header(args.output)

    readme.create_date_info()

    readme.add_line(utility_scripts.main_parameter_text)
    readme.add_line(utility_scripts.main_directory_header)
    readme.add_line(utility_scripts.main_directory_text)

    # Next few lines give settings for the logger

    # Creating a logfile with the current day as part of the log file name
    current_day = datetime.now().strftime("%m_%d_%Y")

    file_name: str = "".join([args.output, current_day, "_run.log"])

    logger: object = utility_scripts.create_logger(__name__, file_name)

    utility_scripts.record_user_arguments(logger, args)

    # creating an object that will ask for user input for the Analysis type,
    # the minimum  centimorgan threshold, the threads to be used and the
    # minor allele frequency threshold
    input_gatherer: object = utility_scripts.Input_Gather()

    # Getting the objects attributes
    object_dict: dict = utility_scripts.get_dict_of_variables(input_gatherer)

    # unpacking the dictionary into the user parameters
    ANALYSIS_TYPE: str = object_dict["ANALYSIS_TYPE"]
    MIN_CM: str = object_dict["MIN_CM"]
    THREADS: str = object_dict["THREADS"]
    MAF_FILTER: str = object_dict["MAF_THRESHOLD"]

    # setting parameters for the ilash file paths and the hapibd file paths
    ILASH_PATH: str = "/data100t1/share/BioVU/shapeit4/Eur_70k/iLash/min100gmap/"

    HAPIBD_PATH: str = "/data100t1/share/BioVU/shapeit4/Eur_70k/hapibd/"

    logger.info(f"iLASH files directory: {ILASH_PATH}")
    logger.info(f"Hapibd files directory: {HAPIBD_PATH}")

    IBD_PATHS_LIST: list = [ILASH_PATH, HAPIBD_PATH]

    # output path that all the plink files will output to
    plink_file_path = "".join([args.output, "plink_output_files/"])
    # TODO: create a quick function that can simplify this repetitive need to check if the path already exist
    # Next line checks to see if a plink output file already exist and delets it if it does
    if os.path.exists("".join(
        [args.output, "plink_output_files/", "plink_log.log"])):

        os.remove("".join(
            [args.output, "plink_output_files/", "plink_log.log"]))
    # TODO: Refactor this next line from line 93
    if args.var_file:
        analysis_type_checker: object = plink_initial_format_scripts.Analysis_Checker(
            ANALYSIS_TYPE,
            args.recode_options,
            args.binary_file,
            args.output,
            plink_file_path,
            var_file=args.var_file,
            maf_filter=MAF_FILTER)

        analysis_type_checker.check_analysis(
            readme_txt=utility_scripts.plink_readme_body_text,
            output=plink_file_path)

        analysis_type_checker.check_missing_var_count()

    else:
        analysis_type_checker: object = plink_initial_format_scripts.Analysis_Checker(
            ANALYSIS_TYPE,
            args.recode_options,
            args.binary_file,
            args.output,
            plink_file_path,
            maf_filter=MAF_FILTER)

        analysis_type_checker.check_analysis(
            range=args.range,
            readme_txt=utility_scripts.plink_readme_body_text,
            output=plink_file_path)

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
    carrier_analysis_scripts.single_variant_analysis(
        recode_filepath=plink_file_path,
        output=args.output,
        pop_info=args.pop_info,
        pop_code=args.pop_code,
        readme_output="".join([args.output, "carrier_analysis_output/"]),
        readme_text=utility_scripts.carrier_analysis_body_text,
    )
    logger.info(
        f"Writing the results of the carrier analysis called: {''.join([args.output, 'carrier_analysis_output/'])}"
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
        logger.warning(
            f"The variants {', '.join(variants_above_threshold)} were above the arbitrary threshold of {THRESHOLD}"
        )

        logger.warning(
            "Ending program so that the user can remove the above variants")
        sys.exit(1)

    elif program_end_code == 1 and variants_above_threshold:
        logger.warning(
            f"The variants {', '.join(variants_above_threshold)} were above the arbitrary threshold of {THRESHOLD}"
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

        convert_ibd_func_param: dict = {
            "ibd_file_path": ibd_file,
            "carrier_file": "".join([args.output, "carrier_analysis_output/"]),
            "ibd_program": program,
            "output_path": IBD_search_output_files,
            "map_files": "".join([args.output, "plink_output_files/"]),
            "ibd_file_suffix": file_suffix,
            "min_CM_threshold": MIN_CM,
            "threads": THREADS
        }

        # pre_shared_segments_analysis_scripts.convert_ibd(
        #     ibd_file, "".join([args.output, "carrier_analysis_output/"]),
        #     program, IBD_search_output_files,
        #     "".join([args.output,
        #              "plink_output_files/"]), file_suffix, MIN_CM, THREADS)
        pre_shared_segments_analysis_scripts.convert_ibd(
            convert_ibd_func_param,
            readme_output=IBD_search_output_files,
            readme_text=utility_scripts.formatted_ibd_dir_body_text_1)

    print("combining segment output...")
    ibd_dir_dict: dict = {"ilash": ILASH_PATH, "hapibd": HAPIBD_PATH}
    pre_shared_segments_analysis_scripts.combine_output(
        "".join([IBD_search_output_files, "reformatted_ibd_output/"]),
        args.ibd_programs, IBD_search_output_files,
        "".join([args.output, "carrier_analysis_output/reformatted/"]),
        ibd_dir_dict, "".join([args.output, "plink_output_files/"]))

    pre_shared_segments_analysis_scripts.reformat_files(
        "".join([args.output, "carrier_analysis_output/"]),
        "".join([args.output, "plink_output_files/"]),
        IBD_search_output_files,
        IBD_search_output_files,
        "".join([IBD_search_output_files, "no_carriers_in_file.txt"]),
    )
    logger.info(
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

    logger.info(
        f"Writing the results of the network analysis to: {''.join([args.output, 'networks/'])}"
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

    logger.info(
        f"Writing the output of the haplotype analysis to: {''.join([args.output, 'haplotype_analysis/'])}"
    )
    print(
        "generating a file contain the genotypes for all IIDs in the provided file"
    )

    full_analysis.get_all_genotypes(
        "".join([args.output, "plink_output_files/"]),
        "".join([IBD_search_output_files, "confirmed_carriers.txt"]),
        args.pop_info, "".join([args.output,
                                "haplotype_analysis/"]), args.pop_code)

    logger.info(
        f"Writing the result of getting all the genotypes for all IIDs in the provided file to: {''.join([args.output, 'haplotype_analysis/'])}"
    )

    # Add the function that gets the haplotype string to this
    logger.info(
        f"Writing the most probable haplotypes for each network to the file 'network_haplotypes.txt' at {''.join([args.output, 'haplotype_analysis/'])}"
    )
    print("Finding the most probable haplotypes")

    haplotype_segments_analysis.gather_haplotypes(
        haplotype_info_path, "".join([args.output, "haplotype_analysis/"]),
        args.binary_file, args.pop_info, args.pop_code)

    logger.info('Analysis finished...')

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
