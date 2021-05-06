#!/usr/bin/env python

import argparse
import os.path
from os import path
import sys
from datetime import datetime
from collections import namedtuple

import parser_arguments
import carrier_analysis_scripts
import create_network_scripts
import run_plink
import pre_shared_segments_analysis_scripts
import pre_shared_segments_analysis_scripts.shared_segment_detection
from pre_shared_segments_analysis_scripts.shared_segment_detection.combine_output.combine_ibd_pairs import combine_output
import utility_scripts
import collect_phenotype_info


@utility_scripts.func_readme_generator
def run(args: list, **kwargs: dict):

    # Next few lines give settings for the logger

    # Creating a logfile with the current day as part of the log file name
    current_day = datetime.now().strftime("%m_%d_%Y")

    file_name: str = "".join([args.output, current_day, "_run.log"])

    logger: object = utility_scripts.create_logger(__name__, file_name)

    utility_scripts.record_user_arguments(logger, args)

    # creating an object that will ask for user input for the Analysis type, the minimum  centimorgan threshold, the threads to be used and the
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

    # need to have it create branch where it will either do a 
    # # phenotype analysis or begin running plink
    if ANALYSIS_TYPE == "phenotype":
        logger.info("Beginning phenotype based analysis")
        # loading all the necessary files into dataframes

        pheno_df, pheno_carriers_df = collect_phenotype_info.load_pheno_file(args.pheno_file, args.phenotype_carriers)

    else:
        # output path that all the plink files will output to
        plink_file_path = "".join([args.output, "plink_output_files/"])
        # TODO: create a quick function that can simplify this repetitive need to check if the path already exist
        # Next line checks to see if a plink output file already exist and delets it if it does
        # TODO: Need tp refactpr the inputs to the plink _output_files part so that it can be decoratored to see if the file exist
        if os.path.exists("".join(
                [args.output, "plink_output_files/", "plink_log.log"])):

            os.remove("".join(
                [args.output, "plink_output_files/", "plink_log.log"]))
        # Refactor
        # TODO: Refactor this next line from line 93
        if args.var_file:
            analysis_type_checker: object = run_plink.Analysis_Checker(
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
            analysis_type_checker: object = run_plink.Analysis_Checker(
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
        arguments_dict = {
            "recode_filepath": plink_file_path,
            "output": args.output,
            "pop_info": args.pop_info,
            "pop_code": args.pop_code,
            "readme_output": "".join([args.output, "carrier_analysis_output/"]),
            "readme_text": utility_scripts.carrier_analysis_body_text,
        }

        carrier_analysis_scripts.single_variant_analysis(
            parameter_dict=arguments_dict
        )
        logger.info(
            f"Writing the results of the carrier analysis called: {''.join([args.output, 'carrier_analysis_output/'])}"
        )

        # The above function outputs files to a subdirectory called
        # "carrier_analysis_output"

        print(
            "determining the minor allele frequency within the provided binary file..."
        )
        carrier_analysis_scripts.determine_maf(
            "".join([args.output, "carrier_analysis_output/"]),
            "".join([args.output, "plink_output_files/"]), args.pop_info,
            args.pop_code, args.output)

        THRESHOLD: float = 0.10

        print(
            f"checking the variant minor allele frequencies against an arbitrary threshold of {THRESHOLD}"
        )
        # This next function will check to see if variants are above a specified threshold
        # If they are then the user has an option to end the program and remove the 
        # variants or just continue on with the program. The function will return a tuple 
        # where the first value is a list containing variants that are above a specified 
        # threshold and the second value is either 0 or 1 where 0 quits the program and 1 
        # continues

        variants_above_threshold_tuple: tuple = carrier_analysis_scripts.check_mafs(
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

    print("identifying pairs within region that shared IBD segments...")

    # If the user selects to run the program on phenotype then the analysis should start at this step effectively

    # # check to make sure the formatted_ibd_output exists
    IBD_search_output_files: str = utility_scripts.check_dir(args.output, "formatted_ibd_output/")


    for program in args.ibd_programs:

        suffix_dict: dict = {
            "ilash": "*.match.gz",
            "hapibd": "*.ibd.gz",
        }

        file_suffix: str = suffix_dict[program]

        # getting the correct ibd_file_path
        ibd_file: str = [
            file for file in IBD_PATHS_LIST if program in file.lower()
        ][0]

        if ANALYSIS_TYPE == "phenotype":

            pre_shared_segments_analysis_scripts.shared_segment_detection.gather_shared_segments(ibd_file, pheno_df, pheno_carriers_df, IBD_search_output_files, program, MIN_CM, file_suffix, THREADS)

        else:
            convert_ibd_func_param: dict = {
                "ibd_file_path": ibd_file,
                "carrier_file": "".join([args.output, "carrier_analysis_output/"]),
                "ibd_program": program,
                "output": IBD_search_output_files,
                "map_files": "".join([args.output, "plink_output_files/"]),
                "ibd_file_suffix": file_suffix,
                "min_CM_threshold": MIN_CM,
                "threads": THREADS,
                "readme_output":IBD_search_output_files,
                "readme_text":utility_scripts.formatted_ibd_dir_body_text_1
            }

        # Forming the file dictionary which is a file that contains the appropriate files for each chromosome
            file_dict: dict = pre_shared_segments_analysis_scripts.shared_segment_detection.collect_files(parameter_dict=convert_ibd_func_param)

            # iterating over this dictionary so that we can get the 
            # variants for each chromosome
            pre_shared_segments_analysis_scripts.shared_segment_detection.iterate_file_dict(file_dict, IBD_search_output_files, THREADS, program, MIN_CM)

    print("Identifying networks of pairs...")

    ibd_dir_dict: dict = {"ilash": ILASH_PATH, "hapibd": HAPIBD_PATH}

    if ANALYSIS_TYPE == "phenotype":
        gathered_file_dict: dict = pre_shared_segments_analysis_scripts.shared_segment_detection.gather_files(
        ibd_dir_dict,
        os.path.join(IBD_search_output_files, "collected_pairs/")
        )
    else:
    # getting a dictionary of all the files
        gathered_file_dict: dict = pre_shared_segments_analysis_scripts.shared_segment_detection.gather_files(
            ibd_dir_dict,
            os.path.join(IBD_search_output_files, "collected_pairs/"),
            map_file_dir=os.path.join(args.output, "plink_output_files/")
            )
    # getting a dictionary of all the files with the ibd files
    ibd_file_dict: dict = pre_shared_segments_analysis_scripts.shared_segment_detection.build_file_dict(gathered_file_dict["ibd_pair_file_list"], args.ibd_programs, ANALYSIS_TYPE)
    
    if ANALYSIS_TYPE == "phenotype":
        analysis_files: dict = {
            "pheno_gmap_df": pheno_df, 
            "pheno_carrier_df": pheno_carriers_df
        }
        combine_output(
            gathered_file_dict, ibd_file_dict, IBD_search_output_files, ANALYSIS_TYPE, THREADS,
            analysis_files)
        
        reformatter =  pre_shared_segments_analysis_scripts.shared_segment_detection.Pheno_Reformatter(
            IBD_search_output_files, 
            pheno_df,
            pheno_carriers_df,
            os.path.join(IBD_search_output_files, "pairs/")
        )

        reformatter.reformat()
    else:
        analysis_files: dict = {
            "carrier_dir": os.path.join(args.output, "carrier_analysis_output/"), 
        }

        combine_output(
            gathered_file_dict, ibd_file_dict, IBD_search_output_files, ANALYSIS_TYPE,THREADS, analysis_files)
        
        reformatter = pre_shared_segments_analysis_scripts.shared_segment_detection.Gene_Reformatter(
            os.path.join(args.output, "carrier_analysis_output/"),
            os.path.join(IBD_search_output_files, "pairs/"),
            os.path.join(args.output, "plink_output_files/"),
            os.path.join(IBD_search_output_files, "nopairs-identified.txt"),
            IBD_search_output_files
        )

        reformatter.reformat()
        
   
    logger.info(
        f"Writing the results from the IBD conversion files to: {IBD_search_output_files}\n"
    )

    print(
        "Identifying networks of individuals who share a segment"
    )

  

    # checking to make sure that the directory that the network files gets written to is real
    network_dir: str = utility_scripts.check_dir(args.output, "networks")


    create_network_scripts.create_networks(
        os.path.join(IBD_search_output_files, "pairs/"), network_dir, ANALYSIS_TYPE, os.path.join(IBD_search_output_files, "confirmed_carriers.txt"))

    logger.info(
        f"Writing the results of the network analysis to: {''.join([args.output, 'networks/'])}"
    )

   

    # if not path.exists("".join([args.output, "haplotype_analysis/"])):

    #     os.mkdir("".join([args.output, "haplotype_analysis/"]))

    # )
    # print(
    #     "generating a file contain the genotypes for all IIDs in the provided file"
    # )

    # full_analysis.get_all_genotypes(
    #     "".join([args.output, "plink_output_files/"]),
    #     "".join([IBD_search_output_files, "confirmed_carriers.txt"]),
    #     args.pop_info, "".join([args.output,
    #                             "haplotype_analysis/"]), args.pop_code)

    # logger.info(
    #     f"Writing the result of getting all the genotypes for all IIDs in the provided file to: {''.join([args.output, 'haplotype_analysis/'])}"
    # )

    

    logger.info('Analysis finished...')

    # TODO": Fix the timing issue so that it gives the correct time
    finishing_time = datetime.utcnow()
    print(
        f"The program successfully finished at {finishing_time.strftime('%H:%M:%S')}"
    )


def main():
    # creating the parser object
    parser = parser_arguments.create_args(run)
    
    # parser the user input arguments
    args = parser.parse_args()

    # getting the parameters for the run function
    readme_text: list = [utility_scripts.main_parameter_text,
                         utility_scripts.main_directory_header, utility_scripts.main_directory_text]

    args.func(args, parameter_dict={
              "readme_text": readme_text, "readme_output": args.output})


if __name__ == "__main__":
    main()
