#!/usr/bin/env python


import os.path
from typing import List, Optional, Dict
from os import path
import typer
import sys
from datetime import datetime


import carrier_analysis_scripts
import create_network_scripts
import run_plink
import pre_shared_segments_analysis_scripts
import pre_shared_segments_analysis_scripts.shared_segment_detection
from pre_shared_segments_analysis_scripts.shared_segment_detection.combine_output.combine_ibd_pairs import (
    combine_output,
)
import utility_scripts
import collect_phenotype_info

from carrier_identification import carrier_id_main
from segment_analysis import segment_analysis_main

# creating a cli object
cli = typer.Typer(
    add_completion=False,
    help="Tool identifies networks of individuals who share an IBD segment around a genomic region of interest",
)
# create a subcommand for the carrier identification step
cli.add_typer(carrier_id_main.carrier_id_app, name="carrier_identification")

cli.add_typer(segment_analysis_main.segment_analysis_app, name="segment_analysis")




@cli.command()
def gene_analysis(
    binary_file: str = typer.Argument(
        ..., help="This is the binary file that has the genotyping information"
    ),
    output: str = typer.Option(
        "./", help="file directory to output all output files into"
    ),
    pop_info_file: Optional[str] = typer.Argument(
        None,
        help="filepath to the file that has all the ancestry information of the individuals in the binary file",
    ),
    pop_code: Optional[str] = typer.Argument(
        None,
        help="Optional population code for the data to be subset into. This is based on 1000 genomes",
    ),
    range_str: Optional[str] = typer.Argument(
        None,
        help="Genomic start position to extract a range of snps from. This is used for a gene approach",
    ),
    variant_file: Optional[str] = typer.Argument(
        None,
        help="This is the file that list variants of interest and there genomic information",
    ),
    range_end: Optional[str] = typer.Argument(
        None,
        help="Genomic end position to extract a range of snps from. This is used for a gene approach",
    ),
):



    

    print("identifying pairs within region that shared IBD segments...")

    IBD_search_output_files: str = os.path.join(output, "formatted_ibd_output/")

    #This section is going through each program
    for program in IBD_programs.split(" "):

        suffix_dict: dict = {
            "ilash": "*.match.gz",
            "hapibd": "*.ibd.gz",
        }

        file_suffix: str = suffix_dict[program]

        # getting the correct ibd_file_path
        ibd_file: str = [file for file in IBD_PATHS_LIST if program in file.lower()][0]

        # TODO: consider adding a readme class since the readme_output and readme_text show up in several spots
        convert_ibd_func_param: dict = {
            "ibd_file_path": ibd_file,
            "map_files": os.path.join([output, "plink_output_files/"]),
            "ibd_file_suffix": file_suffix,
        }

        ibd_readme_info: utility_scripts.Readme_Info = utility_scripts.Readme_Info(IBD_search_output_files, utility_scripts.formatted_ibd_dir_body_text_1)

        # Forming the file dictionary which is a file that 
        # contains the appropriate files for each chromosome
        file_dict: dict = (
            pre_shared_segments_analysis_scripts.shared_segment_detection.collect_files(
                parameter_dict=convert_ibd_func_param
            )
        )

        # iterating over this dictionary so that we can get the
        # variants for each chromosome
        pre_shared_segments_analysis_scripts.shared_segment_detection.iterate_file_dict(
            file_dict, IBD_search_output_files, THREADS, program, MIN_CM
        )

    print("Identifying networks of pairs...")

    ibd_dir_dict: dict = {"ilash": ILASH_PATH, "hapibd": HAPIBD_PATH}

    # getting a dictionary of all the files
    gathered_file_dict: dict = (
        pre_shared_segments_analysis_scripts.shared_segment_detection.gather_files(
            ibd_dir_dict,
            os.path.join(IBD_search_output_files, "collected_pairs/"),
            map_file_dir=os.path.join(output, "plink_output_files/"),
        )
    )

    ibd_file_dict: dict = (
        pre_shared_segments_analysis_scripts.shared_segment_detection.build_file_dict(
            gathered_file_dict["ibd_pair_file_list"], IBD_programs.split(" "), "gene"
        )
    )

    analysis_files: dict = {
        "carrier_dir": os.path.join(output, "carrier_analysis_output/"),
    }

    combine_output(
        gathered_file_dict,
        ibd_file_dict,
        IBD_search_output_files,
        "gene",
        THREADS,
        analysis_files,
    )

    reformatter = (
        pre_shared_segments_analysis_scripts.shared_segment_detection.Gene_Reformatter(
            os.path.join(output, "carrier_analysis_output/"),
            os.path.join(IBD_search_output_files, "pairs/"),
            os.path.join(output, "plink_output_files/"),
            os.path.join(IBD_search_output_files, "nopairs-identified.txt"),
            IBD_search_output_files,
        )
    )

    reformatter.reformat()

    logger.info(
        f"Writing the results from the IBD conversion files to: {IBD_search_output_files}\n"
    )

    print("Identifying networks of individuals who share a segment")

    # checking to make sure that the directory that the network files gets written to is real
    network_dir: str = utility_scripts.check_dir(args.output, "networks")

    create_network_scripts.create_networks(
        os.path.join(IBD_search_output_files, "pairs/"),
        network_dir,
        "gene",
        os.path.join(IBD_search_output_files, "confirmed_carriers.txt"),
    )

    logger.info(
        f"Writing the results of the network analysis to: {os.path.join(output, 'networks/')}"
    )

    logger.info("Analysis finished...")

    # TODO": Fix the timing issue so that it gives the correct time
    finishing_time = datetime.utcnow()
    print(f"The program successfully finished at {finishing_time.strftime('%H:%M:%S')}")


@cli.command()
def phenotype_analysis(
    binary_file: str = typer.Argument(
        ..., help="This is the binary file that has the genotyping information"
    ),
    output: str = typer.Option(
        "./", help="file directory to output all output files into"
    ),
    pop_info_file: Optional[str] = typer.Argument(
        None,
        help="filepath to the file that has all the ancestry information of the individuals in the binary file",
    ),
    pop_code: Optional[str] = typer.Argument(
        None,
        help="Optional population code for the data to be subset into. This is based on 1000 genomes",
    ),
    pheno_gmap: str = typer.Argument(
        None,
        help="filepath to the file that contains information about the chromosome of interest and the start and end position of the gene of interest",
    ),
    pheno_carriers: str = typer.Argument(
        ...,
        help="Filepath to a file that contains information about grids identified as carriers based on phenotype",
    ),
    MIN_CM: int = typer.Option(
        3,
        help="minimum centimorgan threshold that will be used to find the minimum size of shared IBD segments between pairs",
        prompt=True,
        callback=utility_scripts.check_provided_integer,
    ),
    THREADS: int = typer.Option(
        3,
        help="number of threads to parallelize the process to. This value has to be at least 3",
        prompt=True,
        callback=utility_scripts.check_provided_integer,
    ),
    IBD_programs: str = typer.Option(
        ...,
        help="Enter the ibd programs whose output was used in the analysis. Enter each program separated by a space",
        prompt=True,
    ),
):

    # remove files from a prior run
    utility_scripts.remove_dir(output)

    logger: object = utility_scripts.create_logger(output, __name__)

    logger.info("Beginning phenotype based analysis")

    utility_scripts.record_user_arguments(
        logger,
        {
            "genotype_input_file": binary_file,
            "output_directory": output,
            "demographic info file": pop_info_file,
            "population code used": pop_code,
            "minimum centimorgan threshold": MIN_CM,
            "phenotype genome information file": pheno_gmap,
            "affected individual files": pheno_carriers,
            "Thread count": THREADS,
            "IBD program output used": IBD_programs,
        },
    )

    # location the filepath to the ibd files
    logger.info(f"iLASH files directory: {ILASH_PATH}")

    logger.info(f"Hapibd files directory: {HAPIBD_PATH}")

    # loading all the necessary files into dataframes
    pheno_df, pheno_carriers_df = collect_phenotype_info.load_pheno_file(
        pheno_gmap, pheno_carriers
    )

    print("identifying pairs within region that shared IBD segments...")

    # If the user selects to run the program on phenotype then the analysis should start at this step effectively

    # # check to make sure the formatted_ibd_output exists
    IBD_search_output_files: str = os.path.join(output, "formatted_ibd_output/")

    for program in IBD_programs.split(" "):

        suffix_dict: dict = {
            "ilash": "*.match.gz",
            "hapibd": "*.ibd.gz",
        }

        file_suffix: str = suffix_dict[program]

        # getting the correct ibd_file_path
        ibd_file: str = [file for file in IBD_PATHS_LIST if program in file.lower()][0]

        pre_shared_segments_analysis_scripts.shared_segment_detection(
            ibd_file,
            pheno_df,
            pheno_carriers_df,
            IBD_search_output_files,
            program,
            MIN_CM,
            file_suffix,
            THREADS,
        )

    print("Identifying networks of pairs...")

    ibd_dir_dict: dict = {"ilash": ILASH_PATH, "hapibd": HAPIBD_PATH}

    gathered_file_dict: dict = (
        pre_shared_segments_analysis_scripts.shared_segment_detection.gather_files(
            ibd_dir_dict, os.path.join(IBD_search_output_files, "collected_pairs/")
        )
    )

    ibd_file_dict: dict = (
        pre_shared_segments_analysis_scripts.shared_segment_detection.build_file_dict(
            gathered_file_dict["ibd_pair_file_list"],
            IBD_programs.split(" "),
            "phenotype",
        )
    )

    analysis_files: dict = {
        "pheno_gmap_df": pheno_df,
        "pheno_carrier_df": pheno_carriers_df,
    }

    combine_output(
        gathered_file_dict,
        ibd_file_dict,
        IBD_search_output_files,
        "phenotype",
        THREADS,
        analysis_files,
    )

    reformatter = (
        pre_shared_segments_analysis_scripts.shared_segment_detection.Pheno_Reformatter(
            IBD_search_output_files,
            pheno_df,
            pheno_carriers_df,
            os.path.join(IBD_search_output_files, "pairs/"),
        )
    )

    reformatter.reformat()

    logger.info(
        f"Writing the results from the IBD conversion files to: {IBD_search_output_files}\n"
    )

    print("Identifying networks of individuals who share a segment")

    # checking to make sure that the directory that the network files gets written to is real
    network_dir: str = utility_scripts.check_dir(output, "networks")

    create_network_scripts.create_networks(
        os.path.join(IBD_search_output_files, "pairs/"),
        network_dir,
        "phenotype",
        os.path.join(IBD_search_output_files, "confirmed_carriers.txt"),
    )

    logger.info(
        f"Writing the results of the network analysis to: {os.path.join(output, 'networks/')}"
    )

    logger.info("Analysis finished...")

    # TODO": Fix the timing issue so that it gives the correct time
    finishing_time = datetime.utcnow()
    print(f"The program successfully finished at {finishing_time.strftime('%H:%M:%S')}")


if __name__ == "__main__":
    cli()
