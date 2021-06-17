import typer
import os
from typing import List, Dict, Optional
from datetime import datetime
import pandas as pd

import utility_scripts
import segment_analysis.pre_shared_segments_analysis_scripts.shared_segment_detection as shared_segment_detection
from segment_analysis import pre_shared_segments_analysis_scripts
from segment_analysis import collect_phenotype_info
from segment_analysis import create_network_scripts

segment_analysis_app = typer.Typer()

def determine_segments(output: str, carrier_file_dir: Optional[str], IBD_dir_dict: Dict[str, str], IBD_programs: str, pheno_gmap: Optional[str], pheno_carriers: Optional[str], MIN_CM: int, THREADS: int, plink_file_path: Optional[str], logger) -> str:
    """Function to determine the shared segments between pairs
    Parameters
    __________
    output : str
        string that list the directory to output all the files into
    
    IBD_dir_dict: Dict[str, str]
        Dictionary that has the name of ilash and hapibd 
        for keys and the filepath for values the 
        filepath to the directory for the hapibd files 
        and the ilash files

    IBD_programs : str
        string that list the IBD programs that IBD 
        information comes from for example hapibd and 
        ilash. These should be separaated by a space
    
    returns
    _______
    str
        returns a string which is the filepath for the 
        main directory of the IBD output
    """
    # making sure the directory exists that the ibd 
    # output will be put into and returning that file 
    # path

    IBD_search_output_files: str = utility_scripts.check_dir(output, "formatted_ibd_output/")
    
    #This section is going through each program
    for program in IBD_programs.split(" "):

        suffix_dict: dict = {
            "ilash": "*.match.gz",
            "hapibd": "*.ibd.gz",
        }

        file_suffix: str = suffix_dict[program]

        # getting the correct ibd_file_path
        ibd_file: str = [file for file in IBD_dir_dict.values() if program in file.lower()][0]

        # if the program provides the files necessary for the phenotype 
        # analysis then they will be loaded into dataframes
        if pheno_gmap and pheno_carriers:
            pheno_df, pheno_carriers_df = collect_phenotype_info.load_pheno_file(pheno_gmap, pheno_carriers)

            shared_segment_detection.gather_shared_segments(
            ibd_file,
            pheno_df,
            pheno_carriers_df,
            IBD_search_output_files,
            program,
            MIN_CM,
            file_suffix,
            THREADS
            )
        else: 
            convert_ibd_func_param: Dict = {
            "ibd_file_path": ibd_file, 
            "map_files": plink_file_path,
            "ibd_file_suffix": file_suffix,
            }

            ibd_readme_info: utility_scripts.Readme_Info = utility_scripts.Readme_Info(IBD_search_output_files, utility_scripts.formatted_ibd_dir_body_text_1)

            # Forming the file dictionary which is a file that 
            # contains the appropriate files for each chromosome
            file_dict: dict = (
                shared_segment_detection.collect_files(
                    convert_ibd_func_param, os.path.join(carrier_file_dir, "single_variant_carriers.csv")
                )
            )

            # iterating over this dictionary so that we can get the
            # variants for each chromosome
            shared_segment_detection.iterate_file_dict(
                file_dict, IBD_search_output_files, THREADS, program, MIN_CM
            ) 


    if pheno_gmap and pheno_carriers:
        gathered_file_dict: dict = (
        shared_segment_detection.gather_files(
            IBD_dir_dict, os.path.join(IBD_search_output_files, "collected_pairs/")
            )
        )

        ibd_file_dict: dict = (
            shared_segment_detection.build_file_dict(
                gathered_file_dict["ibd_pair_file_list"],
                IBD_programs.split(" "),
                "phenotype",
            )
        )

        analysis_files: dict = {
            "pheno_gmap_df": pheno_df,
            "pheno_carrier_df": pheno_carriers_df,
        }

        shared_segment_detection.combine_output(
            gathered_file_dict,
            ibd_file_dict,
            IBD_search_output_files,
            "phenotype",
            THREADS,
            analysis_files,
        )

        reformatter = (
            shared_segment_detection.Pheno_Reformatter(
                IBD_search_output_files,
                pheno_df,
                pheno_carriers_df,
                os.path.join(IBD_search_output_files, "pairs/"),
            )
        )

        reformatter.reformat()
    else:
        # getting a dictionary of all the files
        gathered_file_dict: dict = (
            shared_segment_detection.gather_files(
                IBD_dir_dict,
                os.path.join(IBD_search_output_files, "collected_pairs/"),
                map_file_dir=os.path.join(output, "plink_output_files/"),
            )
        )

        ibd_file_dict: dict = (
            shared_segment_detection.build_file_dict(
                gathered_file_dict["ibd_pair_file_list"], IBD_programs.split(" "), "gene"
            )
        )

        analysis_files: Dict[str, pd.DataFrame] = {
            "carrier_df": pd.read_csv(os.path.join(output, "carrier_analysis_output/single_variant_carriers.csv"), sep="\t"),
        }

        shared_segment_detection.combine_output(
            gathered_file_dict,
            ibd_file_dict,
            IBD_search_output_files,
            "gene",
            THREADS,
            analysis_files,
        )

        reformatter = (
            shared_segment_detection.Gene_Reformatter(
                pd.read_csv(os.path.join(output, "carrier_analysis_output/single_variant_carriers.csv"), sep="\t"),
                os.path.join(IBD_search_output_files, "pairs/"),
                os.path.join(output, "plink_output_files/"),
                os.path.join(output, "nopairs-identified.txt"),
                IBD_search_output_files,
            )
        )

        reformatter.reformat()

    logger.info(
        f"Writing the results from the IBD conversion files to: {output}\n"
    )

    return IBD_search_output_files

def determine_networks(output: str, ibd_file_dir: str, is_phenotype_analysis: bool, logger):
    print("Identifying networks of individuals who share a segment")

    # checking to make sure that the directory that the network files gets written to is real
    network_dir: str = utility_scripts.check_dir(output, "networks")
    
    if is_phenotype_analysis:
        create_network_scripts.create_networks(
        os.path.join(ibd_file_dir, "pairs/"),
        network_dir,
        "phenotype",
        os.path.join(ibd_file_dir, "confirmed_carriers.txt"))
    else: 
        create_network_scripts.create_networks(
            os.path.join(ibd_file_dir, "pairs/"),
            network_dir,
            "gene",
            os.path.join(ibd_file_dir, "confirmed_carriers.txt"),
        )

    logger.info(
        f"Writing the results of the network analysis to: {network_dir}"
    )

    logger.info("Analysis finished...")

    # TODO": Fix the timing issue so that it gives the correct time
    finishing_time = datetime.utcnow()
    print(f"The program successfully finished at {finishing_time.strftime('%H:%M:%S')}")

@segment_analysis_app.command()
def main(
    output: str = typer.Option(
        ..., 
        help="Filepath that the output files will be written to"
    ),
    ped_file_path: str = typer.Option(
        None,
        help="Filepath to the directory that has the ped and map files from PLINK"
    ),
    carrier_files: str = typer.Option(
        None,
        help="Filepath to the the csv that list all the individuals who carry a certain variant."
    ),
    pheno_gmap: str = typer.Option(
        None,
        help="filepath to the file that contains information about the chromosome of interest and the start and end position of the gene of interest",
    ),
    pheno_carriers: str = typer.Option(
        None,
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
    MAF_THRESHOLD: float = typer.Option(
        0.05,
        help="maf frequency threshold where only variants below the threshold will be kept",
        prompt=True,
        callback=utility_scripts.check_provided_maf,
    ),
    IBD_programs: str = typer.Option(
        ...,
        help="Enter the ibd programs whose output was used in the analysis. Enter each program separated by a space",
        prompt=True,
    )):

    #removing all files from a previous run
    # remove files from a prior run
    utility_scripts.remove_dir(os.path.join(output, "formatted_ibd_output"))

    utility_scripts.remove_dir(os.path.join(output, "networks"))

    # creating a logger object
    logger: object = utility_scripts.create_logger(output, __name__)

    # recording the users arguments that were passed to the function
    utility_scripts.record_user_arguments(
        logger,
        {
            "analysis_type": "Shared Segment Analysis",
            "output_directory": output,
            "minimum centimorgan threshold": MIN_CM,
            "number of threads used": THREADS,
            "minor allele frequency threshold": MAF_THRESHOLD,
            "IBD programs used": IBD_programs,
        }
    )
    
    # setting parameters for the ilash file paths and the hapibd file paths
    ILASH_PATH: str = "/data100t1/share/BioVU/shapeit4/Eur_70k/iLash/min100gmap/"

    HAPIBD_PATH: str = "/data100t1/share/BioVU/shapeit4/Eur_70k/hapibd/"
    # This list groups the ibd files together and they will be used later in the program
    ibd_dir_dict: dict = {"ilash": ILASH_PATH, "hapibd": HAPIBD_PATH}

    logger.info(f"iLASH files directory: {ILASH_PATH}")

    logger.info(f"Hapibd files directory: {HAPIBD_PATH}")

    print("identifying pairs within region that shared IBD segments...")


    IBD_search_output_files: str = determine_segments(output, carrier_files, ibd_dir_dict, IBD_programs, pheno_gmap, pheno_carriers, MIN_CM, THREADS, ped_file_path, logger)

    determine_networks(output, IBD_search_output_files, pheno_carriers!=None and pheno_gmap!=None, logger)