import typer
import os
from typing import List, Dict

import utility_scripts
segment_analysis_app = typer.Typer()

@segment_analysis_app.command()
def main(
    output: str = typer.Option(
        ..., 
    help="Filepath that the output files will be written to"),
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

    utility_scripts.remove_dir(output, "networks")

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
    IBD_PATHS_LIST: List[str] = [ILASH_PATH, HAPIBD_PATH]

    logger.info(f"iLASH files directory: {ILASH_PATH}")

    logger.info(f"Hapibd files directory: {HAPIBD_PATH}")

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