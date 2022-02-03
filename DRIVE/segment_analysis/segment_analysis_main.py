import typer
import os
from typing import List, Dict, Optional
from datetime import datetime
from pathlib import Path
import shutil
import pandas as pd
from dataclasses import dataclass
import callbacks
import utility_scripts
import segment_analysis.pre_shared_segments_analysis_scripts.shared_segment_detection as shared_segment_detection
from segment_analysis import pre_shared_segments_analysis_scripts
from segment_analysis import collect_phenotype_info
from segment_analysis import create_network_scripts

segment_analysis_app = typer.Typer()

@dataclass
class Directories:
    output: str
    ibd_programs: str
    ibd_file_dir: Dict[str, str]
    plink_files: Optional[str] = None
    affected_file: Optional[str] = None
    pheno_gene_info: Optional[str] = None
    pheno_carriers: Optional[str] = None 
    
    

def determine_segments(inputs: Directories, MIN_CM: int, THREADS: int, plink_file_path: Optional[str], logger) -> str:
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

    # ibd_readme_info: utility_scripts.Readme_Info = utility_scripts.Readme_Info(IBD_search_output_files, utility_scripts.formatted_ibd_dir_body_text_1, "formatted_ibd_output_README.md")

    # creating a dictionary that can keep track of information about the pairs
    pair_info_dict: Dict[str, Dict] = {}

 
    #This section is going through each program
    for program in inputs.ibd_programs.split(" "):

        suffix_dict: Dict[str, str] = {
            "ilash": "*.match.gz",
            "hapibd": "*.ibd.gz",
        }

        file_suffix: str = suffix_dict[program.lower()]

        # getting the correct ibd file directory based on the program
        ibd_dir: str = inputs.ibd_file_dir[program.lower()]

        # This part may be better done by polymorphism rather than what i have
        # if the program provides the files necessary for the phenotype 
        # analysis then they will be loaded into dataframes
        break
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
                                                            THREADS, pair_info_dict
                                                            )


        else: 
            convert_ibd_func_param: Dict = {
            "ibd_file_path": ibd_file, 
            "map_files": plink_file_path,
            "ibd_file_suffix": file_suffix,
            }


            # Forming the file dictionary which is a file that 
            # contains the appropriate files for each chromosome
            file_dict: Dict = (
                shared_segment_detection.collect_files(
                    convert_ibd_func_param, os.path.join(carrier_file_dir, "single_variant_carriers.csv")
                )
            )

            # iterating over this dictionary so that we can get the
            # variants for each chromosome
            shared_segment_detection.iterate_file_dict(
                file_dict, IBD_search_output_files, THREADS, program, MIN_CM, pair_info_dict
            ) 

    
    if pheno_gmap and pheno_carriers:

        gathered_file_dict: Dict = (
        shared_segment_detection.gather_files(
            IBD_dir_dict, os.path.join(IBD_search_output_files, "collected_pairs/")
            )
        )

        ibd_file_dict: Dict = (
            shared_segment_detection.build_file_dict(
                gathered_file_dict["ibd_pair_file_list"],
                IBD_programs.split(" "),
                "phenotype",
            )
        )

        analysis_files: Dict = {
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
            pair_info_dict
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
        gathered_file_dict: Dict = (
            shared_segment_detection.gather_files(
                IBD_dir_dict,
                os.path.join(IBD_search_output_files, "collected_pairs/"),
                map_file_dir=os.path.join(output, "plink_output_files/"),
            )
        )

        ibd_file_dict: Dict = (
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
            pair_info_dict
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

        create_network_scripts.create_distributions(network_dir, "phenotype")

    else: 
        create_network_scripts.create_networks(
            os.path.join(ibd_file_dir, "pairs/"),
            network_dir,
            "gene",
            os.path.join(ibd_file_dir, "confirmed_carriers.txt"),
        )

        create_network_scripts.create_distributions(network_dir, "gene")

    logger.info(
        f"Writing the results of the network analysis to: {network_dir}"
    )

    logger.info("Analysis finished...")

    # TODO": Fix the timing issue so that it gives the correct time
    finishing_time = datetime.utcnow()
    print(f"The program successfully finished at {finishing_time.strftime('%H:%M:%S')}")

def load_env_vars(verbose: bool, loglevel: str, cM: int, cpus: int) -> None:
    """Function that will load the verbose flag and the loglevel flag into the environmental variables
    
    Parameters

    verbose : bool
        verbose status of the program. This will either be True or False
    
    loglevel : str
        string value that has what loglevel the user wishes to view
    
    cM : int
        minimum centimorgan threshold to use during the analysis
    
    cpus : int
        number of cpu cores to use during the analysis
    """
    os.environ['verbose'] = str(verbose)
    os.environ['loglevel'] = loglevel
    os.environ['cM'] = str(cM)
    os.environ['cpus'] = str(cpus)

def setup_output_directories(output_paths: List[str]) -> None:
    """Function that will see if the output directories exists. If they don't it will 
    create them. If it does then it will make sure they are empty.
    
    Parameters

    output_paths : List[str]
        directories that the files will be written to
    """
    for directory in output_paths:
        path: Path = Path(directory)

        if path.exists():
            shutil.rmtree(directory, ignore_errors=False)
            path.mkdir()
        else:
            path.mkdir()
    
@segment_analysis_app.command()
def main(
    output: str = typer.Option(
        ..., 
        "--output", 
        "-o",
        help="Filepath that the output files will be written to"
    ),
    plink_files: str = typer.Option(
        None,
        help="Filepath to the directory that has the ped and map files from PLINK"
    ),
    carrier_files: str = typer.Option(
        None,
        "--carriers",
        "-c",
        help="Filepath to the the csv that list all the individuals who carry a certain variant."
    ),
    pheno_gmap: str = typer.Option(
        None,
        "--pheno-gene-info",
        "-pg",
        help="filepath to the file that contains information about the chromosome of interest and the start and end position of the gene of interest",
    ),
    pheno_carriers: str = typer.Option(
        None,
        "--pheno_carriers",
        "-pc",
        help="Filepath to a file that contains information about grids identified as carriers based on phenotype",
    ),
    MIN_CM: int = typer.Option(
        3,
        help="minimum centimorgan threshold that will be used to find the minimum size of shared IBD segments between pairs",
        prompt=True,
        callback=callbacks.check_provided_integer,
    ),
    THREADS: int = typer.Option(
        3,
        help="number of threads to parallelize the process to. This value has to be at least 3",
        prompt=True,
        callback=callbacks.check_provided_integer,
    ), 
    IBD_programs: str = typer.Option(
        ...,
        help="Enter the ibd programs whose output was used in the analysis. Enter each program separated by a space",
        prompt=True
    ),
    verbose: bool = typer.Option(
        False, 
        "--verbose", 
        "-v", 
        help="Optional Flag to run the program in verbose mode", 
        is_flag=True
    ),
    loglevel: str = typer.Option(
        "INFO",
        "--loglevel",
        "-l",
        help="This argument sets the logging level for the program"
    )):

    load_env_vars(verbose, loglevel, MIN_CM, THREADS)

    setup_output_directories(
        [os.path.join(output, "formatted_ibd_output"), os.path.join(output, "networks")]
        )

    # creating a logger object
    logger: object = utility_scripts.create_logger(output, __name__)

    # recording the users arguments that were passed to the function
    # utility_scripts.record_user_arguments(
    #     logger,
    #     {
    #         "analysis_type": "Shared Segment Analysis",
    #         "output_directory": output,
    #         "minimum centimorgan threshold": MIN_CM,
    #         "number of threads used": THREADS,
    #         "minor allele frequency threshold": MAF_THRESHOLD,
    #         "IBD programs used": IBD_programs,
    #     }
    # )
    
    # This list groups the ibd files together and they will be used later in the program
    ibd_dir_dict: Dict = {"ilash": os.environ.get("ilash_files"), "hapibd": os.environ.get("hapibd_files")}

    # logger.info(f"iLASH files directory: {ibd_dir_dict['ilash']}")

    # logger.info(f"Hapibd files directory: {ibd_dir_dict['ihapibd']}")

    print("identifying pairs within region that shared IBD segments...")

    input_dirs: Directories = Directories(output, IBD_programs, ibd_dir_dict plink_files, carrier_files, pheno_gmap, pheno_carriers)

    IBD_search_output_files: str = determine_segments(input_dirs, MIN_CM, THREADS, plink_files, logger)

    # IBD_search_output_files: str = utility_scripts.check_dir(output, "formatted_ibd_output/")

    determine_networks(output, IBD_search_output_files, pheno_carriers!=None and pheno_gmap!=None, logger)