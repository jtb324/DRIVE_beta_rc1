import typer
import os
import sys
from typing import List, Optional, Dict, MutableMapping
import toml
import models


import utility_scripts
import carrier_identification.run_plink as plink
from carrier_identification import carrier_analysis_scripts
# making the typer for this first program
carrier_id_app = typer.Typer(add_completion=False)

def main(recode_filepath: str, output: str, pop_info_file: Optional[str], pop_code: Optional[str], logger) -> Dict:
    """Function that will identify the carriers of specific variants
    Parameters
    __________
    recode_filepath : str
        string that list the path to the plink output files

    output : str
        string that list the output filepath that things will be written to
    
    pop_info_file : str
        string that list the filepath to the file that has information about the population demographics for the genotyped individuals
    
    pop_code : str
        string that list the population code that the 
        user wishes to filter for
    
    logger : logging.logger
        logger file that is information is being written to in order to 
        maintain records

    Returns
    _______
    Dict[str, Dict[str, List[str]]]
        returns a dictionary of dictionaries that contains the carriers of each variant for each chromosome
    """
    carrier_dir: str = utility_scripts.check_dir(output, "carrier_analysis_output/")
    # creating a readme that will provide info about the single variant analysis
    _ = utility_scripts.Readme_Info(carrier_dir, utility_scripts.carrier_analysis_body_text, "carrier_analysis_output_README.md")

    # running the single_variant_analysis function to 
    # determine individuals who carrier a variant of interest
    carrier_dict: Dict = carrier_analysis_scripts.single_variant_analysis(
        parameter_dict={
            "recode_filepath": recode_filepath,
            "output": output,
            "pop_info": pop_info_file,
            "pop_code": pop_code,    
        }
    )

    logger.info(
        f"Writing the results of which individuals are called as carrying a variant of interest to the directory at: {''.join([output, 'carrier_analysis_output/'])}"
    )

    print("\ndetermining the minor allele frequency within the provided binary file...")
    maf_dict: Dict = carrier_analysis_scripts.get_allele_frq(
        carrier_dict,
        "".join([output, "plink_output_files/"]),
        pop_info_file,
        pop_code,
        output,
    )

    THRESHOLD: float = 0.10

    print(
        f"checking the variant minor allele frequencies again against an arbitrary threshold of {THRESHOLD}"
    )

    # This next function will check to see if variants are above a specified threshold
    # If they are then the user has an option to end the program and remove the
    # variants or just continue on with the program. The function will return a tuple
    # where the first value is a list containing variants that are above a specified
    # threshold and the second value is either 0 or 1 where 0 quits the program and 1
    # continues

    variants_above_threshold_tuple: tuple = carrier_analysis_scripts.check_mafs(
        maf_dict,
        THRESHOLD,
    )

    variants_above_threshold: list = variants_above_threshold_tuple[0]

    program_end_code: int = variants_above_threshold_tuple[1]

    if program_end_code == 0:
        logger.warning(
            f"The variants {', '.join(variants_above_threshold)} were above the arbitrary threshold of {THRESHOLD}"
        )

        logger.warning("Ending program so that the user can remove the above variants")
        sys.exit(1)

    elif program_end_code == 1 and variants_above_threshold:
        logger.warning(
            f"The variants {', '.join(variants_above_threshold)} were above the arbitrary threshold of {THRESHOLD}"
        )
    return carrier_dict

def plink_runner(variant_file: str, recode_flags: List[str], binary_file: str, output: str, plink_file_path: str, MAF_THRESHOLD: float, range_list: List[int]) -> str:
    """function to run plink based on the inputs that the user provides
    Parameters
    __________
    variant_file : str
        string that list the filepath to either an xlsx or csv file that 
        has the variants of interest
    
    recode_flags : List[str]
        a list of strings that has the recode options passed to plink
    
    binary_file : str
        filepath to the binary file that has the genotyping from the 
        MEGA array
    
    output : str
        filepath to output the information to
    
    plink_file_path : str
        filepath to all of the ped and map files
    
    MAF_THRESHOLD : float
        float that list the threshold to filter variants above out
    
    range_list : List[int]
        list of integers where the first value is where to start 
        searching the genome and the second value is where to end
    
    Returns
    _______
    str
        string that list the direectory that the ped and map files were written to"""
    # running plink based off whether the user provides a file with just variants or wants to do all the variants within a specified range
    plink_output_path: str = os.path.join(output, "plink_output_files/")

    if variant_file:
        analysis_type_checker: object = plink.Analysis_Checker(
            "gene",
            recode_flags,
            binary_file,
            output,
            plink_output_path,
            var_file=variant_file,
            maf_filter=MAF_THRESHOLD,
        )

        analysis_type_checker.check_analysis(
            readme_txt=utility_scripts.plink_readme_body_text, output=plink_file_path
        )

        analysis_type_checker.check_missing_var_count()

    else:
        analysis_type_checker: object = plink.Analysis_Checker(
            "gene",
            recode_flags,
            binary_file,
            output,
            plink_file_path,
            maf_filter=MAF_THRESHOLD,
        )

        analysis_type_checker.check_analysis(
            range=range_list,
            readme_txt=utility_scripts.plink_readme_body_text,
            output=plink_file_path,
        )
    
    # having the function return the directory that the plink files were written to
    return plink_output_path
    
def import_toml(toml_filepath: str) -> models.Carrier_Parameters:
    """Function that will import the toml configurations and then """
    config: MutableMapping = toml.load(toml_filepath)

    # first need to create the InputParameters class object
    return models.InputParams(**config)
    # carrier_params: models.Carrier_Parameters = models.Carrier_Parameters(output, pop_info_file, pop_code, run_plink)


@carrier_id_app.command()
def determine_carriers(
    toml: str = typer.Option(
        "user_config.toml", 
        help="The filepath to the toml file containing the user configurations"
    ),
    run_plink: bool = typer.Option(
        True, help="Optional argument that will "
    )):
    """Function to return the IIDs of individuals carry at least one variant of interest"""
    
    # REFACTOR ###############################################
    
    # Setting a constant for the recode_flags used to run plink
    # removing any files from a previous run

    # using the users configuration file to gather all the inputs
    inputs: models.InputParams = import_toml(toml)

    analyzer: models.Carrier_Analyzer = models.Carrier_Analyzer(inputs.output_path, inputs.pop_info_file, inputs.pop_code, run_plink)

    analyzer.check_dir()
    # NOTE: 1/21/22 NEed to have some logic so that this is only run if run plink is specified
    if inputs.var_file == "": 
        models.Range_Runner(inputs.output_path, inputs.bfile, )
    else:
        models.Gene_Runner

    #NOTE: # Need to expand on these two above classes#NOTE
    

    # next two lines will create the logger and record initial parameters
    # logger: object = utility_scripts.create_logger(output, __name__)

    user_arg_dict: Dict = {
            "analysis_type": "Determining carriers",
            "genotype_input_file": binary_file,
            "output_directory": output,
            "variants_of_interest_file": variant_file,
            "demographic info file": pop_info_file,
            "population code used": pop_code,
            "gene range start position": RANGE_STR,
            "gene range end position": RANGE_END,
        }

    # Need to create a block that will run plink
    if run_plink:
        
        # asking the user to input a maf_threshold
        MAF_THRESHOLD: float = float(typer.prompt("Enter a minor allele frequency threshold to filter variants above the threshold"))

        print("running plink\n")

        # creating a constant for the recode flags passed to PLINK
        RECODE_FLAGS: List[str] = ["recode", "recodeA"]
        
        # recording the above two parameters
        user_arg_dict["Minor Allele Frequency Threshold"] = MAF_THRESHOLD
        user_arg_dict["recode flags used"] = RECODE_FLAGS
        # if plink is run then the ped_directory value 
        # is changed to the output directory provided by plink
        
        ped_directory = plink_runner(variant_file, RECODE_FLAGS, binary_file, output, ped_directory,MAF_THRESHOLD, [RANGE_STR, RANGE_END])
    

    print("\ngenerating list of individuals genotyped as carrying a variant of interest...\n")

    #running the main function that will determine if the 
    _ = main(ped_directory, output, pop_info_file, pop_code, logger)

    utility_scripts.record_user_arguments(
        logger,
        user_arg_dict,
    )




if __name__ == '__main__':
    carrier_id_app()