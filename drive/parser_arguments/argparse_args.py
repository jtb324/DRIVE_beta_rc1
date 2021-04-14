import argparse

def create_args(run_func: object) -> argparse.ArgumentParser:
    """Function to create the parser object
    Parameters
    __________
    run_func : Function
        function that will use the arguments from the parser object
    
    Returns
    _______
    argparse.ArgumentParser
        returns a parser object that will be use in the main function 
        of the DRIVE.py file
    """
    # creating the parser
    parser = argparse.ArgumentParser(
        description="This identifies individuals who have a specific variant in a raw file from PLINK"
    )

    # creating the arguments
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
        help="This flag will pass the analysis type if the user wants to first run plink for the variants. The flag expects the argument to either be 'gene' or 'maf', or phenotyoe. If the maf option is chosen than the user also needs to specify a start and end bp range and a chromosome as well as a minor allele frequency threshold",
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

    parser.add_argument(
        "--range",
        help="This argument will list the  start and end bp of the range you wish to look at. The argument should be formated like '--range START END'.",
        dest="range",
        nargs="+",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--pheno_file",
        help="This argument will list the file path to a file that contains information about the chromsome of interest and the start and end points",
        dest="pheno_file",
        type=str,
        required=False
    )
    parser.add_argument(
        "--pheno_carriers",
        help="This argument will list the file path to a file that contains information about grids identified as carriers based on phenotype",
        dest="phenotype_carriers",
        type=str,
        required=False
    )

    # setting the default run function
    parser.set_defaults(func=run_func)

    return parser