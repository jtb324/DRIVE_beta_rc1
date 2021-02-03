import sys

# This script will ask for the user to input certain initial parameters


# Function for analysis type
def ask_for_analysis_type(logfile) -> str:
    """
    This function will get the analysis type for the run.
    This can be either "maf", "gene", or ""
    """
    ANALYSIS_TYPE: str = input("Please input an analysis type: ").strip(" ")

    logfile.add_newline("INFO", f"Using the analysis type: {ANALYSIS_TYPE}\n")

    if ANALYSIS_TYPE not in ["gene", "maf", ""]:

        print(
            "Please choose one of the allowed analysis types: 'gene'/'maf'/'' "
        )

        logfile.add_newline("ERROR", "unrecognized analysis type choosen")

        sys.exit(1)

    return ANALYSIS_TYPE


def ask_for_min_cm(logfile) -> int:
    """
    This function will ask the user for the minimum centimorgan
    threshold to be used in analysis
    """
    MIN_CM: str = input(
        "Please input a value for the minimum centimorgan threshold: (Default is 3 cM) "
    )

    if not MIN_CM:
        MIN_CM = 3
    else:
        MIN_CM = int(MIN_CM)

    logfile.add_newline(
        "INFO", f"using a minimum centimorgan threshold of {MIN_CM}\n")

    return MIN_CM


def ask_for_thread_count(logfile) -> int:
    """
    This function will ask the user for the number of threads to be used 
    during the analysis. It returns an integer. The function defaults to 
    3 if the user provides no input
    """
    THREADS: str = input(
        "Please enter the number of threads you wish to use during this process. The default value is 3. (Bear in mind that this number will be used for all parallelized steps): "
    )

    if not THREADS:

        THREADS = 3

    else:

        THREADS = int(THREADS)

    logfile.add_newline("INFO", f"setting the thread count to be {THREADS}\n")

    return THREADS


def ask_for_maf_filter(logfile) -> str:
    """
    This function will ask for a minor allele frequency threshold from the 
    user. If the user provides no value then it will default to 0.05. The
    function returns a string
    """
    MAF_FILTER: str = input(
        "Please input a minor allele frequency threshold to be used. (The default is 0.05): "
    )

    if not MAF_FILTER:

        MAF_FILTER: str = '0.05'

    logfile.add_newline(
        "INFO", f"Using a minor allele frequency threshold of {MAF_FILTER}\n")

    return MAF_FILTER
