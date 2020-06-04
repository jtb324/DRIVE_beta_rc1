# This is the main file that uses the VariantID.py and the MultipleVariantID.py to determine
# the IDs of individuals who contain a variant

from MEGA_CLI_project.VariantID import individualVariantID
from MEGA_CLI_project.MultipleVariantID import multiVariantID


def IdentifyID():
    with open("/data100t1/home/james/MEGA_CFTR_test_case/CFTR_snps_extracted/low_maf_CFTR_recodeA.raw") as geno_file:

        headerLine = next(geno_file)  # This skips the 1st row

        # This next two lines create lists for the total variants and the multivariants ids
        totalVariantList = []

        multiVariantList = []

        for row in geno_file:  # This iterates through each row in the file

            row = row.split()  # This will split the row by white space
            # Next two lines run the corresponding functions
            individualVariantID(row, totalVariantList)

            multiVariantID(row, multiVariantList)

        # The next two lines print the length of the two list made which corresponds to number of ids found
        print(len(totalVariantList))

        print(len(multiVariantList))
