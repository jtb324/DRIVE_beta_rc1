# This file contains text that will be used in making the readmes of each file

# This next section contains text that will be put in the main readme file within the parent directory
main_parameter_text: str = "* *Information about parameters used during the run can be found in the file mega_run.log within the main directory*\n"

main_directory_header: str = "## Direectory Structure:\n"

main_directory_txt: str = "* **carrier_analysis_output**: This directory contains output that list iids that are identified as carrying a variant of interest on the MEGA_Ex array. These files usually end in single_variant_carrier.csv. This directory also has a file 'allele_frequencies.txt' that describes the allele frequencies for each variant based on Piper's dataset.\n\n* **formated_ibd_output**: This directory contains the output from the steps of the program that convert the ibd_files to a more human readable format. These files are the .small.txt.gz files. There are also output files from the steps that determine which pairs share segments. These files are the allpair.txt files. This directory also contains a file titled confirmed_carriers.txt that list which of the carriers identified on the MEGA array are actual confirmed carriers based on the shared segment analysis.\n\n* **networks**: This directory contains information about the percent of confirmed carriers for each variant (this is found in the carriers_in_networks.csv file) and about which network the confirmed carriers belong to (This is found in the network_groups.csv file within the network_imgs subdirectory).\n\n* **plink_output_files**: This directory contains the output from running the recode, recodeA, and freq option for plink.\n\n* **shell_scripts**:This directory contains the script that was used to run the program"

# The next section contains text that will be inserted into the different
plink_readme_body_text: str = "* This directory contains the output from plink while running the initial steps of the program. There should be a ped, map, raw, freq file, and a log file output by plink.\nThis files will be split up per chromosome. The main three types of files used during the analysis are the ped, map, and the raw files.\n"

# This next section contains text that will be placed within the readme of the carrier_analysiss_output directory
carrier_analysis_body_txt: str = "* This directory contains the output that describes which IIDs were detected as carrier specific variants of interest"

# This section contains text for the README that will be placed in the
# formated_ibd_output directory
formated_ibd_dir_body_text_1: str = "* This directory contains the output from the mid steps where grids that share segments are identified from the ILASH and hapibd files\n"

formated_ibd_dir_body_text_2: str = "* *More infomation about the location of the iLASH and hapid files cn be found in the mega_run.log file*"
