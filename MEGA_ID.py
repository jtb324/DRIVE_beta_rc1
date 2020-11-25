#! /usr/bin/env python
import logging
import argparse
import os.path
from os import path
import sys

import carrier_analysis_scripts
import create_network_scripts
import file_creator_scripts
import allele_frequency_analysis_scripts
import pre_shared_segments_analysis_scripts
import plink_initial_format_scripts
import haplotype_segments_analysis



def run(args):

    analysis_arguments: list = [
        "multi", "single", "draw_networks", "maf_determination", "shared_segments_preformat",
        "get_haplotype", "split_input_file"]

    if args.analysis.lower() not in analysis_arguments:
        print("The analysis method selected was not recognized. Please use one of the following analysis methods: \
            *multi\
            *single\
            *draw_networks\
            *plink_splitter\
            *maf_determination\
            *shared_segment_preformat\
            *get_haplotype\
            *split_input_file")

        # This line stops the program if the analysis method is not recognized
        sys.exit(1)

    if args.analysis == "multi":

        log_format = '%(asctime)s - %(levelname)s : %(message)s'

        logging.basicConfig(
            filename="".join([args.output, '/multi_variant_analysis.log']), level=logging.INFO,
            format=log_format)

        logging.info(
            'Generating a list of IIDs who carry multiple variants...')

        print("generating list of individuals carrying multiple variants....")

        carrier_analysis_scripts.multiVariantAnalysis(args.input, args.output,
                                                      args.compatible_format, 'multi_variant_list.csv')

        logging.info(
            'Finished creating a list of individuals carrying multiple variants')

    elif args.analysis == "single":

        log_format = '%(asctime)s - %(levelname)s : %(message)s'

        logging.basicConfig(filename="".join([args.output,
                                              '/single_variant_analysis.log']), level=logging.INFO,
                            format=log_format)

        logging.info("generating list of individuals at each probe id...")

        print("generating list of individuals at each probe id...")
        print(args.input)
        print(args.input[0])
        # The args.input should be a directory indicating where the raw files are located
        carrier_analysis_scripts.singleVariantAnalysis(args.input, args.output, args.compatible_format,
                                                       'single_variant_list.csv', args.pop_info, args.pop_code)

        logging.info(
            'Finished creating a file of individuals for each probe id')

    elif args.analysis == "draw_networks":

        print("generating pdf files of networks of individuals who share segments...")

        # This dictionaru keeps track of how many carriers are actually in the network. It needs to be a global variable so that it is just extended for each variant instead of recreated
        carrier_in_network_dict = dict()

        output = "".join([args.output, "network_imgs"])

        if not path.exists(output):
            os.mkdir(output)

        log_format = '%(asctime)s - %(levelname)s : %(message)s'

        logging.basicConfig(filename="".join([output,
                                              '/drawing_networks.log']), level=logging.INFO,
                            format=log_format)

        logging.info(
            "Drawing networks for all individuals identified as shared segments...")

        carrier_in_network_dict = create_network_scripts.create_networks(args.allpair_file, args.var_file, carrier_in_network_dict,
                                                  output)

        logger = logging.getLogger(args.output+'/carriers_in_network.log')

        # Writing the dictionary to a csv file
        csv_writer = file_creator_scripts.Csv_Writer_Object(
            carrier_in_network_dict, output, "carriers_in_networks.csv", logger)

        csv_writer.log_file_path()

        csv_writer.write_to_csv()

    elif args.analysis.lower() == "maf_determination":

        print("determine the minor allele frequency within the provided binary file")
        allele_frequency_analysis_scripts.determine_maf(args.carrier_file, args.raw_file,
                     args.pop_info, args.pop_code, args.output)

    elif args.analysis.lower() == "shared_segment_preformat":
        print("converting the ibd output to a human readable version...")

        pre_shared_segments_analysis_scripts.convert_ibd(args.ibd_files, args.carrier_file, args.ibd_programs,
                                                         args.output, args.map_file, args.ibd_file_suffix, args.min_CM,
                                                         args.threads, args.var_input_file)

        pre_shared_segments_analysis_scripts.combine_output(
            args.segment_files, args.ibd_programs, args.output, args.reformat_files)

        pre_shared_segments_analysis_scripts.reformat_files(args.carrier_file, args.plink_dir, args.allpair_file,
                                                            args.no_carriers_file)

    elif args.analysis.lower() == "split_input_file":
        print("splitting the input excel or csv file and then extracting the variants through PLINK")

        plink_initial_format_scripts.split_input_and_run_plink(args.input, args.output, args.recode_options,
                                                               args.binary_file, args.plink_dir)

    elif args.analysis.lower() == "get_haplotypes":
        print("getting information about the haplotypes for the confirmed carriers")

        haplotype_segments_analysis.get_segment_lengths(args.input, args.output, args.ilash_dir, args.hapibd_dir,
                                                        args.threads, args.plink_dir, args.reformat_files, args.allpair_file)
def main():
    parser = argparse.ArgumentParser(
        description="This identifies individuals who have a specific variant in a raw file from PLINK")

    parser.add_argument("--input", "-i", help="This is the pathway for the PLINK recoded input file. "
                                              "If you use the matchPED analysis argument then you should provide two input paths. The first is to the list of all variants. This should be a csv file and will have a list of the variant index and then a list of individuals who carry that variant. The second path is to the Pedigree file.At current development, this should be a .fam file.",
                        dest="input", nargs="+", type=str)

    parser.add_argument("--bfile", "-b", help="This argument will list the directory to the bim file which give them"
                                              "genotype information that PLINK uses",
                        dest="binary_file", type=str, required=False)

    parser.add_argument("--cfile", "-c", help="This argument supplies the path for the carrier files that commonly end in '.single_var_list.csv",
                        dest="carrier_file", type=str, required=False)

    parser.add_argument("--recode_options",
                        help="TThis argument list the recode option used to run plink",
                        dest="recode_options", nargs="+", type=str, required=False)

    parser.add_argument("--ilash_files",
                        help="This argument list the directory for all ilash files ending in .match.gz",
                        dest="ilash_dir", type=str, required=False)

    parser.add_argument("--hapibd_files",
                        help="This argument list the directory for all the hapibd files ending in .ibd.gz",
                        dest="hapibd_dir", type=str, required=False)

    parser.add_argument("--nfiles",
                        help="This argument list the path for the 'network_groups.csv file which list what network pairs are a part of'",
                        dest="network_files", type=str, required=False)

    parser.add_argument("--rfile", "-r", help="This argument supplies the path for the raw output files from PLINK. These files will be found in the directory ./variants_of_interest/ if the split_file was run",
                        dest="raw_file", type=str, required=False)

    parser.add_argument("--converted_ibd", help="This argument list the directory for all of the ibd files that were converted into a human readbale path. THese files should end in small.txt.gz",
                        dest="segment_files", type=str, required=False)

    parser.add_argument("--reformated_carriers", help="This argument list the directory for all of the reformated carrier files list. These files should be in a subdirected called reformat",
                        dest="reformat_files", type=str, required=False)

    parser.add_argument("--output", "-o", help="This is the directory that text files containing the ids of the individuals who have desired variants will be written to.",
                        dest="output", type=str, required=True)

    parser.add_argument("--analysis", help="This tag indicates that the multiVariantAnalysis function will be called "
                                           "to analyze how many individuals carry multiple variants. Two csv files are "
                                           "made which contain the indices of the variants and a list of the "
                                           "individuals that contain those variants. This accepts single, total, multi"
                                           ", matchPED, allele_counts, draw_networks",
                        dest="analysis", type=str, required=True)

    parser.add_argument("--thread", "-t", help="This argument list how many workers the program can use for multiprocessing",
                        dest="threads", type=str, required=False, default=1)

    parser.add_argument("--ifiles", "-ibd", help="This argument list the directory that contains the output from the ibd files",
                        dest="ibd_files", type=str, required=False)

    parser.add_argument("--plink_dir", "-pd", help="This argument list the directory that contains the plink files for the run. These should be map and ped files.",
                        dest="plink_dir", type=str, required=False)

    parser.add_argument("--no_carrier_file", "-nc",
                        help="This argument list the file path to the no_carriers.txt file",
                        dest="no_carriers_file", type=str, required=False)

    parser.add_argument("--ibd_programs", help="This argument list which ibd programs were used",
                        dest="ibd_programs", type=str, required=False)

    parser.add_argument("--file_suffix", "-fs", help="This argument list the suffix for the different ibd files because it could change",
                        dest="ibd_file_suffix", type=str, required=False)

    parser.add_argument("--var_input", "-vi", help="This argument list the path to the input file that list all the variants. This should be a csv or xlsx file",
                        dest="var_input_file", type=str, required=False)

    parser.add_argument("--mfile", "-mf", help="This argument list the directory containing all the map files for each chromosome",
                        dest="map_file", type=str, required=False)

    parser.add_argument("--minCM", "-m", help="This argument list the minimum cM threshold",
                        dest="min_CM", type=str, required=False, default=3)

    parser.add_argument("--drop_var", help="This functionality is used to drop variants from a file if needed to for some reason. This is passed into the searchPedigree function incase maybe a certain variant is too common and can be removed",
                        dest="drop_var", type=str, nargs="+")

    parser.add_argument("--format", help="This argument will enable several additional output files to be created that are easier for non python programs to interact with",
                        dest="compatible_format", type=bool, default=False)

    parser.add_argument("--allpair_files", "-a", help="This argument provides a path to the list of shared segments that can be used to form networks",
                        dest="allpair_file", type=str, default=False)

    parser.add_argument("--variant_file", help="This argument provides a path to a file that list all individuals that "
                                               "carry a specific variant",
                        dest="var_file", type=str, default=False)

    parser.add_argument("--var_of_interest", help="This argument passes a variant of interest that can filter down "
                                                  "dataframes when trying to draw networks.",
                        dest="var", type=str, nargs="+", default=False)

    parser.add_argument("--pop_info", help="This argument provides the file path to a file containing the population "
                                           "distribution of a dataset for each grid. This file should be a text file "
                                           "and at least contain two columns, where one column is 'Pop', the "
                                           "population code for each grid based on 1000 genomes, and then the second "
                                           "column is 'grid', which list the IIDs for each each grid.",
                        dest="pop_info", type=str, required=False)

    parser.add_argument("--pop_code", help="This argument can be a single population code or a list of population "
                                           "codes that someone is looking for. Population codes have to match the "
                                           "1000 genomes.",
                        dest="pop_code", type=str, required=False)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
