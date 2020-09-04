import argparse
import pandas as pd

from IBD_compare_output_classes import Output_Comparer, Gather_IBD_Output

# run function


def run(args):
    print("running")
    # First step gets all the file paths
    ibd_output_collector = Gather_IBD_Output(
        args.input_dir, args.ibd, args.map_file)

    map_file_list = ibd_output_collector.get_map_files()

    ibd_files_dict = ibd_output_collector.return_dict()

    # Needs to be a function that can get the variant bp and make the output name
    for chr_variant_tuple in ibd_files_dict.keys():

        file_list = list(ibd_files_dict[chr_variant_tuple])

        chr_num = chr_variant_tuple[0]
        # Removing the ".""
        chr_num = chr_num[:len(chr_num)-1]

        # adding a "_"
        chr_num = "".join([chr_num[1:], "_"])

        variant = chr_variant_tuple[1]

        map_file = [
            map_file for map_file in map_file_list if chr_num in map_file][0]

        print(map_file)
        print(chr_num)
        variant_bp = ibd_output_collector.get_variant_bp(
            variant, args.var_list, map_file)

        if variant_bp == -1:
            print("There was no base position found for the variant")
            continue

        output_file_name = "".join(["IBD_", variant, "_", chr_num])

    # for variant_info_tuple in zip(output_file_name_list, variant_bp_list, variant_id_list):

        full_output_path = "".join([args.output, output_file_name])

        output_comparer = Output_Comparer(full_output_path, int(
            variant_bp), variant, args.input_dir, file_list)

        output_comparer.check_arguments(file_list)

        file_dict = output_comparer.create_file_dict(file_list, variant)

        first_line_dict = output_comparer.read_first_line(file_dict)

        allcomb, combtab = output_comparer.define_input_combinations(file_dict)

        output_comparer.write_to_file(allcomb, combtab, first_line_dict)


def main():
    parser = argparse.ArgumentParser(
        description="This script compares the output of the Shared_segment_Formatter_CLI.py to each other and combines them into one")

    parser.add_argument("-o", '--output', help="This argument just list the directory that files are being written to",
                        dest="output", type=str, required=True)

    parser.add_argument("-d", "--input_dir", help="This argument list the directory that the input files are found in.",
                        dest="input_dir", type=str, required=True)

    parser.add_argument("-c", '--var_list', help="This argument just list the filepath to the original file containing the variants and information about each variant",
                        dest="var_list", type=str, required=True)

    parser.add_argument("-s", '--ibd', help="This argument just list the different programs used for the ibd program",
                        dest="ibd", type=str, nargs="+", required=True)

    parser.add_argument("-m", '--map_file', help="This argument just passes the directory to the map files",
                        dest="map_file", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
