import argparse
import pandas as pd

from IBD_compare_output_classes import Output_Comparer

# run function


def run(args):
    print("running")

    variant_info_df = pd.read_csv(args.var_info, sep=" ", header=None, names=[
        "output_file_name", "variant_bp", "variant_id"])

    output_file_name_list = variant_info_df.output_file_name.values.tolist()

    variant_bp_list = variant_info_df.variant_bp.values.tolist()

    variant_id_list = variant_info_df.variant_id.values.tolist()

    for variant_info_tuple in zip(output_file_name_list, variant_bp_list, variant_id_list):

        output = "".join([args.output, variant_info_tuple[0]])

        var_pos = variant_info_tuple[1]

        var_of_interest = str(variant_info_tuple[2])

        x = Output_Comparer(output, var_pos, var_of_interest, args.input_dir)

        x.check_arguments(args.input)

        file_dict = x.create_file_dict(args.input, var_of_interest)

        first_line_dict = x.read_first_line(file_dict)

        allcomb, combtab = x.define_input_combinations(file_dict)

        x.write_to_file(allcomb, combtab, first_line_dict)


def main():
    parser = argparse.ArgumentParser(
        description="This script compares the output of the Shared_segment_Formatter_CLI.py to each other and combines them into one")

    parser.add_argument("-o", '--output', help="This argument just list the directory that files are being written to",
                        dest="output", type=str, required=True)

    parser.add_argument("-d", "--input_dir", help="This argument list the directory that the input files are found in.",
                        dest="input_dir", type=str, required=True)

    parser.add_argument("-i", '--input', help="This argument just list the different files as input",
                        dest="input", type=str, nargs="+", required=True)

    parser.add_argument("-v", "--var_info", help="This argument list the directory to a file that has three columns,'output_file_name', 'variant_bp', 'variant_id'. Respectively, these columns list the filenames for the output file for each variant, the location of the variant, and the variant id.",
                        dest="var_info", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
