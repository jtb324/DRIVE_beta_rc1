import pandas as pd
import argparse
import multiprocessing as mp
from functools import partial
import os

from get_iid_pairs import get_pair_ids
from plink_haplotype import Plink_Haplotype_Getter


def get_uniq_variants(dataframe: pd.DataFrame) -> list:
    '''This function gets the unique variants within a dataframe'''

    return list(set(dataframe.variant_id.tolist()))


def create_filter_file(pair_1: str, pair_2: str, output: str) -> str:
    '''This function will create the filter file that is need to get only the IIDs of interest. It will return a string of the path to the file'''

    with open("".join([output, "filter_file_", pair_1, "-", pair_2, ".txt"]), "w") as iid_file:

        iid_file.write(f"{pair_1}\t{pair_1}")

        iid_file.write(f"{pair_2}\t{pair_2}")

        iid_file.close()

    return "".join([output, "filter_file_", pair_1, "-", pair_2, ".txt"])


def feed_to_plink_haplotype(binary_file: str, chr_num: str, start_bp: str, end_bp: str, output_path: str, filter_file_path: str) -> str:
    '''This function will run the plink file to get the genotype. It forms a ped file for output and then will return the path to the file'''

    plink_object = Plink_Haplotype_Getter(
        binary_file, chr_num, start_bp, end_bp, output_path, filter_file_path)

    output: str = plink_object.run_plink()

    return "".join([output, ".ped"])


def get_haplotype_str(file_path: str) -> list:
    '''This function will return the haplotype string'''

    haplotype_list: list = []

    with open(file_path, "r") as ped_file:

        for row in ped_file:

            haplotype_str: str = row.split(" ", maxsplit=6)[7]

            haplotype_list.append(haplotype_str.strip("\n"))

    return haplotype_list


def get_haplotype(dataframe: pd.DataFrame, binary_file: str, output_path: str, que_object: str, variant_id: str):
    '''This function will get the haplotype for each pair'''

    dataframe_subset: pd.DataFrame = dataframe[dataframe.variant_id == variant_id]

    # iterating through each row in the dataframe
    for row in dataframe_subset.itertuples():

        print(row)
        # Breaking up the row into separate parts
        pair_1: str = row[1]
        pair_2: str = row[2]
        chr_num: str = str(row[3])
        network_id: str = str(row[5])
        hapibd_start: str = str(row[6])
        hapibd_end: str = str(row[7])
        ilash_start: str = str(row[9])
        ilash_end: str = str(row[10])

        filter_file_path: str = create_filter_file(pair_1, pair_2, output_path)

        full_output_path: str = "".join(
            [output_path, "plink_haplotype_files", pair_1, "-", pair_2, "_haplotypes"])

        # getting info for hapibd
        hapibd_output_path: str = feed_to_plink_haplotype(
            binary_file, chr_num, hapibd_start, hapibd_end, full_output_path, filter_file_path)

        ilash_output_path: str = feed_to_plink_haplotype(
            binary_file, chr_num, ilash_start, ilash_end, full_output_path, filter_file_path)


def parallelize(dataframe: str, workers: int, variant_list: list, output: str, binary_file: str):
    '''This will be the main function to parallelize the analysis'''

    print("attempting to run in parallel...")

    # creating a manager que that can be used to line up how it write to the file
    manager = mp.Manager()

    # creating the que object for the haplotypes
    segment_que = manager.Queue()

    pool = mp.Pool(workers)

    header: str = f"pair_1\tpair_2\tchr\tnetwork_id\tvariant_id\tprogram_used\tsegment\n"

    watcher = pool.apply_async(
        listener, (segment_que, "".join([output, "shared_segment_info.txt"]), header))

    func = partial(get_haplotype, binary_file, dataframe, output)

    pool.map(func, variant_list)

    segment_que.put("kill")

    pool.close()

    pool.join()


def listener(que_object, output: str, header: str):
    '''This function will listen to the que and then write the element of the que to a file'''

    # opening the output file to write to
    with open(output, "a+") as output_file:

        # checking if the file size is zero
        if os.path.getsize(output) == 0:

            output_file.write(header)

        while 1:

            m = que_object.get()

            if m == "kill":

                break

            output_file.write(m)
            output_file.flush()


def run(args):
    "function to run"

    for chunk in pd.read_csv(args.hap_file, sep="\t", chunksize=100000):

        uniq_var_list: list = get_uniq_variants(chunk)

        parallelize(chunk, args.workers, uniq_var_list,
                    args.output, args.bfile)


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument("--bfile", help="This argument takes the input path for the binary file",
                        dest="bfile", type=str, required=True)

    parser.add_argument("--input", help="This argument takes the path for the haplotype_info.txt file",
                        dest="hap_file", type=str, required=True)

    parser.add_argument("--output", help="This argument takes the output path for the analysis",
                        dest="hap_file", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
