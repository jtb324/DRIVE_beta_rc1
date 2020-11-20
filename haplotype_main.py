import pandas as pd
import argparse
import multiprocessing as mp
from functools import partial
import os

from get_iid_pairs import get_pair_ids


def get_uniq_variants(dataframe: pd.DataFrame) -> list:
    '''This function gets the unique variants within a dataframe'''

    return list(set(dataframe.variant_id.tolist()))


def get_haplotype():


def parallelize(dataframe: str, workers: int, variant_list: list, output: str):
    '''This will be the main function to parallelize the analysis'''

    print("attempting to run in parallel...")

    # creating a manager que that can be used to line up how it write to the file
    manager = mp.Manager()

    # creating the que object for the haplotypes
    segment_que = manager.Queue()

    pool = mp.Pool(workers)

    header: str = f"pair_1\tpair_2\tchr\tprogram_used\tsegment\n"

    watcher = pool.apply_async(
        listener, (segment_que, "".join([output, "shared_segment_info.txt"]), header))

    func = partial(get_haplotype,

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

            m=que_object.get()

            if m == "kill":

                break

            output_file.write(m)
            output_file.flush()

def run(args):
    "function to run"

    for chunk in pd.read_csv(args.hap_file, sep="\t", chunksize=100000):

        uniq_var_list: list=get_uniq_variants(chunk)

        parallelize(chunk, args.workers, uniq_var_list, args.output)


def main():
    parser=argparse.ArgumentParser(
        description="")

    parser.add_argument("--bfile", help="This argument takes the input path for the binary file",
                        dest="bfile", type=str, required=True)

    parser.add_argument("--input", help="This argument takes the path for the haplotype_info.txt file",
                        dest="hap_file", type=str, required=True)

    parser.add_argument("--output", help="This argument takes the output path for the analysis",
                        dest="hap_file", type=str, required=True)

    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
