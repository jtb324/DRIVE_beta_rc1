import pandas as pd
import subprocess
import os


def get_pair_haplotype_str(binary_file: str, var_to_keep_file: str,
                           start_bp: str, end_bp: str, output_dir_path: str,
                           ibd_program: str, chr_num: int) -> str:
    '''This function will run plink to extract the haplotype of each pair'''

    # Creating an output path
    output_dir: str = "".join([output_dir_path, "haplotype_analysis/"])

    try:
        os.mkdir(output_dir)

    except FileExistsError:
        pass

    full_output_path: str = "".join(
        [output_dir, "plink_", ibd_program, "_", start_bp, "_", end_bp])

    subprocess.run(
        [
            "plink",
            "--bfile",
            binary_file,
            "--keep",
            var_to_keep_file,
            "from-bp",
            start_bp,
            "--to_bp",
            end_bp,
            "--chr"
            chr_num,
            "--recode",
            "--out",
            full_output_path,
        ],
        check=False,
    )

    return "".join([full_output_path, ".ped"])
