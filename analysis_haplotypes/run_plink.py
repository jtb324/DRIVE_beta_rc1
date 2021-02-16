import pandas as pd
import subprocess
import os


def get_plink_haplotype_str(binary_file: str, start_bp: str, end_bp: str,
                            output_dir_path: str, ibd_program: str,
                            chr_num: str, variant_id: str,
                            network_id: str) -> str:
    '''This function will run plink to extract the haplotype of each pair'''

    # Creating an output path
    output_dir: str = output_dir_path
    try:
        os.mkdir(output_dir)

    except FileExistsError:
        pass

    full_output_path: str = "".join([
        output_dir, variant_id, "_", network_id, "_plink_", ibd_program, "_",
        str(start_bp), "_",
        str(end_bp)
    ])

    subprocess.run(
        [
            "plink",
            "--bfile",
            binary_file,
            "--from-bp",
            start_bp,
            "--to-bp",
            end_bp,
            "--chr",
            chr_num,
            "--recode",
            "--out",
            full_output_path,
        ],
        check=False,
    )

    return "".join([full_output_path, ".ped"])
