import subprocess
import pandas as pd
import os
import sys


class Plink_Haplotype_Getter:

    def __init__(self, bfile: str, chr_num: int, start_bp: str, end_bp: str, output: str, filter_iid_path: str):
        self.binary_file: str = bfile
        self.chr_num: str = str(chr_num)
        self.start_bp: str = str(start_bp)
        self.end_bp: str = str(end_bp)
        self.output_path: str = output
        self.filter_iid_path: str = filter_iid_path

    def run_plink(self) -> str:
        '''This function will run PLINK to get the haplotype'''

        subprocess.run(["plink",
                        "--bfile",
                        self.binary_file,
                        "--chr",
                        self.chr_num,
                        "--from-bp",
                        self.start_bp,
                        "--to-bp",
                        self.end_bp,
                        "--keep",
                        self.filter_iid_path,
                        "--out",
                        self.output_path,
                        "--recode",
                        ])

        return self.output_path
