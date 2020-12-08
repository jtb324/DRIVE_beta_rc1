import re


def get_chr_num(file: str) -> str:
    '''This function will get the chr_num from the file name'''

    match = re.search(r'chr\d\d_', file)

    # find chromosome number
    if match:

        chr_num = match.group(0)

    else:
        match = re.search(r'chr\d_', file)

        chr_num = match.group(0)

    chr_num = chr_num.strip(".").strip("_")

    return chr_num


def get_alt_chr_num(chr_num: str) -> str:
    '''This function is used to convert the chromosome numbers 1-9 from
    a format of .chr01_ to .chr1_ incase the file needs that. Some older
    files may use this format and this will allow the program to keep
    running'''

    return "".join(chr_num.split("0"))
