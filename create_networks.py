import logging
import re

from network_creator_class import Network_Img_Maker


def create_networks(segments_file: str, variant_file: str, ind_in_network_dict: dict, variant_of_interest: str, output_path: str) -> dict:

    logger = logging.getLogger(output_path+'/drawing_networks.log')

    network_drawer = Network_Img_Maker(
        segments_file, variant_file, output_path, logger)

    allpair_file_list = network_drawer.gather_allpair_files(
        segments_file, "*.allpair.new.txt")

    # write allpair file variants to a dictionary
    # Each key will be a variant and the value will be the proper file
    # have to get index of 1st and 3rd "_"
    allpair_var_dict = dict()

    for file in allpair_file_list:

        name_list = file.split("_")
        variant_id = "".join([name_list[1], name_list[2]])
        print(variant_id)

        # also have to get chromosome number
        match = re.search(r'chr\d\d', file)

        chr_num = match.group(0)

        print(chr_num)
        # build dictionary

        allpair_var_dict[(variant_id, chr_num)] = file

    carrier_list = network_drawer.gather_allpair_files(
        variant_file, "*.single_variant_list.csv")

    if not allpair_var_dict:
        print("There was an error that occurred when trying to gather all the different variants on all of the chromosomes")
        return None

    # looping through all the variants in the dictionary
    for item in allpair_var_dict.items():

        variant = item[0][0]

        chr_num = item[0][1]

        segment_file = item[1]

        print(segment_file)

        iid_list = network_drawer.isolate_variant_list(variant)

    # Need to adjust the rest of the file so that it can run for all the variants.

    # this line loads the provided dataframe using a tab separator with no header value. it also skips the first row
    # which only contains information about the chr, the number of pairs. The second row and onwards list all the pairs
        loaded_segments_file = network_drawer.check_file_exist(
            separator="\t", header_value=None, skip_rows=1)

        # This just renames the first row to Pairs which can be used to identify the column later
        loaded_segments_file = loaded_segments_file.rename(columns={
                                                           0: "Pairs"})

        # This just drops any empty rows in the dataframe
        if loaded_segments_file["Pairs"].isnull().any():
            loaded_segments_file = network_drawer.drop_empty_rows(
                loaded_segments_file)

        network_drawer_df = network_drawer.isolate_ids(
            loaded_segments_file, iid_list)

        carrier_in_network_dict = network_drawer.carriers_in_network(
            iid_list, network_drawer_df, ind_in_network_dict)

        network_drawer.draw_networks(network_drawer_df)

    return carrier_in_network_dict
