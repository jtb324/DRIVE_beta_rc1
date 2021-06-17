# This script is designed to take the output and provide something more like davids ideal format
import utility_scripts
import pandas as pd
import os

# Next three classes will be involved in reformatting the information in the allpair.
# txt file to the expected format in the confirmed_carriers.txt files
# The resulting confirmed_carriers.txt file needs the columns are
    # IID variant_id gene genotype confirmed_status chr

class Base_Reformatter:
    
    def __init__(self, output_path: str) -> None:
        """This is the base reformatter that will hold attributes that the other classes have in common"""
        self.output_path: str = os.path.join(output_path, "confirmed_carriers.txt")
    
    def write_to_file(self, dataframe: pd.DataFrame):
        """This function will write the dataframe to a file
        Parameters
        __________
        dataframe: pd.Dataframe
            dataframe that contains information about the iid, the variant or gene name, the genotype, the confirmed carrier status, and the chromosome number
        """

        # writing the dataframe to a file

        dataframe.to_csv(self.output_path, sep="\t", index=None)

    @staticmethod
    def get_file(file_list: list, identifier: str) -> str:
        """Function to get the specified file out of a list based on the identifier
        
        Parameters
        __________
        file_list : list
            This is a list of files that the specified file will be pulled out of
        
        identifier : str
            This is a string that can used to identify the file of interest
        
        Returns
        _______
        string
            returns a string that is the specified filepath that matches the 
            identifier
        """
        
        return [
            file for file in file_list if identifier in file
        ][0]
    
    @staticmethod
    def get_chr_num(pattern_list: list, file) -> str:
        """Function that will match the specified pattern within the provided file
        Parameters
        __________
        pattern_list : list
            list of regular expressions to match against
        
        file : str
            This string is the filename to match against
        
        Returns 
        _______
        str
            returns the chromosome number as a string
        """
        # [r"chr\d_", r"chr\d\d_"] pattern for the gene analysis

        return utility_scripts.match_chr(pattern_list, file)
    
    @staticmethod
    def get_confirmed_carriers(allpair_file: str, carrier_list: list) -> list:
        """Function to get the list of confirmed carriers from the allpairs.txt file
        Parameters
        __________
        allpair_file : str
            string that list the file path to the allpairs.txt file\
        
        Returns
        _______
        list
            returns a list of iids that where the iid is in the carrier list
        """
        # loading only the columns that contain pairs from the allpair file
        allpair_df: pd.DataFrame = pd.read_csv(allpair_file,
                                            sep="\t",
                                            usecols=["pair_1", "pair_2"])
        # filtering the dataframe for only cases where both pairs are in 
        # the carrier list
        allpair_df_carriers: pd.Dataframe = allpair_df[
            (allpair_df.pair_1.isin(carrier_list))
            & (allpair_df.pair_2.isin(carrier_list))]

        # getting all the grids for pair one
        pair_1_list: list = allpair_df_carriers.pair_1.values.tolist()

        # getting all the grids for pair two into a list
        pair_2_list: list = allpair_df_carriers.pair_2.values.tolist()

        # combining the two pair list into a set so that repeating values get dropped
        confirmed_carrier_set: set = set(pair_1_list + pair_2_list)

        # returning the list of individuals that are carriers for this variant and also share segments
        return list(confirmed_carrier_set)
    
    @staticmethod
    def get_confirmed_status(iid: str, confirmed_carriers: list) -> int:
        """Function that will return either a 1 or a zero depending on whether the iid is a confirmed carrier or not
        Parameters
        __________
        iid : str
            string that has the grids iid
        
        confirmed_carriers : list
            list of iids that were confirmed by the combined_ibd_pairs script
            
        Returns
        _______
            returns an integer of either 1 or 0 depending on whether the iid is a confirmed carrier or not
        """
        if iid in confirmed_carriers:
            return 1
        else:
            return 0

class Gene_Reformatter(Base_Reformatter):

    def __init__(self,carrier_df: pd.DataFrame, allpair_file_dir: str, plink_file_dir:str, no_carrier_file: str, output_path: str) -> None:
        """Function to reformat the allpair.txt files for the genotype approach
        Parameters
        __________
        carrier_file_dir : str
            This string is the path to the directory that has the allpair.txt files
        
        allpair_file_dir : str
            This string is the path to the directory that has the allpair.txt files
            
        plink_file_dir : str
            This string is the path to the directory that has the map and ped files 
            from running plink
        
        no_carrier_file : str
            string that list the filepath to the file that list variants that have no carriers
        
        output_path : str
            This is the path to write the output files to"""
        self.carrier_df: pd.DataFrame = carrier_df
        self.allpair_files: str = allpair_file_dir
        self.plink_files: str = plink_file_dir
        self.no_carrier_file: str = no_carrier_file
        super(Gene_Reformatter, self).__init__(output_path)
    
    @staticmethod
    def form_genotype_df(map_file_path: str, ped_file_path: str) -> pd.DataFrame:
        """ This function will form a dictionary of dictionaries where the outer key is the iid,
        the inner key is the variant and the inner value is the genotype string """

        # form three list for the iid, the variant, and the genotype string
        iid_list: list = []
        variant_list: list = []
        geno_list: list = []

        # opening the ped file
        with open(ped_file_path, "r") as ped_file:

            # iterating through row
            for row in ped_file:

                # split row
                row: list = row.split(" ")

                # getting only the genotype portion of the row

                geno_row: list = row[6:]

                # getting the iid
                iid: str = row[1]

                # using enumerate to iterate through each row of the file and get the index of the row
                with open(map_file_path, "r") as map_file:

                    for index, map_row in enumerate(map_file):

                        map_row: list = map_row.split("\t")

                        # getting the variant id
                        variant_id: str = map_row[1]

                        # getting the genotype index positions
                        start_indx: int = index * 2

                        second_indx: int = start_indx + 1

                        # getting the genotypes from the geno_row
                        allele_1: str = geno_row[start_indx]

                        allele_2: str = geno_row[second_indx]

                        # appending to all the list
                        iid_list.append(iid)

                        variant_list.append(variant_id)

                        geno_list.append("".join([allele_1, allele_2]).strip("\n"))

        # inserting values into the genotyped dict
        genotype_dict = {
            "IID": iid_list,
            "Variant": variant_list,
            "Genotype": geno_list
        }

        genotype_df = pd.DataFrame.from_dict(genotype_dict)

        return genotype_df

    @staticmethod
    def check_no_carrier(no_carrier_file: str, variant_id: str) -> int:
        """function that will check if the variant_id has no carriers"""

        # if this file is not present than the function needs to return 0
        if not os.path.isfile(no_carrier_file):

            return 0

        else:
            # loading the file into a dataframe
            no_carrier_df: pd.DataFrame = pd.read_csv(no_carrier_file,
                                                    sep="\t",
                                                    names=["variant", "chr"])

            # getting a list of all variants that had no carriers
            variant_list: list = no_carrier_df.variant.values.tolist()

            # figuring out if the current id is not in the no_carrier list
            if variant_id in variant_list:

                # returns 1 if it is
                return 1

            else:

                print(f"The variant {variant_id} is not a valid variants")

                # returns 0 if the variant is not found
                return 0
                
    def get_files(self) -> dict:
        """Function gather files into a dictionary
        Parameters
        __________
        analysis_type : str
            string telling what analysis type is used during the program.
            This value will be either gene, phenotype, "", or another

        allpair_files : str
            directory containing all the *allpair.txt files that have 
            the pairs where one iid is a carrier and whether or not the 
            second pair is also a carrier

        carrier_files : str
            directory containing all the *single_var_carrier.csv files

        plink_files : str
            directory containing the plink files that were output after running plink
            
        Returns
        _______
        dict
            a dictionary that has all of the list of the files of interest
        """
        # Getting list of the carrier files, the map files, the ped 
        # files and the allpair_files
        # carrier_files_list: list = utility_scripts.get_file_list(self.carrier_files, "*single_variant_carrier.csv")

        map_files_list: list = utility_scripts.get_file_list(self.plink_files, "*.map")

        ped_files_list: list = utility_scripts.get_file_list(self.plink_files, "*.ped")

        allpair_files_list: list = utility_scripts.get_file_list(self.allpair_files, "*allpair.txt")

        # returning the output as a dictionary
        self.file_dict: dict = {
            "carrier_df": self.carrier_df,
            "map_files": map_files_list,
            "ped_files": ped_files_list,
            "allpair_files": allpair_files_list
            }

    @staticmethod
    def get_genotype_list(geno_df: pd.DataFrame, carrier_list: list,
                    variant_id: str) -> list:
        """This function will get the genotypes out for each variant
        Parameters
        __________
        geno_df : pd.DataFrame
            dataframe that has all of the genotypes for each iid and each variant
        
        carrier_list : list
            list of iids that are identified as carriers for the variant 
        
        variant_id : str
            string that has the variant_id
        
        Returns 
        _______
        list
            returns a list of genotypes"""

        # subsetting the dataframe to only the list of identified carriers
        geno_df_subset: pd.DataFrame = geno_df[geno_df["IID"].isin(carrier_list)]

        geno_list = geno_df_subset[geno_df_subset.Variant == variant_id]["Genotype"].values.tolist()

        return geno_list

    def reformat(self):
        """Function that will write the information from the allpairs.txt file about 
        the confirmed carriers to the confirmed_carriers.txt file"""

        # creating a dictionary to keep track of parameters
        output_dict: dict = {
                "IID":[],
                "variant_id": [],
                "gene_name": [],
                "genotype": [],
                "confirmed_status":[],
                "chr":[]
            }
        
        # gathering all the files into a dictionary
        self.get_files()

        # Iterating through the ped files
        for file in self.file_dict["ped_files"]:

            # This gets the chromosome number for the files
            chr_num: str = utility_scripts.match_chr([r"chr\d_", r"chr\d\d_"], file)

            # removing the _ from the end of the chr_num
            chr_num = chr_num[:len(chr_num) - 1]
            
            map_file: str = Base_Reformatter.get_file(self.file_dict["map_files"], chr_num)

            genotype_df: pd.DataFrame = self.form_genotype_df(map_file, file)

            # getting the correct carrier_df subset based on chromosome based off of the chromosome
            car_df: pd.DataFrame = self.file_dict["carrier_df"]
            
            car_df_subset: pd.DataFrame = car_df[car_df.chr == chr_num]



            # getting a list of all the unique variants in the dataframe
            variant_list: list = list(set(car_df_subset["variant_id"].values.tolist()))

            # Iterating through each variant and then getting a list of carriers
            # for each variant
            for variant in variant_list:
                
                # getting a list of carriers from the carrier_df for the specific variant
                iid_list: list = car_df_subset[car_df_subset["variant_id"] == variant ]["iid"].values.tolist()

                # using list comprehension to get the allpair file for a specific chromosome and variant

                # There is an error if it does not find an allpair_file. Some of these files don't exist because there are no carriers
                try:
                    allpair_file: str = [
                        file for file in self.file_dict["allpair_files"]
                        if "".join([chr_num, "."]) in file and variant in file
                    ][0]

                except IndexError:

                    # returns either a 1 or 0 if the variant is in the no_carrier_file.txt list
                    carrier_int: int = self.check_no_carrier(self.no_carrier_file, variant)

                    if carrier_int == 1:

                        # print statement that explains that ther is no variant
                        print(f"The variant, {variant}, has no carriers")
                        # skips to the next iteration of the for loop
                        continue

                    elif carrier_int == 0:

                        print(f"The variant, {variant}, failed")

                        # writing the variant that failed to a file
                        with open("".join([self.output_path, "failed_variants.txt"]),
                                "a+") as file:

                            file.write(f"{variant}\n")

                        continue

                # getting a list of carriers that share segments and therefore are "confirmed"

                confirmed_carrier_list: list = Base_Reformatter.get_confirmed_carriers(allpair_file, iid_list)

                # need to get the genotype for each carrier_list

                # This will subset the dataframe for the carriers
                genotype_list: list = self.get_genotype_list(genotype_df,
                                    iid_list, variant[:-2])
                
                output_dict["IID"].extend(iid_list)
                output_dict["variant_id"].extend([variant]*len(iid_list))
                output_dict["gene_name"].extend(["N/A"]*len(iid_list))
                output_dict["genotype"].extend(genotype_list)
                output_dict["confirmed_status"].extend([Base_Reformatter.get_confirmed_status(iid, confirmed_carrier_list) for iid in iid_list])
                output_dict["chr"].extend([chr_num.strip(".")[-2:]]*len(iid_list))

        output_df: pd.DataFrame = pd.DataFrame.from_dict(output_dict)
        
        # Combining the dataframes
        self.write_to_file(output_df)
            

class Pheno_Reformatter(Base_Reformatter):

    def __init__(self, output_path: str, pheno_gmap: pd.DataFrame, pheno_carrier: pd.DataFrame, allpairs_file_dir: str) -> None:
        """Child class of the Base Reformatter that will reformat the allpairs.txt file for the phenotype analysis
        Parameters
        __________
        output_path : str
            string that the output file will be written to
        
        pheno_gmap : pd.Dataframe
            dataframe that has the gene name, the chromosome, and the gene start and 
            end location
        
        pheno_carrier : pd.DataFrame
            dataframe that has two columns where one is the IID and the other is the 
            gene name. This file list the suspected carriers based on phenotype
        
        allpairs_file_dir : str
            string that list the path to the directory that has the allpair.txt files
        """
        self.allpair_files: str = allpairs_file_dir
        self.pheno_gmap: pd.DataFrame = pheno_gmap
        self.pheno_carrier: pd.DataFrame = pheno_carrier
        super(Pheno_Reformatter, self).__init__(output_path)
    
    
    def get_iids(self, gene_name: str) -> list:
        """method to get the list of identified carriers for a certain gene
        Parameters
        __________
        gene_name : str
            string that has all the identified carriers for a certain gene
        
        Returns
        _______
        list
            returns a list of carriers for that gene
        """
        # filtering the dataframe for only the specific gene
        filtered_df: pd.DataFrame = self.pheno_carrier[self.pheno_carrier.gene == gene_name]

        # returning the list of suscepted carriers for that gene
        return filtered_df.IID.values.tolist()
    
    def get_gene_list(self):
        """Function to get the list of genes in the pheno_gmap file"""

        return self.pheno_gmap[0].values.tolist()
    
    

    def reformat(self):
        """Function that will reformat the allpair.txt file to the confimed_carriers.
        txt"""
        # getting a list of genes to iterate through
        gene_list: list = self.get_gene_list()

        # creating a dictionary to keep track of parameters
        output_dict: dict = {
                "IID":[],
                "variant_id": [],
                "gene_name": [],
                "genotype": [],
                "confirmed_status":[],
                "chr":[]
            }

        for gene in gene_list:

            # getting the list of supposed carriers for that gene
            iid_list: list = self.get_iids(gene)

            # get the allpair files from allpair file directory
            allpair_file_list: list = utility_scripts.get_file_list(self.allpair_files, "*allpair.txt")

            # getting the allpair.txt file
            allpair_file: str = Base_Reformatter.get_file(allpair_file_list, gene)

            # getting a list carriers which means that the iid is paired with
            # another iid that is in the iid_list
            confirmed_carriers: list = Base_Reformatter.get_confirmed_carriers(allpair_file, iid_list)
            
            chr_num: str = self.get_chr_num([r'chr\d\d.', r'chr\d.'], allpair_file)

            if len(chr_num.strip(".")) == 4:

                chr_num = "".join([chr_num[:3], "0", chr_num[-1]])
            



            output_dict["IID"].extend(iid_list)
            output_dict["variant_id"].extend(["N/A"]*len(iid_list))
            output_dict["gene_name"].extend([gene]*len(iid_list))
            output_dict["genotype"].extend(["N/A"]*len(iid_list))
            output_dict["confirmed_status"].extend([Base_Reformatter.get_confirmed_status(iid, confirmed_carriers) for iid in iid_list])
            output_dict["chr"].extend([chr_num.strip(".")[-2:]]*len(iid_list))
            
        # converting the dictionary to a dataframe
        confirmed_carrier_df: pd.DataFrame = pd.DataFrame.from_dict(output_dict)

        self.write_to_file(confirmed_carrier_df)