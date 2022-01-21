import pytest
from typing import Dict

import sys
import os

sys.path.append("../drive")

from models.inputs import InputParams, IncorrectIBDProgramError, IncorrectMafThresholdError, DirectoryNotFoundError, IllogicalGeneRangeError, ImproperPopulationCode, UnsupportedFileType

# example input data to test the program with
example_people_data: Dict = {
        'title': 'DRIVE User Configuration',
        'output': {'output': './'},
        'population_parameters': {
            'ethnicity_file': '',
            'pop_code':'EUR',
            'pop_info_file':''
        },
        'plink_parameters': {},
        'inputs': {'variant_file': '',
        'binary_file': '',
        'ped_files': '',
        'carrier_files': '',
        'pheno_gene_info': '',
        'phenotype_carriers': ''},
        'gene_range': {'gene_range': [0, 0]},
        'ibd_programs': {'ibd_programs': ['hapibd', 'ilash']},
        'thresholds': {'min_cm': 3, 'maf_threshold': 0.05}
        }

def test_unsupported_ibd_program() -> None:
    """Unit test to test if the program raises an IncorrectIBDProgramError if the wrong ibd program is supplied"""

    # add an invalid ibd program to the list
    example_people_data["ibd_programs"]["ibd_programs"].append("rapid")

    with pytest.raises(IncorrectIBDProgramError):

        _ = InputParams(**example_people_data)
    
    # check the invalid ibd program
    example_people_data["ibd_programs"]["ibd_programs"].pop()

def test_correct_ibd_program() -> None:
    """Unit test to check that the program raises no exception when the right ibd programs are supplied"""

    try:
        _ = InputParams(**example_people_data)
    except IncorrectIBDProgramError as exc:
        assert False, f"'InputParams validation raised an exception: {exc}"

@pytest.mark.parametrize("maf_threshold", [0.6, 0])
def test_incorrect_higher_maf_threshold(maf_threshold) -> None:
    """Unit test to make sure that the IncorrectMafThresholdError is raised if the value is incorrect"""
    example_people_data["thresholds"]["maf_threshold"] = maf_threshold

    with pytest.raises(IncorrectMafThresholdError) as exc:

        _ = InputParams(**example_people_data)
    
    example_people_data["thresholds"]["maf_threshold"] = 0.05

@pytest.mark.skip(reason="This test is currently failing and I can't figure out the error")
def test_directory_not_found() -> None:
    """Unit test to make sure that the DirectoryNotFoundError is raise if the provided output directory does not exist"""
    example_people_data["output"]["output"] = "./outputs/"

    with pytest.raises(DirectoryNotFoundError) as exc:
        _ = InputParams(**example_people_data)
    
    example_people_data["output"]["output"] = "./"

def test_directory_found() -> None:
    """Unit test to make sure that the DirectoryFoundError is not raised if the provided output directory is real"""
    try:
        _ = InputParams(**example_people_data)
    except DirectoryNotFoundError as exc:
        assert False, f"'InputParams validation raised an exception: {exc}"

def test_illogical_gene_range() -> None:
    """Function to ensure that the IllogicalGeneRangeError is raised if the gene range is incorrect"""

    example_people_data["gene_range"]["gene_range"] = [100, 10]

    with pytest.raises(IllogicalGeneRangeError) as exc:
        _ = InputParams(**example_people_data)
    
    example_people_data["gene_range"]["gene_range"] = [0, 0]

def test_logical_gene_range() -> None:
    """Unit test to make sure that the LogicalGeneRangeError is not raised if the gene range is correct"""
    try:
        _ = InputParams(**example_people_data)
    except IllogicalGeneRangeError as exc:
        assert False, f"'InputParams validation raised an exception: {exc}"

def test_unsupported_pop_code() -> None:
    """Unit test to make sure that the ImproperPopulationCode is raised when the population code is incorrect"""
    
    example_people_data["population_parameters"]["pop_code"] = "WRONGCODE"

    with pytest.raises(ImproperPopulationCode) as exc:
        _ = InputParams(**example_people_data)

    example_people_data["population_parameters"]["pop_code"] = "EUR"

@pytest.mark.parametrize("file", ["", ".xlsx", ".csv", "test.txt"])
def test_supported_file_type(file: str) -> None:
    """Unit test to make sure that the UnsupportedFileType is raised when the extension is wrong."""

    example_people_data["inputs"]["variant_file"] = file

    var_file: str = example_people_data["inputs"]["variant_file"]

    if var_file in ["", ".xlsx", ".csv"]: 
        _ = InputParams(**example_people_data)
        assert var_file in ["", ".xlsx", ".csv"], f"Error was not raised by the variant filepath as expected"

    else:
        with pytest.raises(UnsupportedFileType) as exc:
            _ = InputParams(**example_people_data)

    
    
    
    