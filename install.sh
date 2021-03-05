#!/bin/bash
###############################################################################
#Purpose:
#
#Install script to create the necessary conda environment for the DRIVE program
###############################################################################
#The script expects one argument to be passed along with the script name.
#This argument will be whatever the user wants to name the conda environment.

#Creating colors to be used during the program
RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m'
#This next block checks to make sure that argument is passed
echo "checking for the correct number of arguments"
if (($# != 1)); then
    echo -e "${RED}ERROR:${NC}--Illegal number of parameters passed"
    echo -e "${YELLOW}FIX:${NC}--The script was expecting two  following the format: install.sh '{conda environment name}'"
else 
    echo -e "${GREEN}SUCCESS:${NC}--Correct number of arguments identified"
    ENV_NAME="$1"

    #checking to see if a conda environment is installed
    echo "checking to see if anaconda is within the users PATH"
    if command -v conda &> /dev/null;then
        echo -e "${GREEN}SUCCESS:${NC}--found anaconda installation"
        echo "Creating the conda virtual environment"
        conda create -n "$ENV_NAME" python=3.7
        echo -e "${GREEN}SUCCESS:${NC}--conda environment created. Now activating the environment..."
        source activate "$ENV_NAME"
        echo "Fetching the newest version of pip..."
        conda install pip
        echo -e "${GREEN}SUCCESS:${NC}--newest version of pip successfully installed"
        echo "installing all required packages from the requirements.txt file"
        pip install -r requirements.txt
        echo -e "${GREEN}SUCCESS:${NC}--installed all necessary packages for the DRIVE program."
        echo "Activate the new conda environment using conda activate $ENV_NAME"
    else
        echo -e "${RED}ERROR:${NC}--no anaconda installation was found within the users PATH"
        echo -e "${YELLOW}FIX:${NC}--please visit https://docs.continuum.io/anaconda/install/ to download and install the anaconda package manager and make sure to add the program to you \$PATH"
        echo -e "${YELLOW}FIX:${NC}--or visit the miniconda documentation at  https://docs.conda.io/en/latest/miniconda.html for a more lightweight installation"
    fi
fi
