# README for the gnomad_webscrapper directory:

## Goal

The goal of this directory was to automate the process of finding minor allele frequencies in gnomad

## Files within directory:

- _gnomad_scrapper.py_: This file contains the main functions used in the webscrapper as well as the CLI

- _webdriver.py_: This file contains a function that determines the appropriate webdriver to use or exits the system if the user tries to use a non-supported driver

- _requirements.txt_: This file describes the modules required to run the program

## Directories to remember:

- _shell_scripts_: This directory contains example scripts I used in test runs

## Running the program:

create a simple shell script following the examples. The program has four current flags:

_-v_: This provides the file path to the xlsx/csv/txt file containing a column named RS Name which contains all the rsIDs of the variants of interest

_-o_: This argument provides an name for the output file

_-b_: This argument indicates which browser is being used. Either firefox or chrome has to be used to get the appropriate webdriver

_-t_: This argument specifies the wait time that the webdriver will wait to receive a http response from the gnomad website. This value defaults to 20 seconds. If a timeoutexception is thrown then increase the time.
