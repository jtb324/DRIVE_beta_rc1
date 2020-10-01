# This script is going to be responsible for  figuring out which webdriver to use
from selenium import webdriver
import sys


def determine_webdriver(browser: str):

    if browser == "firefox":

        driver = webdriver.Firefox()

    elif browser == "chrome":

        driver = webdriver.Chrome()

    else:

        print("Please use one of the two supported browsers: Firefox or Google Chrome")

        sys.exit(1)

    return driver
