# Signifies our desired python version
# Makefile macros (or variables) are defined a little bit differently than traditional bash, keep in mind that in the Makefile there's top-level Makefile-only syntax, and everything else is bash script syntax.
PYTHON = python3


# Defines the default target that `make` will to try to make, or in the case of a phony target, execute the specified commands
# This target is executed whenever we just type `make`
.DEFAULT_GOAL = help


# The @ makes sure that the command itself isn't echoed in the terminal
help:
	@echo "------------------HELP------------------"
	@echo "To test the project type make test"
	@echo "To run the project type make run"
	@echo "To build the project type make build"
	@echo "To clean the __pycache__ type make clean"
	@echo "----------------------------------------"

unittest:
	@echo "Running tests to ensure the integrity of all modules..."
	@cd test/ && ${PYTHON} -m pytest -v -s


build:
	@echo "Building the executable file"
	@echo "Make sure to have activated the conda environment"
	@echo please delete any previous build or dist directories you may have and the *.spec file
	@pyinstaller --onefile ./drive/DRIVE.py

clean:
	@find . | grep -E "(__pycache__|\.pytest_cache)" | xargs rm -rf

remove_prior_build:
	@echo "removing the files from the prior build"
	@rm -r build/ dist/ DRIVE.spec
