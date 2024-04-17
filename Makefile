.PHONY:  test

SHELL := /bin/bash

repo=$(shell basename $(CURDIR))
client="mock"

test: 
	source $(CURDIR)/venv/bin/activate && \
		python3 -m coverage run -m pytest -k "test_cli or managed_data" --client $(client) && \
		python3 -m coverage html && \
		deactivate
	-open htmlcov/index.html
	


# set up jupyter dev kernel
jupyter:
	-deactivate
	-yes | jupyter kernelspec uninstall $(repo)
	source $(CURDIR)/venv/bin/activate && \
		python3 -m ipykernel install --user --name $(repo) && \
		deactivate

# install in a virtual env
install:
	python3 -m venv venv
	source $(CURDIR)/venv/bin/activate && \
		pip install -e .[test,jupyter] && \
		deactivate


