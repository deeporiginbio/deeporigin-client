.PHONY:  test test-github jupyter install

SHELL := /bin/bash
uname=$(shell uname -s)

repo=$(shell basename $(CURDIR))


chosen_tests=""
org_key="deeporigin"

test: 
	venv/bin/ruff format .
	venv/bin/ruff check --select I . --fix
	venv/bin/interrogate -c pyproject.toml -vv . -f 100 --omit-covered-files
	venv/bin/pytest -x -n "auto" --failed-first -k $(chosen_tests) --mock --org_key $(org_key) --dist loadfile
	venv/bin/pytest -x docs --markdown-docs --markdown-docs-syntax=superfences




# set up jupyter dev kernel
jupyter:
	-deactivate
	-yes | jupyter kernelspec uninstall $(repo)
	@source $(CURDIR)/venv/bin/activate && \
		python3 -m ipykernel install --user --name $(repo) && \
		deactivate

# install in a virtual env with all extras
install:
	@echo "Checking python version, venv and pip availability..."
	@python3 -c 'import sys; exit(0 if sys.version_info >= (3,10) else 1)' \
		|| (echo "You need Python 3.10 or higher. Please install it, e.g.: \n  sudo apt install python3.10 python3.10-venv python3.10-dev" && exit 1)	
	@(python3 -m venv --help >/dev/null 2>&1 && python3 -m pip --version >/dev/null 2>&1) \
		|| (echo "Python venv and/or pip are not available. Please install them, e.g.: \n  sudo apt install python3.10-venv python3.10-dev python3-pip" && exit 1)

	@echo "Installing deeporigin in editable mode in a venv..."
	@python3 -m venv venv
	@source $(CURDIR)/venv/bin/activate && \
		pip install --upgrade pip && \
	    pip install -v --no-cache-dir -e .[lint,test,dev,docs,plots,tools] && \
	    deactivate
	@-mkdir -p ~/.deeporigin
	@test -f ~/.deeporigin/deeporigin || ln -s $(CURDIR)/venv/bin/deeporigin ~/.deeporigin/deeporigin


docs-build:
	bash scripts/build_docs.sh

docs-serve:
	@echo "Serving docs locally..."
	@source $(CURDIR)/venv/bin/activate && \
	    mkdocs serve && \
	    deactivate

docs-deploy: 
	@echo "Deploying to live environment..."
	@source $(CURDIR)/venv/bin/activate && \
	    mkdocs gh-deploy && \
	    deactivate


test-github-live:
	pytest -v --ignore=tests/test_config.py --ignore=tests/test_context.py -n "auto" --dist loadfile


notebooks-html:
	@echo "Making marimo notebooks..."
	@source $(CURDIR)/venv/bin/activate && \
	rm -f docs/notebooks/*.html && \
	marimo export html notebooks/docking.py -o docs/notebooks/docking.html && \
	deactivate