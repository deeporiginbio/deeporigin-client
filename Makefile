.PHONY:  test test-github jupyter install

SHELL := /bin/bash
uname=$(shell uname -s)

repo=$(shell basename $(CURDIR))

# determines how tests are run. if "mock", tests are run
# using mocked responses. if "default", tests run against
# a live instance 
client="mock"
chosen_tests=""


test: 
ifeq ($(client), "mock")
	source $(CURDIR)/venv/bin/activate && \
		interrogate -c pyproject.toml -v . -f 100 && \
		python3 -m coverage run -m pytest -x --failed-first -k $(chosen_tests) --client $(client) && \
		python3 -m coverage html && \
		deactivate
	if [ "$(uname)" = "Linux" ]; then \
		xdg-open htmlcov/index.html; \
	else \
		open htmlcov/index.html; \
	fi
else 
	source $(CURDIR)/venv/bin/activate && \
		interrogate -c pyproject.toml -v . -f 100 && \
		python3 -m coverage run -m pytest --failed-first -k $(chosen_tests) --client $(client) -n auto && \
		python3 -m coverage html && \
		deactivate
	if [ "$(uname)" = "Linux" ]; then \
		xdg-open htmlcov/index.html; \
	else \
		open htmlcov/index.html; \
	fi
endif 

# set up jupyter dev kernel
jupyter:
	-deactivate
	-yes | jupyter kernelspec uninstall $(repo)
	source $(CURDIR)/venv/bin/activate && \
		python3 -m ipykernel install --user --name $(repo) && \
		deactivate

# install in a virtual env with all extras
install:
	@echo "Installing deeporigin in editable mode in a venv..."
	@python3 -m venv venv
	@source $(CURDIR)/venv/bin/activate && \
		pip install --upgrade pip && \
	    pip install -e .[test,jupyter,docs] && \
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

test-github:
	python3 -m coverage run -m pytest -k $(chosen_tests) --client $(client)
