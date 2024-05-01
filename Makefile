.PHONY:  test test-github jupyter install

SHELL := /bin/bash

repo=$(shell basename $(CURDIR))

# determines how tests are run. if "mock", tests are run
# using mocked responses. if "default", tests run against
# a live instance 
client="mock"
chosen_tests=""

test: 
	source $(CURDIR)/venv/bin/activate && \
		python3 -m coverage run -m pytest --failed-first -k $(chosen_tests) --client $(client) && \
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

# install in a virtual env with all extras
install:
	@python3 -m venv venv
	@source $(CURDIR)/venv/bin/activate && \
	    pip install -e .[test,jupyter,docs] && \
	    deactivate
	@-mkdir -p ~/.deeporigin
	if [ ! -f ~/.deeporigin/deeporigin ]; then
		@-ln -s $(CURDIR)/venv/bin/deeporigin ~/.deeporigin/deeporigin
	fi
	@echo 'export PATH="$$HOME/.deeporigin:$$PATH"' >> ~/.bash_profile

build-docs:
	bash scripts/build_docs.sh

serve-docs:
	@echo "Serving docs locally..."
	@source $(CURDIR)/venv/bin/activate && \
	    mkdocs serve && \
	    deactivate

test-github:
	python3 -m coverage run -m pytest -k $(chosen_tests) --client $(client)


deploy-docs: 
	@echo "Deploying to live environment..."
	@source $(CURDIR)/venv/bin/activate && \
	    mkdocs gh-deploy && \
	    deactivate