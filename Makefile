.PHONY: test test-github jupyter install coverage docs-build docs-serve docs-deploy test-github-live notebooks-html

SHELL := /bin/bash
uname=$(shell uname -s)

repo=$(shell basename $(CURDIR))

# determines how tests are run. if "mock", tests are run
# using mocked responses. if "default", tests run against
# a live instance 
client="mock"
chosen_tests=""
responses="pass"

help:       ## Show available commands
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z0-9_-]+:.*?## / { \
		if ($$2 ~ /^─+/) print "\n" $$2; \
		else printf "  \033[36m%-22s\033[0m %s\n", $$1, $$2 \
	}' $(MAKEFILE_LIST)
%:

test:       ## Run tests with coverage and linting
	ruff check .
	ruff check --select I --fix
ifeq ($(client), "mock")
	$(eval n_workers=1)
else 
	$(eval n_workers="auto")
endif 
	@source $(CURDIR)/venv/bin/activate && \
	interrogate -c pyproject.toml -vv . -f 100 --omit-covered-files && \
	python3 -m coverage run --source="src" -m pytest -x -n $(n_workers) --failed-first -k $(chosen_tests) --client $(client) --responses $(responses) --dist loadfile && \
	python3 -m coverage html && \
	deactivate

coverage:   ## Run coverage and open HTML report
	@source $(CURDIR)/venv/bin/activate && \
	python3 -m coverage run -m pytest -x --client $(client)  && \
	python3 -m coverage html && \
	open htmlcov/index.html && \
	deactivate

jupyter:    ## Set up Jupyter kernel for this repo
	-deactivate
	-yes | jupyter kernelspec uninstall $(repo)
	@source $(CURDIR)/venv/bin/activate && \
		python3 -m ipykernel install --user --name $(repo) && \
		deactivate

install:    ## Install dependencies in editable mode inside venv
	@echo "Installing deeporigin in editable mode in a venv..."
	@python3 -m venv venv
	@source $(CURDIR)/venv/bin/activate && \
		pip install --upgrade pip && \
	    pip install --no-cache-dir -e .[lint,test,dev,docs,plots,tools] && \
	    deactivate
	@-mkdir -p ~/.deeporigin
	@test -f ~/.deeporigin/deeporigin || ln -s $(CURDIR)/venv/bin/deeporigin ~/.deeporigin/deeporigin

docs-build: ## Build documentation
	bash scripts/build_docs.sh

docs-serve: ## Serve docs locally with MkDocs
	@echo "Serving docs locally..."
	@source $(CURDIR)/venv/bin/activate && \
	    mkdocs serve && \
	    deactivate

docs-deploy:## Deploy documentation to GitHub Pages
	@echo "Deploying to live environment..."
	@source $(CURDIR)/venv/bin/activate && \
	    mkdocs gh-deploy && \
	    deactivate

test-github: ## Run GitHub CI-style test with default client
	python3 -m coverage run -m pytest --client $(client)

test-github-live: ## Run GitHub CI-style test without mocks and ignoring config/context tests
	python3 -m coverage run -m pytest --ignore=tests/test_config.py --ignore=tests/test_context.py --client default -n "auto" --dist loadfile

notebooks-html: ## Export Marimo notebooks to HTML
	@echo "Making marimo notebooks..."
	@source $(CURDIR)/venv/bin/activate && \
	rm -f docs/notebooks/*.html && \
	marimo export html notebooks/docking.py -o docs/notebooks/docking.html && \
	deactivate