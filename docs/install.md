<<<<<<< HEAD
<<<<<<< HEAD
## Quickstart

Installation is as simple as running:
=======

=======
>>>>>>> feat(docs): refined docs
# Installation

## For users

Run the commands below to use your favorite package manager to install the Deep Origin CLI and Python client.

=== "pip"

    Run the following command:

    ```bash
    pip install deeporigin
    ```

    !!! warning
        We recommend installing this package into a virtual environment, using a tool such as venv, pyenv, Pipenv, Poetry, or conda.

=== "Pipenv"
    Run the following command:

    ```bash
    pipenv install deeporigin
    ```

=== "Poetry"
    Run the following command:

    ```bash
    poetry add deeporigin
    ```

## For developers

!!! warning "Developers only"
    If you intend to contribute this package, we recommend installing this package via the instructions below. If you merely intend to use this package, we recommend following the instructions above.

### Get the code
>>>>>>> fix(docs): wrote docs for _api

First, run the following command:
```bash
git clone git@github.com:formiclabs/deeporigin-client.git
```

### Install the development dependencies for this package

Second, install Python 3.9+ and
[make](https://www.gnu.org/software/make//).

To verify that both are installed, run the following commands:

```bash
python3 --version
# Python 3.12.3

make --version
# GNU Make 4.4.1
```

### Install this package locally

<<<<<<< HEAD
<<<<<<< HEAD
If you want to install this package to develop it further, run:
=======
You can then install the API locally by navigating to the directory you downloaded the code in, and running:
>>>>>>> fix(docs): wrote docs for _api
=======
Third, navigate to the directory you downloaded the code to, and run the following command:
>>>>>>> feat(docs): refined docs

```bash
make install
```

<<<<<<< HEAD
<<<<<<< HEAD
This will create a virtual environment and install `deeporigin` in editable mode into the environment. These commands will also add the CLI your PATH so that it can be run
without first activating the virtual environment.
=======
This installs in "editable" mode, so changes you made take effect
immediately. 

>>>>>>> fix(docs): wrote docs for _api
=======
This will install this package in an "editable" mode. In this mode, your changes will take effect
immediately.
>>>>>>> feat(docs): refined docs

### Running the tests for this package

Once installed, you can test this package by running the following command:

```bash
make test
<<<<<<< HEAD
<<<<<<< HEAD
```
=======
```
>>>>>>> fix(docs): wrote docs for _api
=======
```
>>>>>>> feat(docs): refined docs
