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

Third, navigate to the directory you downloaded the code to, and run the following command:

```bash
make install
```

This will install this package in an "editable" mode. In this mode, your changes will take effect
immediately.

### Running the tests for this package

Once installed, you can test this package by running the following command:

```bash
make test
```
