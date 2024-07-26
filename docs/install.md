# Installation

## For users

Run the commands below to use your favorite package manager to install the Deep Origin CLI and Python client.

Run the following command:


!!! warning
    We recommend installing this package into a virtual environment. Some tools create their own virtual environments, like poetry.


=== "pip"

    ```bash
    pip install -q deeporigin
    ```

=== "pixi"

    ```bash
    pixi add deeporigin
    ```

=== "pipx"

    ```bash
    pipx install deeporigin
    ```

=== "pipenv"

    ```bash
    pipenv install deeporigin
    ```

=== "poetry"

    ```bash
    poetry add deeporigin
    ```

=== "uv"

    ```bash
    uv pip install deeporigin
    ```

=== "micromamba"

    ```bash
    micromamba install -c https://repo.prefix.dev/deeporigin-public deeporigin
    ```

=== "flit"

    ```bash
    flit install deeporigin
    ```

## Supported Python versions

`deeporigin` is [tested](https://github.com/deeporiginbio/deeporigin-client/actions/workflows/main.yml) against these versions of Python using GitHub Actions:

- 3.9
- 3.10
- 3.11
- 3.12

## For developers

!!! warning "Developers only"
    If you intend to contribute this package, we recommend installing this package via the instructions below. If you merely intend to use this package, we recommend following the instructions above.

### Get the code

First, run the following command:
```bash
git clone git@github.com:deeporiginbio/deeporigin-client.git
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
