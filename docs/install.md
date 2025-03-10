# Installation


## On a Deep Origin workstation

No installation needed!

The Deep Origin CLI and Python client are installed on every workstation.

## On your local computer

Run the command below for your favorite package manager to install the Deep Origin CLI and Python client.

=== "pip"

    ```bash
    pip install deeporigin
    ```

=== "pixi"

    ```bash
    pixi add --pypi deeporigin
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

=== "flit"

    ```bash
    flit install deeporigin
    ```

## Using `uv` to set up `deeporigin` on your computer

We recommend using the following:

```bash
# make a folder for your project
# this can be any name (except deeporigin)
mkdir do-client
cd do-client

# install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

uv init
uv python install
uv add deeporigin

# if you're using Deep Origin tools (like FEP), use:
uv add deeporigin --extra tools

# run jupyter using this environment
uv run --with jupyter jupyter lab
```

## For developers

!!! warning "Developers only"
    If you intend to contribute this package, we recommend installing this package via the instructions below. If you simply intend to use this package, we recommend following the instructions above.

### Get the code

First, run the following command to clone the source code repository:

```bash
git clone git@github.com:deeporiginbio/deeporigin-client.git
```

### Install the development dependencies for this package

Second, install Python 3.9+ and
[make](https://www.gnu.org/software/make/).

To verify that both are installed, run the following commands:

```bash
python3 --version
# Python 3.12.3

make --version
# GNU Make 4.4.1
```

### Create an editable installation of this package from the cloned source code

Third, navigate to the directory you cloned the code to, and run the following command:

```bash
make install
```

This will install this package in an "editable" mode. In this mode, changes to the source code will take effect
immediately.

### Run the linting and tests for this package

Once installed, you can lint and test this package by running the following commands:

```bash
make lint
make test
```

### Compiling and serving the documentation for this package

You can compile and serve the documentation for this package by running:

```bash
make docs-serve
```

## Supported Python versions

`deeporigin` is [tested](https://github.com/deeporiginbio/deeporigin-client/actions/workflows/main.yml) against the following versions of Python in Ubuntu Linux and Windows:

- 3.10
- 3.11
- 3.12
- 3.13
