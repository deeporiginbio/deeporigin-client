# Installation


## On a Deep Origin workstation

No installation needed!

The Deep Origin CLI and Python client are installed on every [Deep Origin Workstation :octicons-link-external-16:](https://docs.deeporigin.io/docs/os/compute-hub/workstations).

## Recommended installation

This sections describes how to install the Deep Origin CLI and Python client on your computer using our recommendations. We recommend:

- Using the [uv :octicons-link-external-16:](https://docs.astral.sh/uv/) package and environment manager 
- Using [Jupyter Lab :octicons-link-external-16:](https://jupyter.org/) to run the Deep Origin Python client

=== "macOS/Linux"

    First, create a folder to contain your project. This can be any name (except deeporigin).

    ```bash
    mkdir do-client
    cd do-client
    ```

    Now use this one liner to install deeporigin using the [uv :octicons-link-external-16:](https://docs.astral.sh/uv/) package manager. 

    ```bash
    curl -LsSf https://client-docs.deeporigin.io/install.sh | sh
    ```

    Now start a Jupyter Lab instance that uses this environment, allowing you to use the `deeporigin` python client in it:


    ```bash
    # run jupyter using this environment
    uv run --with jupyter jupyter lab

    ```


=== "Windows"

    ```bash
    # make a folder for your project
    # this can be any name (except deeporigin)
    mkdir do-client
    cd do-client

    # install uv if needed
    powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

    # If uv tells you to add something to PATH,
    # you may need to do that first
    uv init

    # we use the latest stable version of python
    uv python install 3.13

    # install deeporigin in the current environment
    uv add deeporigin

    # if you're using Deep Origin tools (like FEP), use:
    uv add deeporigin --extra tools

    # also install as a uv tool to run from the command line
    uv tool install deeporigin

    # run jupyter using this environment
    uv run --with jupyter jupyter lab
    ```
## Upgrading to a new version

If you followed the recommended installation steps, you can upgrade to the latest version of the Deep Origin Python client by running:

```bash
uv add --upgrade deeporigin
```

in your project root. 

## For developers

!!! warning "Developers only"
    If you intend to contribute this package, we recommend installing this package via the instructions below. If you simply intend to use this package, use the instructions above.

### Get the code

First, run the following command to clone the source code repository:

```bash
git clone git@github.com:deeporiginbio/deeporigin-client.git
```

### Install the development dependencies for this package

Second, install Python 3.10+ and
[GNU make :octicons-link-external-16:](https://www.gnu.org/software/make/).

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

Once installed, you can test this package by running the following commands:

```bash
make test
```

### Compiling and serving the documentation for this package

You can compile and serve the documentation for this package by running:

```bash
make docs-serve
```

## OS and python support

`deeporigin` is [tested](https://github.com/deeporiginbio/deeporigin-client/actions/workflows/main.yml) against the following versions of Python on Ubuntu Linux and Windows:

| Python | macOS | Windows | Ubuntu |
| -- | -- | -- | -- | 
| 3.10| | ✅ | ✅ |
| 3.11| | ✅ | ✅ |
| 3.12| | ✅ | ✅ |
| 3.13 | ✅ | ✅ | ✅ |


