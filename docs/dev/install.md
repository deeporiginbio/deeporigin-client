
# Installation Guide for Developers

!!! warning "Developers only"
    If you intend to contribute this package, we recommend installing this package via the instructions below. If you simply intend to use this package, use [these instructions](../install.md).

## Get the code

First, run the following command to clone the source code repository:

```bash
git clone git@github.com:deeporiginbio/deeporigin-client.git
```

## Install the development dependencies for this package

Second, install Python 3.10+ and
[GNU make :octicons-link-external-16:](https://www.gnu.org/software/make/).

To verify that both are installed, run the following commands:

```bash
python3 --version
# Python 3.12.3

make --version
# GNU Make 4.4.1
```

## Create an editable installation of this package from the cloned source code

Third, navigate to the directory you cloned the code to, and run the following command:

```bash
make install
```

This will install this package in an "editable" mode. In this mode, changes to the source code will take effect
immediately.

## Run the linting and tests for this package

Once installed, you can test this package by running the following commands:

```bash
make test
```

## Compiling and serving the documentation for this package

You can compile and serve the documentation for this package by running:

```bash
make docs-serve
```