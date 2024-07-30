# deeporigin

![PyPI](https://img.shields.io/pypi/v/deeporigin)

This repository contains the `deeporigin` CLI and
Python client, which allows you to interact with
Deep Origin from the command line and Python.

> [!WARNING]  
> The `deeporigin` client is under active development. Features
> may change or be removed.

## Installing

> [!CAUTION]
> As a best practice, we recommend installing this package in a virtual environment.

To install this package, run the following:

```bash
pip install deeporigin
```

## Configuration

To run this package outside of a Deep Origin workstation (for example, on your own computer), first you need to configure this package. After installing this package, run the following to configure your organization, replacing `org-id` with the ID of the Deep Origin organization that you would like to work with.

```bash
deeporigin config set organization_id [org-id]
```

## Developing

First, download the source code from GitHub:

```bash
git clone git@github.com:deeporiginbio/deeporigin-client.git
cd deeporigin-client
```

Second, run the code below to create a virtual environment and install this package into it. This requires make v4.4 or higher.

```bash
make install
```

## Testing

### Running the tests locally

To run the tests locally, execute the following:

```bash
make test
```

By default, the tests are run using mocked responses. To run the tests against the live Deep Origin API, execute the following:

```bash
make test client=default
```

### Automated tests on GitHub Actions

The tests are automatically run on GitHub Actions on every commit to every pull request.

## License

MIT
