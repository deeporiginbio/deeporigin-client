# deeporigin 

This repository contains the `deeporigin` CLI and 
python client, which allows you to interact with 
Deep Origin from the command line. 

## Installation 

> [!CAUTION]
> You are strongly advised to install this in a virtual environment. 

### For end users

```bash
pip install deeporigin
```

### For developers

First, download from Github:

```bash
git clone git@github.com:formiclabs/cli.git
cd cli
```
Set up your virtual environment using your preferred method. 
Then,

```bash
pip install -e .[test,jupyter]
```

## Configuration

## Testing 

### Running tests locally using mocked responses

Tests can be run against a fully mocked service. This allows
you to test the CLI without requiring a connection to Deep Origin,
and tests run much faster because they're not waiting for a 
network response. 

To run using mocking, set the following in your `~/.deeporigin/config.yml` file:

```bash
nucleus_api_endpoint: https://deeporigin.mock/
```

and then run tests using:

```bash
make test
```

### Running tests locally against a live instance

Tests can be run against a live instance of Deep Origin. Point
the CLI to a URL that resolves to a live instance of Deep Origin
using the `nucleus_api_endpoint` field in your `~/.deeporigin/config.yml` file.

Then, run tests using:

```bash
make test
```

### Tests on Github Actions

## License 

MIT