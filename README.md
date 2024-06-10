# deeporigin 

![PyPI](https://img.shields.io/pypi/v/deeporigin)

This repository contains the `deeporigin` CLI and 
python client, which allows you to interact with 
Deep Origin from the command line. 

> [!WARNING]  
> The `deeporigin` client is under active development. Features
> may change, or be removed. 

## Installing 

> [!CAUTION]
> As a best practice, we recommend  installing this in a virtual environment. 


```bash
pip install deeporigin
```

## Developing 

First, download from Github:

```bash
git clone git@github.com:deeporiginbio/deeporigin-client.git
cd deeporigin-client
```
Using `make > v4.4`, e.g.:

```bash
make --version # 4.4.1
```

Install in a virtual env using:

```bash
make install
```

## Configuration

## Testing 

### Running tests locally 

By default, tests are run using mocked responses. To run tests locally:

```bash
make test
```

If you want to run tests against a live instance:

```bash
make test client=default
```

### Tests on GitHub Actions

Tests are run on GitHub Actions on every commit to every pull request. 


## License 

MIT