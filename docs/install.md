
# Installation

## For Users

### Install latest version of code 

If you plan to use this in your own project, use your 
favorite package manager to install this in your project.



=== "pip"

    ```bash
    # you are strongly encouraged to install in a virtual envrionment
    pip install deeporigin
    ```

=== "poetry"


    ```bash
    poetry add deeporigin
    ```



## For Developers

!!! danger "Developers only"
    The instructions below are only if you intend to develop `deeporigin`. If you merely want to use this API to read data, you do not have to do this.

### Get the code

```bash
git clone git@github.com:formiclabs/deeporigin-client.git
```

### Prerequisites 

Make sure you have Python 3.9+ and
[make](https://www.gnu.org/software/make//) installed. 
Verify that both are installed:

```bash
python3 --version
# Python 3.12.3

make --version
# GNU Make 4.4.1
```

### Install locally

You can then install the API locally by navigating to the directory you downloaded the code in, and running:

```bash
make install
```

This installs in "editable" mode, so changes you made take effect
immediately. 


### Running tests

Once installed, tests can be run using:

```bash
make test
```