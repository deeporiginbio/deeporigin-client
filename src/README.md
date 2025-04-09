# Source Directory Contents

Here's a brief summary each sub-module functionality and some implementation details.

* `cli` - Contains the command line interface source, supporting commands such as `deeporigin config show` or `deeporigin authenticate`. There is a broad set of such commands, registered with help of the Python `cement` library Controller classes and decorators on functions. The `__init__.py` registers all of the controllers, suppports authentication and provides the `main()` entry point. Additional files provide controllers for configuration settings, data hub interactions, and variable managament for workstations.

* `config` - Defines the template for the .yml configuration file, supports reading it and setting values in it.  Configuration stores authentication, API urls and other info.

* `data` - Stores SDF and PDB files used for sample drug_discovery calls.

* `data_hub` - Privides the main public API interface to data hub, exposed as methods on the `data_hub.api` module, such as `data_hub.api.list_files(...)`. Functionality includes creating and deleting database tables and columns and rows, downloading and uploading files, writing data and accesisng it through DataFrame, etc. Most functions here provide validation and massage data before delegating to the low-level back end methods of `deeporigin_data.DeeporiginData` client, dynamically exposed through internal `data_hub._api`.

* `drug_discovery` - Defines the main interface to DeepOrign drig discovery tools, exposed though the Complex class, and its related Docking, RBFE and ABFE child clasess. Complex operates by working with a user local directory for accessing ligand and protein files, and also connecting to the platfrom and datahub. Connection is required before tool runs. Child class discovery tools expose their own pubcic API run methods such as `complex.abfe.run_end_to_end(...)`, etc. Currently, inputs are and outputs for tool runs are stored in the datahub tables, before being passed to the back end, so the apropriate tables columns are created and updated as needed. Runs deledate to `tools.run._process_job` which creates local run job records, and ultimately to `platform.tools.execute_tool` which delegates to auto-generated back end stubs from the `deeporigin-data-sdk` repo.

* `json` - Contains JSON files that provide input parameters to different back end tools, such as abfe_end_to_end.json

* `platform` - Defines a number of sub-modules that provide access to the DeepOrigin platform APIs, including clusters, users, workstations and volumes. These are wrapped through `add_functions_to_module` and dynamically pulled from the `do_sdk_platform` module, coming from the `generated-sdks` repo.

* `tools` - Defines both general tool/job run support in the `utils` module (`run_tool`, `wait_for_job`, `run_tool`, methods etc) and tool-specific functionality in the `run` module. `run` module provides methods to run some of the open-source tools we host on our platform, such as autodock vina, pdb_pdbqt_converter, and the related ligand/receptor prep. Internal logic also includes synchronizing with the data_hub as it's used for argument passing to the back end.

* `utils` - Utility methods and constants used primarily from data hub, DataFrame, and CLI. Includes methods for url construction, checking for updates, and constants such as database and urls prefixes, attribute keys and data hub supported DataType. Also `core` submodule defines PrettyDict and  a numebr of utility functions used in drug_disovery. (Would be nice to clean this up).

* `variables` - Provides a framework for managing various types of credentials, secrets, and configuration variables. It allows users retrieve, install, and uninstall different types of variables. Intended to be used on DeepOrigin workstations (?).

* `variables\types` - Defines many types of variables and secrets, implemented in individual files modules that support installing or uninstalling them. Includes support for AWS profile API keys, SSH/GPG keys, and secret files, vertain license files, etc.

