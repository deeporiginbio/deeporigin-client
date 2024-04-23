"""this implements controllers and hooks to connect to
managed_data.py"""

import json
from typing import Union

import cement
from beartype import beartype
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data import _api
from deeporigin.managed_data.api import (
    download,
    get_dataframe,
    get_row_data,
    get_tree,
    upload,
)
from deeporigin.managed_data.client import DeepOriginClient
from deeporigin.utils import PREFIX
from tabulate import tabulate


@beartype
def _print_tree(tree: dict, offset: int = 0) -> None:
    """helper function to pretty print a tree"""
    print(" " * offset + tree["hid"])

    if "children" not in tree.keys():
        return
    for child in tree["children"]:
        _print_tree(child, offset + 2)


def _print_dict(
    data: dict,
    *,
    json: bool = True,
    transpose: bool = True,
) -> None:
    """helper function to pretty print a dict as a table"""

    if json:
        _show_json(data)
    else:
        if transpose:
            data = data.items()
            headers = ["Name", "Value"]
        else:
            headers = "keys"
        print(
            tabulate(
                data,
                headers=headers,
                tablefmt="rounded_outline",
            )
        )


class DataController(cement.Controller):
    class Meta:
        label = "data"
        stacked_on = "base"
        stacked_type = "nested"
        help = "explore and fetch data from Deep Origin managed data"
        description = """
List data in managed data on Deep Origin, and save
databases to CSV files. 
            """

    def _get_client(self):
        """helper method that returns an authenticated
        client if the app has no client configured"""
        try:
            return self.app.client
        except Exception:
            client = DeepOriginClient()
            client.authenticate(refresh_tokens=False)

            return client

    @cement.ex(
        help="Describe and get metadata of file uploaded to database in your Deep Origin data management system",
        arguments=[
            (
                ["file_id"],
                {"help": "File ID", "action": "store"},
            ),
            (
                ["--json"],
                {
                    "action": "store_true",
                    "help": "Whether to return JSON formatted data [default: False]",
                },
            ),
        ],
    )
    def describe_file(self):
        """describe file"""

        data = _api.describe_file(
            self.app.pargs.file_id,
            client=self._get_client(),
        )

        _print_dict(data, json=self.app.pargs.json)

    @cement.ex(
        help="Describe row",
        arguments=[
            (
                ["row_id"],
                {"help": "Row or Database ID", "action": "store"},
            ),
            (
                ["--json"],
                {
                    "action": "store_true",
                    "help": "Whether to return JSON formatted data [default: False]",
                },
            ),
        ],
    )
    def describe_row(self):
        """describe row"""

        data = _api.describe_row(
            self.app.pargs.row_id,
            client=self._get_client(),
        )

        _print_dict(data, json=self.app.pargs.json)

    @cement.ex(
        help="List rows in a database or row",
        arguments=[
            (
                ["row_id"],
                {"help": "Row or Database ID", "action": "store"},
            ),
            (
                ["--json"],
                {
                    "action": "store_true",
                    "help": "Whether to return JSON formatted data [default: False]",
                },
            ),
        ],
    )
    def list_rows(self):
        """list rows"""

        rows = _api.list_rows(
            parent_id=self.app.pargs.row_id,
            client=self._get_client(),
        )

        if self.app.pargs.json:
            _show_json(rows)
            return

        # convert a list of dicts to a dict of lists
        data = {}
        keys = rows[0].keys()
        for key in keys:
            if key == "parentId" or key == "type":
                continue
            data[key] = [row[key] for row in rows]

        _print_dict({"Parent ID": rows[0]["parentId"]}, json=False)
        _print_dict(data, json=False, transpose=False)

    @cement.ex(
        help="List rows in a database",
        arguments=[
            (
                ["db_id"],
                {"help": "Database ID", "action": "store"},
            ),
            (
                ["--json"],
                {
                    "action": "store_true",
                    "help": "Whether to return JSON formatted data [default: False]",
                },
            ),
        ],
    )
    def show_db(self):
        """list database row"""

        data = get_dataframe(self.app.pargs.db_id, return_type="dict")

        _print_dict(data, json=self.app.pargs.json, transpose=False)

    @cement.ex(
        help="List child workspaces and databases in given workspace",
        arguments=[
            (
                ["object_id"],
                {
                    "help": "Workspace ID",
                    "action": "store",
                    "nargs": "?",
                    "const": None,
                },
            )
        ],
    )
    def ls(self):
        """list rows in db"""

        tree = get_tree(
            client=self._get_client(),
        )
        _print_tree(tree)

    @cement.ex(
        help="Show column names and values for a given row",
        arguments=[
            (
                ["row_id"],
                {"help": "Row ID", "action": "store"},
            ),
            (
                ["--json"],
                {
                    "action": "store_true",
                    "help": "Whether to return JSON formatted data [default: False]",
                },
            ),
        ],
    )
    def row(self):
        """list the columns of the row and their values, where applicable"""
        row_data = get_row_data(
            self.app.pargs.row_id,
            client=self._get_client(),
        )

        _print_dict(row_data, json=self.app.pargs.json)


class CopyController(cement.Controller):
    class Meta:
        label = "cp"
        stacked_on = "data"
        stacked_type = "nested"
        help = "Copy cells, columns, rows, files or databases"
        description = """Copy data from cells, columns, rows,
files or databases to or from a local filesystem. Databases
are downloaded as CSV files; included files are downloaded 
as-is. """
        arguments = [
            (
                ["source"],
                {
                    "type": str,
                    "metavar": "<source>",
                    "help": "Identifier of Remote Resource: (Workspace, Database, or Row), or path to local folder or file",
                },
            ),
            (
                ["destination"],
                {
                    "type": str,
                    "metavar": "<destination>",
                    "help": "Identifier of Remote Resource: (Workspace, Database, or Row), or path to local folder or file",
                },
            ),
            (
                ["--include-files"],
                {
                    "action": "store_true",
                    "help": "Whether to download files in database [default: False]",
                },
            ),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs

        if PREFIX in args.source and PREFIX not in args.destination:
            download(
                args.source,
                args.destination,
                include_files=args.include_files,
            )
        elif PREFIX in args.destination and PREFIX not in args.source:
            upload(args.source, args.destination)
        else:
            raise DeepOriginException(
                f"Exactly one of <source> and <destination> should be prefixed with `{PREFIX}`"
            )


@beartype
def _show_json(data: Union[list, dict]) -> None:
    """utility for pretty printing JSON"""
    print(json.dumps(data, indent=2))


CONTROLLERS = [
    DataController,
    CopyController,
]
