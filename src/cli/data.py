"""this implements controllers and hooks to connect to
managed_data.py"""

import json
from typing import Union

import cement
from beartype import beartype
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data import _api
from deeporigin.managed_data.api import download, get_children, get_row_data, upload
from deeporigin.utils import PREFIX


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

    @cement.ex(
        help="Describe and get metadata of file uploaded to database in your Deep Origin data management system",
        arguments=[
            (
                ["file_id"],
                {"help": "File ID", "action": "store"},
            )
        ],
    )
    def describe_file(self):
        """describe row"""
        _show_json(_api.describe_file(self.app.pargs.file_id))

    @cement.ex(
        help="Describe row",
        arguments=[
            (
                ["row_id"],
                {"help": "Row or Database ID", "action": "store"},
            )
        ],
    )
    def describe_row(self):
        """describe row"""
        _show_json(_api.describe_row(self.app.pargs.row_id))

    @cement.ex(
        help="List rows in a database or row",
        arguments=[
            (
                ["row_id"],
                {"help": "Row or Database ID", "action": "store"},
            )
        ],
    )
    def list_rows(self):
        """describe row"""
        _show_json(_api.list_rows(self.app.pargs.row_id))

    @cement.ex(
        help="List rows in a database",
        arguments=[
            (
                ["db_id"],
                {"help": "Database ID", "action": "store"},
            )
        ],
    )
    def list_database_rows(self):
        """describe row"""
        _show_json(_api.list_database_rows(self.app.pargs.db_id))

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
        _show_json(get_children(self.app.pargs.object_id))

    @cement.ex(
        help="Show column names and values for a given row",
        arguments=[
            (
                ["row_id"],
                {"help": "Row ID", "action": "store"},
            ),
        ],
    )
    def row(self):
        """list the columns of the row and their values, where applicable"""
        row_data = get_row_data(self.app.pargs.row_id)
        _show_json(row_data)


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
