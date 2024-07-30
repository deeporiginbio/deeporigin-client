"""this implements controllers and hooks to connect to
data_hub.py"""

import os

import cement
from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils import PREFIXES, _print_dict, _print_tree, _show_json, _truncate


class DataController(cement.Controller):
    """Controller for data subcommand of CLI"""

    class Meta:
        label = "data"
        stacked_on = "base"
        stacked_type = "nested"
        help = "Explore and fetch data from the Deep Origin data hub"
        description = """List data in the data hub on Deep Origin, and save databases to CSV files."""

    def _get_client(self):
        """helper method that returns an authenticated
        client if the app has no client configured"""
        try:
            return self.app.client
        except Exception:
            client = api._get_default_client()

            return client  # pragma: no cover

    @cement.ex(
        help="Merge two databases into a single one, integrating cross-references",
        arguments=[
            (
                ["--databases"],
                {
                    "help": "List of databases to merge",
                    "action": "store",
                    "nargs": "*",
                    "required": True,
                },
            ),
            (
                ["--destination"],
                {
                    "type": str,
                    "required": True,
                    "metavar": "<destination>",
                    "help": "Folder on local disk to save to",
                },
            ),
            (
                ["--include-files"],
                {
                    "action": "store_true",
                    "help": "Whether to download files in database [default: False]",
                },
            ),
        ],
    )
    def merge_db(self):
        """Merge databases and save as CSV"""

        databases = self.app.pargs.databases
        destination = self.app.pargs.destination
        dfs = []
        save_name = "-".join(databases)
        for db in databases:
            df = api.get_dataframe(
                db,
                client=self._get_client(),
            )
            dfs.append(df)
        df = api.merge_databases(dfs)

        # save to CSV
        df.to_csv(os.path.join(destination, f"merged-{save_name}.csv"))

        if not self.app.pargs.include_files:
            return

        for df in dfs:
            files = df.attrs["file_ids"]
            for file in files:
                api.download_file(file, destination=destination)

    @cement.ex(
        help="Copy files or databases from or to Deep Origin",
        arguments=[
            (
                ["source"],
                {"help": "Source", "action": "store"},
            ),
            (
                ["destination"],
                {"help": "Destination", "action": "store"},
            ),
            (
                ["--include-files"],
                {
                    "action": "store_true",
                    "help": "Whether to also download files in a database: [False]",
                },
            ),
        ],
    )
    def copy(self):
        """download or upload files or databases"""

        args = self.app.pargs

        if PREFIXES.DO in args.source and PREFIXES.DO not in args.destination:
            api.download(
                args.source,
                args.destination,
                include_files=args.include_files,
            )
        elif PREFIXES.DO in args.destination and PREFIXES.DO not in args.source:
            raise NotImplementedError("Uploading has not been implemented yet")
            # upload(args.source, args.destination)
        else:
            raise DeepOriginException(
                f"Exactly one of <source> and <destination> should be prefixed with `{PREFIXES.DO}`"
            )

    @cement.ex(
        help="List files, rows, databases, folders in Deep Origin",
        arguments=[
            (
                ["--files"],
                {
                    "action": "store_true",
                    "help": "Whether to list files: [False]",
                },
            ),
            (
                ["--rows"],
                {
                    "action": "store_true",
                    "help": "Whether to list rows: [False]",
                },
            ),
            (
                ["--folders"],
                {
                    "action": "store_true",
                    "help": "Whether to list folders: [False]",
                },
            ),
            (
                ["--databases"],
                {
                    "action": "store_true",
                    "help": "Whether to list databases: [False]",
                },
            ),
            (
                ["--json"],
                {
                    "action": "store_true",
                    "help": "Whether to return JSON formatted data [default: [False]",
                },
            ),
        ],
    )
    def list(self):
        """list files, rows, databases, folders"""

        if self.app.pargs.files:
            # we will only list files
            files = api.list_files(
                client=self._get_client(),
            )
            if not self.app.pargs.json:
                # show a table with file names, ids, status
                pdata = dict(Name=[], Status=[], ID=[])
                for item in files:
                    pdata["Name"].append(item.file.name)
                    pdata["Status"].append(item.file.status)
                    pdata["ID"].append(item.file.id)
                _print_dict(pdata, json=False, transpose=False)
            else:
                files = [file.dict() for file in files]
                _show_json(files)
            return

        if (
            not self.app.pargs.rows
            and not self.app.pargs.folders
            and not self.app.pargs.databases
        ):
            # we want to list everything, show show a tree.
            if not self.app.pargs.json:
                tree = api.get_tree(
                    client=self._get_client(),
                )
                for branch in tree:
                    _print_tree(branch)
            else:
                rows = api.list_rows(client=self._get_client())
                rows = [row.dict() for row in rows]
                _show_json(rows)

            return

        # at this point it is not possible to construct a tree.
        # so we will only show a table, or JSON output
        rows = []
        if self.app.pargs.rows:
            rows += api.list_rows(
                row_type="row",
                client=self._get_client(),
            )
        if self.app.pargs.databases:
            rows += api.list_rows(
                row_type="database",
                client=self._get_client(),
            )
        if self.app.pargs.folders:
            rows += api.list_rows(
                row_type="workspace",
                client=self._get_client(),
            )

        if not self.app.pargs.json:
            pdata = dict(Name=[], Type=[], ID=[])

            for item in rows:
                pdata["Name"].append(item.name)
                pdata["Type"].append(item.type)
                pdata["ID"].append(item.hid)
            _print_dict(pdata, json=False, transpose=False)
        else:
            rows = [row.dict() for row in rows]
            _show_json(rows)

    @cement.ex(
        help="Describe and get metadata of about a file, row, or database in your Deep Origin data hub",
        arguments=[
            (
                ["object_id"],
                {"help": "File ID or row ID", "action": "store"},
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
    def describe(self):
        """describe file or row or database"""

        if PREFIXES.FILE in self.app.pargs.object_id:
            key_label = "Property"

            data = api.describe_file(
                file_id=self.app.pargs.object_id,
                client=self._get_client(),
            )
            data = data.dict()
        else:
            # not a file
            key_label = "Column"

            data = api.describe_row(
                row_id=self.app.pargs.object_id,
                client=self._get_client(),
            )

            data = dict(data)

            data.pop("row_json_schema", None)
            data.pop("rowJsonSchema", None)
            data.pop("editor", None)

            if "cols" in data.keys() and data["cols"] is not None:
                col_names = [col["name"] for col in data["cols"]]
                col_keys = [col["key"] for col in data["cols"]]

                col_names = ", ".join(col_names)
                col_keys = ", ".join(col_keys)

                if data["type"] == "database" and not self.app.pargs.json:
                    data["Column Names"] = _truncate(col_names)
                    data["Column Keys"] = _truncate(col_keys)

                    data.pop("cols", None)

        _print_dict(data, json=self.app.pargs.json, key_label=key_label)

    @cement.ex(
        help="Show a row or a database",
        arguments=[
            (
                ["object_id"],
                {"help": "Row ID or database ID", "action": "store"},
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
    def show(self):
        """show database or row in Deep Origin"""

        data = api.describe_row(
            row_id=self.app.pargs.object_id,
            client=self._get_client(),
        )
        row_type = data.type

        if row_type == "database":
            data = api.get_dataframe(
                self.app.pargs.object_id,
                return_type="dict",
                client=self._get_client(),
            )
            _print_dict(data, json=self.app.pargs.json, transpose=False)
        elif row_type == "row":
            data = api.get_row_data(
                self.app.pargs.object_id,
                client=self._get_client(),
            )
            _print_dict(
                data, json=self.app.pargs.json, transpose=True, key_label="Column"
            )

    @cement.ex(
        help="Upload a file to database",
        arguments=[
            (
                ["file_path"],
                {"help": "File path to upload", "action": "store"},
            ),
            (
                ["--database"],
                {
                    "type": str,
                    "required": False,
                    "metavar": "<database_id>",
                    "help": "ID of database to assign to",
                },
            ),
            (
                ["--column"],
                {
                    "type": str,
                    "required": False,
                    "metavar": "<column_id>",
                    "help": "ID of column to assign to",
                },
            ),
            (
                ["--row"],
                {
                    "type": str,
                    "required": False,
                    "metavar": "<row_id>",
                    "help": "ID of row to assign to",
                },
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
    def upload(self):
        """upload file to database in Deep Origin"""

        data = api.upload_file(
            file_path=self.app.pargs.file_path,
            client=self._get_client(),
        )

        if not self.app.pargs.database:
            # we are not making an assignment, so abort
            _print_dict(data.dict(), json=self.app.pargs.json, key_label="Property")
            return

        if self.app.pargs.column and self.app.pargs.database:
            data = api.assign_files_to_cell(
                file_ids=[data.id],
                database_id=self.app.pargs.database,
                column_id=self.app.pargs.column,
                row_id=self.app.pargs.row,
            )

            data = data.rows

            data = [row.dict() for row in data][0]
            data.pop("fields", None)

            _print_dict(
                data,
                json=self.app.pargs.json,
                transpose=True,
                key_label="Property",
            )

    @cement.ex(
        help="Write data to database cell",
        arguments=[
            (
                ["data"],
                {"help": "Data to set", "action": "store"},
            ),
            (
                ["--database"],
                {
                    "type": str,
                    "required": True,
                    "metavar": "<database_id>",
                    "help": "ID of database to write to",
                },
            ),
            (
                ["--column"],
                {
                    "type": str,
                    "required": True,
                    "metavar": "<column_id>",
                    "help": "ID of column to write to",
                },
            ),
            (
                ["--row"],
                {
                    "type": str,
                    "required": True,
                    "metavar": "<row_id>",
                    "help": "ID of row to write to",
                },
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
    def write(self):
        """Write data to a cell in a database"""

        api.set_cell_data(
            self.app.pargs.data,
            database_id=self.app.pargs.database,
            column_id=self.app.pargs.column,
            row_id=self.app.pargs.row,
        )

        print(f"✔︎ Wrote {self.app.pargs.data} to database")


CONTROLLERS = [
    DataController,
]
