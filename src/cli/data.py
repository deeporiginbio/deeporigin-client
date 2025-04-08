"""this implements controllers and hooks to connect to
data_hub.py"""

import cement

from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils.constants import (
    PREFIXES,
    DataType,
)
from deeporigin.utils.core import (
    _print_dict,
    _print_tree,
    _show_json,
    _truncate,
    humanize_file_size,
)


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
        help="Download files from Deep Origin to local computer",
        arguments=[
            (
                ["--file-ids"],
                {
                    "help": "IDs of files to download",
                    "nargs": "+",
                },
            ),
            (
                ["--assigned-row-ids"],
                {
                    "help": "IDs of rows that files are assigned to to download.",
                    "nargs": "+",
                },
            ),
        ],
    )
    def download_files(self):
        """Download multiple files from Deep Origin"""

        print(self.app.pargs.assigned_row_ids)

        files = api.list_files(
            file_ids=self.app.pargs.file_ids,
            assigned_row_ids=self.app.pargs.assigned_row_ids,
        )

        total_size = humanize_file_size(sum(file.file.content_length for file in files))

        print(f"{len(files)} files ({total_size}) will be downloaded")

        api.download_files(files)

    @cement.ex(
        help="Copy a file or database to or from your data hub",
        arguments=[
            (
                ["source"],
                {
                    "help": "Source: ID of a file or database, or local path for a file or database",
                    "action": "store",
                },
            ),
            (
                ["destination"],
                {
                    "help": "Destination: ID of a file or database, or local path for a file or database",
                    "action": "store",
                },
            ),
            (
                ["--include-files"],
                {
                    "action": "store_true",
                    "help": "Whether to also download the files in the database: [False]",
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
                f"Exactly one of <source> or <destination> should be prefixed with `{PREFIXES.DO}`"
            )

    @cement.ex(
        help="List the files, rows, databases, and/or folders (workspaces) in your Deep Origin data hub",
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
                    "help": "Whether to return data in JSON format [default: [False]",
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
                pdata["Name"].append(getattr(item, "name", None))

                pdata["Type"].append(item.type)
                pdata["ID"].append(item.hid)
            _print_dict(pdata, json=False, transpose=False)
        else:
            _show_json(rows)

    @cement.ex(
        help="Show the metadata about a file, row, or database in your data hub",
        arguments=[
            (
                ["object_id"],
                {"help": "ID for the file, row, or database", "action": "store"},
            ),
            (
                ["--json"],
                {
                    "action": "store_true",
                    "help": "Whether to return data in JSON format [default: False]",
                },
            ),
        ],
    )
    def describe(self):
        """describe file or row or database"""

        key_label = "Property"

        col_data = None

        if PREFIXES.FILE in self.app.pargs.object_id:
            data = api.describe_file(
                file_id=self.app.pargs.object_id,
                client=self._get_client(),
            )

        else:
            # not a file

            data = api.describe_row(
                row_id=self.app.pargs.object_id,
                client=self._get_client(),
                fields=False,
            )

            data.pop("row_json_schema", None)
            data.pop("rowJsonSchema", None)
            data.pop("editor", None)

            if "cols" in data.keys() and data["cols"] is not None:
                col_names = [col["name"] for col in data["cols"]]
                col_types = [col["type"] for col in data["cols"]]
                col_ids = [col["id"] for col in data["cols"]]

                col_names_str = ", ".join(col_names)

                col_data = dict(
                    Name=col_names,
                    Type=col_types,
                    ID=col_ids,
                )

                if data["type"] == "database" and not self.app.pargs.json:
                    data["Column Names"] = _truncate(col_names)
                    data["Column Keys"] = _truncate(col_names_str)

                    data.pop("cols", None)

        _print_dict(data, json=self.app.pargs.json, key_label=key_label)

        if not self.app.pargs.json and col_data is not None:
            print("Column information:")
            _print_dict(
                col_data,
                json=False,
                transpose=False,
            )

    @cement.ex(
        help="Show a row or a database from your data hub",
        arguments=[
            (
                ["object_id"],
                {"help": "ID for the row or database", "action": "store"},
            ),
            (
                ["--json"],
                {
                    "action": "store_true",
                    "help": "Whether to return data in JSON format [default: False]",
                },
            ),
            (
                ["--notebook"],
                {
                    "action": "store_true",
                    "help": "Whether to show the notebook entry for each row [default: False]",
                },
            ),
        ],
    )
    def show(self):
        """show database or row in Deep Origin"""

        if self.app.pargs.notebook:
            # show notebook
            document = api.get_body_document(
                row_id=self.app.pargs.object_id,
                client=self._get_client(),
            )

            print(document)
            return

        data = api.describe_row(
            row_id=self.app.pargs.object_id,
            client=self._get_client(),
        )
        hid = data.hid
        row_type = data.type

        if row_type == "database":
            data = api.get_dataframe(
                self.app.pargs.object_id,
                return_type="dict",
                reference_format="system-id",
                client=self._get_client(),
            )

            _print_dict(data, json=self.app.pargs.json, transpose=False)
        elif row_type == "row":
            data = api.get_row_data(
                self.app.pargs.object_id,
                client=self._get_client(),
            )

            # insert HID as the first column
            data["ID"] = hid

            _print_dict(
                data, json=self.app.pargs.json, transpose=True, key_label="Column"
            )

    @cement.ex(
        help="Upload a file to a database, column of a database, or a cell of a database in your data hub",
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
                    "help": "ID of the database to assign the file to",
                },
            ),
            (
                ["--column"],
                {
                    "type": str,
                    "required": False,
                    "metavar": "<column_id>",
                    "help": "ID of the column to assign the file to",
                },
            ),
            (
                ["--row"],
                {
                    "type": str,
                    "required": False,
                    "metavar": "<row_id>",
                    "help": "ID of the row to assign the file to",
                },
            ),
            (
                ["--json"],
                {
                    "action": "store_true",
                    "help": "Whether to return data in JSON format [default: False]",
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
            _print_dict(data, json=self.app.pargs.json, key_label="Property")
            return

        if self.app.pargs.column and self.app.pargs.database:
            data = api.assign_files_to_cell(
                file_ids=[data.id],
                database_id=self.app.pargs.database,
                column_id=self.app.pargs.column,
                row_id=self.app.pargs.row,
            )

            data = data.rows[0]

            data.pop("fields", None)

            _print_dict(
                data,
                json=self.app.pargs.json,
                transpose=True,
                key_label="Property",
            )

    @cement.ex(
        help="Write data to a cell of a database in your data hub",
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
                    "help": "ID of the database to write to",
                },
            ),
            (
                ["--column"],
                {
                    "type": str,
                    "required": True,
                    "metavar": "<column_id>",
                    "help": "ID of the column to write to",
                },
            ),
            (
                ["--row"],
                {
                    "type": str,
                    "required": True,
                    "metavar": "<row_id>",
                    "help": "ID of the row to write to",
                },
            ),
            (
                ["--json"],
                {
                    "action": "store_true",
                    "help": "Whether to return data in JSON format [default: False]",
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

    @cement.ex(
        help="Create a new folder (workspace), database, or database column",
        arguments=[
            (
                ["object_type"],
                {"help": "Type of resource to create: folder, database, or column "},
            ),
            (
                ["--name"],
                {
                    "type": str,
                    "required": True,
                    "help": "Name of database, folder, or column to create",
                },
            ),
            (
                ["--parent-id"],
                {
                    "type": str,
                    "required": False,
                    "help": "ID of the parent folder to create the folder or database in",
                },
            ),
            (
                ["--database"],
                {
                    "type": str,
                    "required": False,
                    "help": "ID of database to create the column in",
                },
            ),
            (
                ["--key"],
                {
                    "type": str,
                    "required": False,
                    "help": "Programmatic key of the column to create",
                },
            ),
            (
                ["--type"],
                {
                    "type": str,
                    "required": False,
                    "help": "Type of the column to create",
                },
            ),
            (
                ["--json"],
                {
                    "action": "store_true",
                    "help": "Whether to return data in JSON format [default: False]",
                },
            ),
        ],
    )
    def new(self):
        """Create a new database, column, or row in your data hub"""

        if self.app.pargs.object_type not in ["database", "folder", "column"]:
            raise DeepOriginException(
                "First argument should be one of [database, folder, column]"
            )

        if self.app.pargs.object_type == "database":
            api.create_database(
                name=self.app.pargs.name,
                client=self._get_client(),
                parent_id=self.app.pargs.parent_id,
            )

        elif self.app.pargs.object_type == "folder":
            api.create_workspace(
                name=self.app.pargs.name,
                client=self._get_client(),
                parent_id=self.app.pargs.parent_id,
            )
        elif self.app.pargs.object_type == "column":
            if self.app.pargs.database is None:
                raise DeepOriginException(
                    "You must specify a database to create a column in using --database"
                )

            if self.app.pargs.type is None:
                raise DeepOriginException(
                    f"You must specify a type for a column from one of {DataType} using --type"
                )
            api.add_database_column(
                database_id=self.app.pargs.database,
                type=self.app.pargs.type,
                name=self.app.pargs.name,
            )

        print(
            f"✔︎ Created a new {self.app.pargs.object_type} with name: {self.app.pargs.name}"
        )

    @cement.ex(
        help="Delete rows, columns, databases and/or folders (workspaces)",
        arguments=[
            (
                ["--database", "-d"],
                {
                    "type": str,
                    "required": False,
                    "help": "ID of the database to delete",
                },
            ),
            (
                ["--folder", "--workspace", "-w", "--ws", "-f"],
                {
                    "type": str,
                    "required": False,
                    "help": "ID of folder to delete",
                },
            ),
            (
                ["--column", "-c"],
                {
                    "type": str,
                    "required": False,
                    "help": "Column ID to delete",
                },
            ),
            (
                ["--row", "-r"],
                {
                    "type": str,
                    "required": False,
                    "help": "Row ID to delete",
                },
            ),
        ],
    )
    def delete(self):
        """Delete rows, columns, databases and/or folders"""

        if self.app.pargs.column:
            if self.app.pargs.database is None:
                raise DeepOriginException(
                    "Use the --database argument to specify the parent database. To delete a column, the parent database must be specified."
                )
            api.delete_database_column(
                column_id=self.app.pargs.column,
                database_id=self.app.pargs.database,
                client=self._get_client(),
            )
            print(
                f"✔︎ Deleted column: {self.app.pargs.column} in database: {self.app.pargs.database}"
            )
        elif self.app.pargs.row:
            if self.app.pargs.database is None:
                raise DeepOriginException(
                    "Use the --database argument to specify the parent database. To delete a row, the parent database must be specified."
                )
            api.delete_rows(
                row_ids=[self.app.pargs.row],
                database_id=self.app.pargs.database,
                client=self._get_client(),
            )
            print(
                f"✔︎ Deleted row: {self.app.pargs.row} in database: {self.app.pargs.database}"
            )
        elif self.app.pargs.database:
            api.delete_database(
                database_id=self.app.pargs.database,
                client=self._get_client(),
            )
            print(f"✔︎ Deleted database: {self.app.pargs.database}")
        elif self.app.pargs.folder:
            api.delete_workspace(
                workspace_id=self.app.pargs.folder,
                client=self._get_client(),
            )
            print(f"✔︎ Deleted folder: {self.app.pargs.folder}")


CONTROLLERS = [
    DataController,
]
