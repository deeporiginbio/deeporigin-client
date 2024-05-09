"""this implements controllers and hooks to connect to
managed_data.py"""

import os

import cement
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data import _api, api
from deeporigin.managed_data.client import DeepOriginClient
from deeporigin.utils import PREFIXES, _print_dict, _print_tree, _show_json, _truncate


class DataController(cement.Controller):
    """Controller for data subcommand of CLI"""

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
            client = DeepOriginClient()  # pragma: no cover
            client.authenticate(refresh_tokens=False)  # pragma: no cover

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
                _api.download_file(file, destination=destination)

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
        help="List files, rows, databases, workspaces in Deep Origin",
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
                ["--workspaces"],
                {
                    "action": "store_true",
                    "help": "Whether to list workspaces: [False]",
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
        """list files, rows, databases, workspaces"""

        if self.app.pargs.files:
            # we will only list files
            data = _api.list_files(
                client=self._get_client(),
            )
            if not self.app.pargs.json:
                # show a table with file names, ids, status
                pdata = dict(Name=[], Status=[], ID=[])
                for item in data:
                    pdata["Name"].append(item["file"]["name"])
                    pdata["Status"].append(item["file"]["status"])
                    pdata["ID"].append(item["file"]["id"])
                _print_dict(pdata, json=False, transpose=False)
            else:
                _show_json(data)
            return

        if (
            not self.app.pargs.rows
            and not self.app.pargs.workspaces
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
                rows = _api.list_rows(client=self._get_client())
                _show_json(rows)

            return

        # at this point it is not possible to construct a tree.
        # so we will only show a table, or JSON output
        rows = []
        if self.app.pargs.rows:
            rows += _api.list_rows(
                row_type="row",
                client=self._get_client(),
            )
        if self.app.pargs.databases:
            rows += _api.list_rows(
                row_type="database",
                client=self._get_client(),
            )
        if self.app.pargs.workspaces:
            rows += _api.list_rows(
                row_type="workspace",
                client=self._get_client(),
            )

        if not self.app.pargs.json:
            pdata = dict(Name=[], Type=[], ID=[])

            for item in rows:
                pdata["Name"].append(item["name"])
                pdata["Type"].append(item["type"])
                pdata["ID"].append(item["hid"])
            _print_dict(pdata, json=False, transpose=False)
        else:
            _show_json(rows)

    @cement.ex(
        help="Describe and get metadata of file uploaded to database in your Deep Origin data management system",
        arguments=[
            (
                ["object_id"],
                {"help": "File ID or Row ID", "action": "store"},
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
        """describe file or row"""

        if PREFIXES.FILE in self.app.pargs.object_id:
            data = _api.describe_file(
                self.app.pargs.object_id,
                client=self._get_client(),
            )
        else:
            # not a file
            data = _api.describe_row(
                self.app.pargs.object_id,
                client=self._get_client(),
            )

            data.pop("rowJsonSchema", None)

            if "cols" in data.keys():
                col_names = [col["name"] for col in data["cols"]]
                col_keys = [col["key"] for col in data["cols"]]

                col_names = ", ".join(col_names)
                col_keys = ", ".join(col_keys)

                if data["type"] == "database" and not self.app.pargs.json:
                    data["Column Names"] = _truncate(col_names)
                    data["Column Keys"] = _truncate(col_keys)

                    data.pop("cols", None)

        _print_dict(data, json=self.app.pargs.json)

    @cement.ex(
        help="Show database or row",
        arguments=[
            (
                ["object_id"],
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
    def show(self):
        """show database or row in Deep Origin"""

        data = _api.describe_row(
            row_id=self.app.pargs.object_id,
            client=self._get_client(),
        )
        row_type = data["type"]

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
            _print_dict(data, json=self.app.pargs.json, transpose=True)


CONTROLLERS = [
    DataController,
]
