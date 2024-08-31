"""
This module defines a class called DataFrame that is a drop-in
replacement for a pandas DataFrame, but also allows for easy
updating of Deep Origin databases.
"""

from datetime import datetime, timezone
from typing import Optional

import humanize
import pandas as pd
from deeporigin.data_hub import api
from deeporigin.utils import (
    DatabaseReturnType,
    IDFormat,
    construct_resource_url,
)


class DataFrame(pd.DataFrame):
    """A subclass of pandas DataFrame that allows for easy updating of a Deep Origin database. This can be used as a drop-in replacement for a pandas DataFrame, and should support all methods a pandas DataFrame supports.

    The primary method of creating an object of this type is to use the [from_deeporigin][src.data_hub.dataframe.DataFrame.from_deeporigin] class method.
    """

    auto_sync: bool = False
    """When `True`, changes made to the dataframe will be automatically synced to the Deep Origin database this dataframe represents."""

    _modified_columns: set = set()
    """if data is modified in a dataframe, and auto_sync is False, this list will contain the columns that have been modified so that the Deep Origin database can be updated. If an empty list, the Deep Origin database will not be updated, and the dataframe matches the Deep Origin database at the time of creation."""

    class AtIndexer:
        """this class override is used to intercept calls to at indexer of a pandas dataframe"""

        def __init__(self, obj):
            self.obj = obj

        def __getitem__(self, key):
            """intercept for the set operation"""

            return self.obj._get_value(*key)

        def __setitem__(self, key, value):
            """intercept for the set operation""" ""

            old_value = self.obj._get_value(*key)
            if value == old_value:
                # noop
                return

            rows = [key[0]]
            columns = [key[1]]

            # Perform the actual setting operation
            self.obj._set_value(*key, value)

            # now update the DB. note that self is an AtIndexer
            # object, so we need to index into the pandas object
            if self.obj.auto_sync:
                self.obj.to_deeporigin(columns=columns, rows=rows)
            else:
                self.obj._modified_columns.add(key[1])

    @property
    def at(self):
        """Override the `at` property to return an AtIndexer"""
        return self.AtIndexer(self)

    def __setitem__(self, key, value):
        """Override the __setitem__ method to update the Deep Origin database when changes are made to the local
        dataframe"""

        # first, call the pandas method
        super().__setitem__(key, value)

        # now, update the Deep Origin database with the changes
        # we just made
        if self.auto_sync:
            self.to_deeporigin(columns=[key])
        else:
            self._modified_columns.add(key)

    def _repr_html_(self):
        """method override to customize printing in a Jupyter notebook"""

        name = self.attrs["metadata"]["name"]
        hid = self.attrs["metadata"]["hid"]
        url = construct_resource_url(
            name=hid,
            row_type="database",
        )

        # extract org name from url
        try:
            url_parts = url.split("/")
            org_name = url_parts[url_parts.index("org") + 1]
        except Exception:
            org_name = ""

        # Convert the string to a datetime object
        date_str = self.attrs["metadata"]["dateCreated"]
        date_obj = datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S.%f").replace(
            tzinfo=timezone.utc
        )

        now = datetime.now(timezone.utc)

        # Convert the time difference into "x time ago" format
        created_time_ago = humanize.naturaltime(now - date_obj)

        date_str = self.attrs["last_updated_row"].date_updated
        date_obj = datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S.%f").replace(
            tzinfo=timezone.utc
        )
        edited_time_ago = humanize.naturaltime(now - date_obj)

        header = f'<h4 style="color: #363636;">Deep Origin / {org_name} / <a href = "{url}">{name} </a></h4>'
        txt = f'<p style="font-size: 12px; color: #808080;">Created {created_time_ago}. Row {self.attrs["last_updated_row"].hid} was last edited {edited_time_ago}'
        try:
            txt += (
                "  by "
                + self.attrs["last_updated_row"].edited_by_user_drn.split("|")[1]
                + ".</p>"
            )
        except Exception:
            # give up. this should not cause the dataframe to
            # not print.
            txt += ".</p>"

        if self._modified_columns:
            txt += '<p style="color: #808080; font-size: 12px">⚠️ This dataframe contains changes that have not been written back to the Deep Origin database.</p>'
        elif self.auto_sync:
            txt += '<p style="color: #808080; font-size: 12px">🧬 This dataframe will automatically write changes made to it back to Deep Origin.</p>'
        df_html = super()._repr_html_()
        return header + txt + df_html

    def __repr__(self):
        """method override to customize printing in an interactive session"""

        header = f'{self.attrs["metadata"]["hid"]}\n'
        df_representation = super().__repr__()
        return header + df_representation

    @classmethod
    def from_deeporigin(
        cls,
        database_id: str,
        *,
        use_file_names: bool = True,
        reference_format: IDFormat = "human-id",
        return_type: DatabaseReturnType = "dataframe",
        client=None,
    ):
        """Create a local Deep Origin DataFrame from a Deep Origin database.

        Args:
            database_id (str): The ID of the Deep Origin database.
            use_file_names (bool, optional): Whether to use the file names in the Deep Origin database. Defaults to True.
            reference_format (IDFormat, optional): The format of the IDs in the Deep Origin database. Defaults to "human-id".
            return_type (DatabaseReturnType, optional): The type of return value. Defaults to "dataframe".

        """

        return api.get_dataframe(
            database_id=database_id,
            use_file_names=use_file_names,
            reference_format=reference_format,
            return_type=return_type,
            client=client,
        )

    def to_deeporigin(
        self,
        *,
        columns: Optional[list] = None,
        rows: Optional[list] = None,
    ):
        """Write data in dataframe to Deep Origin

        !!! tip "Deep Origin DataFrames can automatically synchronize"
            To automatically save changes to local DataFrames to Deep Origin databases, set the `auto_sync` attribute of the dataframe `True`.


        Args:
            columns (list, optional): The columns of the dataframe to update. When None, all modified columns are updated.
            rows (list, optional): The rows to update. Defaults to None. When None, all rows in the relevant columns are updated.

        """

        if columns is None:
            columns = self._modified_columns.copy()

        for column in columns:
            if column in ["Validation Status", "ID"]:
                continue

            column_metadata = [
                col for col in self.attrs["metadata"]["cols"] if col["name"] == column
            ]

            if len(column_metadata) == 0:
                raise NotImplementedError(
                    "Column metadata not found. This is likely because it's a new column"
                )

            column_metadata = column_metadata[0]

            if column_metadata["type"] == "file":
                continue

            if rows is None:
                # we're updating a whole column
                rows = list(self.index)

            api.set_data_in_cells(
                values=self[column],
                row_ids=rows,
                column_id=column,
                database_id=self.attrs["id"],
            )

            print(f"✔︎ Wrote {len(rows)} rows in {column} to Deep Origin database.")
            self._modified_columns.discard(column)

    @property
    def _constructor(self):
        """this method overrides the _constructor property to return a DataFrame and is required for compatibility with a pandas DataFrame"""

        return DataFrame
