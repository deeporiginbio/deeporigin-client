"""
This module defines a class called DataFrame that is a drop-in
replacement for a pandas DataFrame, but also allows automatic
updating of Deep Origin databases.
"""

from datetime import datetime
from typing import Optional

import humanize
import pandas as pd
from deeporigin.data_hub import api
from deeporigin.utils import construct_resource_url


class DataFrame(pd.DataFrame):
    """A subclass of pandas DataFrame that allows for automatic updates to Deep Origin databases. This can be used as a drop-in replacement for a pandas DataFrame, and should support all methods a pandas DataFrame supports.

    The primary method of creating an object of this type is to use the [api.get_dataframe][src.data_hub.api.get_dataframe] function.
    """

    auto_sync: bool = False
    """When True, changes made to the dataframe will be automatically synced to the Deep Origin database this dataframe represents."""

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
                self.obj.sync(columns=columns, rows=rows)

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
            self.sync(columns=[key])

    def _repr_html_(self):
        """method override to customize printing in a Jupyter notebook"""

        name = self.attrs["metadata"]["name"]
        url = construct_resource_url(
            name=name,
            row_type="database",
        )

        # Convert the string to a datetime object
        date_str = self.attrs["metadata"]["dateCreated"]
        date_obj = datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S.%f")

        now = datetime.now()

        # Calculate the difference
        time_diff = now - date_obj

        # Convert the time difference into "x time ago" format
        time_ago = humanize.naturaltime(time_diff)

        header = f'<h4>{name} <a href = "{url}">🔗</a></h4>'
        txt = f'<p style="font-size: 12px; color: #808080;">Created {time_ago}.</p>'
        df_html = super()._repr_html_()
        return header + txt + df_html

    def __repr__(self):
        """method override to customize printing in an interactive session"""

        header = f'{self.attrs["metadata"]["hid"]}\n'
        df_representation = super().__repr__()
        return header + df_representation

    def sync(
        self,
        *,
        columns: Optional[list] = None,
        rows: Optional[list] = None,
    ):
        """Manually synchronize data in the dataframe to the underlying Deep Origin database.

        !!! tip "Deep Origin DataFrames automatically synchronize"
            Typically, you do not need to manually synchronize. If the `auto_sync` attribute of the dataframe is set to `True`, the dataframe will automatically synchronize when changes are made to the dataframe.


        Args:
            columns (list, optional): The columns of the dataframe to update. Defaults to None.
            rows (list, optional): The rows to update. Defaults to None. When None, all rows in the relevant columns are updated.

        """

        if columns is None:
            columns = self.columns

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

    @property
    def _constructor(self):
        """this method overrides the _constructor property to return a DataFrame and is required for compatibility with a pandas DataFrame"""

        return DataFrame
