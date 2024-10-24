"""
This module defines a class called DataFrame that is a drop-in
replacement for a pandas DataFrame, but also allows for easy
updating of Deep Origin databases.
"""

from datetime import datetime, timezone

import humanize
import pandas as pd
from beartype import beartype
from dateutil.parser import parse
from deeporigin.data_hub import api
from deeporigin.platform.api import get_last_edited_user_name
from deeporigin.utils.config import construct_resource_url
from deeporigin.utils.constants import IDFormat
from deeporigin.utils.network import check_for_updates

check_for_updates()


__NO_NEW_ROWS_MSG__ = "Adding rows is not allowed, because this dataframe corresponds to a subset of the rows in the corresponding database."


class DataFrame(pd.DataFrame):
    """A subclass of pandas DataFrame that allows for easy updating of a Deep Origin database. This can be used as a drop-in replacement for a pandas DataFrame, and should support all methods a pandas DataFrame supports.

    The primary method of creating an object of this type is to use the [from_deeporigin][src.data_hub.dataframe.DataFrame.from_deeporigin] class method.
    """

    auto_sync: bool = False
    """When `True`, changes made to the dataframe will be automatically synced to the Deep Origin database this dataframe represents."""

    _modified_columns: dict = dict()
    """if data is modified in a dataframe, and auto_sync is False, this list will contain the columns that have been modified so that the Deep Origin database can be updated. If an empty list, the Deep Origin database will not be updated, and the dataframe matches the Deep Origin database at the time of creation."""

    _allow_adding_rows: bool = True
    """If `True`, new rows can be added to the dataframe. If `False`, new rows cannot be added to the dataframe."""

    @property
    def loc(self):
        class _LocIndexer:
            def __init__(self, df):
                self.df = df

            def __getitem__(self, key):
                # first call the superclass method
                df = super(DataFrame, self.df).loc[key]

                # inherit attributes
                df.attrs = self.df.attrs

                # disallow adding rows
                df._allow_adding_rows = False

                df._modified_columns = self.df._modified_columns

                return df

            def __setitem__(self, key, value):
                """callback for adding a new row, typically"""

                if not self.df._allow_adding_rows:
                    # adding rows is not allowed
                    if isinstance(key, (list, pd.Index)):
                        if not all(k in self.df.index for k in key):
                            raise ValueError(__NO_NEW_ROWS_MSG__)
                    elif key not in self.df.index:
                        raise ValueError(__NO_NEW_ROWS_MSG__)
                super(DataFrame, self.df).loc[key] = value

                # TODO we need to mark that every column has been modified, but only this row

        # Return the custom _LocIndexer instance
        return _LocIndexer(self)

    class AtIndexer:
        """this class override is used to intercept calls to at indexer of a pandas dataframe"""

        def __init__(self, obj):
            self.obj = obj

        def __getitem__(self, key):
            """intercept for the set operation"""

            return self.obj._get_value(*key)

        def __setitem__(self, key, value) -> None:
            """intercept for the set operation"""

            if (
                isinstance(value, pd.Series)
                and len(value) > len(self)
                and not self._allow_adding_rows
            ):
                raise ValueError(__NO_NEW_ROWS_MSG__)

            old_value = self.obj._get_value(*key)
            if value == old_value:
                # noop
                return

            rows = [key[0]]
            column = key[1]

            # Perform the actual setting operation
            self.obj._set_value(*key, value)

            # now update the DB. note that self is an AtIndexer
            # object, so we need to index into the pandas object
            if self.obj.auto_sync:
                # auto sync enabled, simply write ASAP
                self.obj.to_deeporigin()
            else:
                # auto sync not enabled, so we need to
                # keep track of changes in _modified_columns
                if column not in self.obj._modified_columns.keys():
                    # this is the first time we're modifying this column
                    self.obj._modified_columns[column] = set(rows)
                else:
                    # we've already modified this column before, so update the rows we're touched
                    self.obj._modified_columns[column].update(set(rows))

    @property
    def at(self):
        """Override the `at` property to return an AtIndexer"""
        return self.AtIndexer(self)

    def __setitem__(self, key, value):
        """Override the __setitem__ method to update the Deep Origin database when changes are made to the local
        dataframe. This method is called when an entire column is being updated"""

        # first, call the pandas method
        super().__setitem__(key, value)

        # now, update the Deep Origin database with the changes
        # we just made
        if self.auto_sync:
            self.to_deeporigin()
        else:
            # an empty set means "all "
            self._modified_columns[key] = set()

    def head(self, n=5):
        """Override the `head` method so that we don't display a spurious modified warning"""

        df = super().head(n)
        df._modified_columns = None
        return df

    def tail(self, n=5):
        """Override the `tail` method so that we don't display a spurious modified warning"""

        df = super().tail(n)
        df._modified_columns = None
        return df

    def append(
        self,
        other,
        ignore_index=False,
        verify_integrity=False,
        sort=False,
    ):
        """Override the `append` method"""
        if self._allow_adding_rows:
            return super().append(other, ignore_index, verify_integrity, sort)
        else:
            raise ValueError(__NO_NEW_ROWS_MSG__)

    def _repr_html_(self):
        """method override to customize printing in a Jupyter notebook"""

        # the entirety of this code is in a try/catch
        # block because pretty printing may fail,
        # and failure to pretty print should never
        # block other more mission critical code
        try:
            name = self.attrs["metadata"]["name"]
            hid = self.attrs["metadata"]["hid"]

            # placeholder URL, only used if something goes wrong
            url = "https://os.deeporigin.com/"

            # extract org name from url
            try:
                url = construct_resource_url(
                    name=hid,
                    row_type="database",
                )
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

            header = f'<h4 style="color: #808080;">Deep Origin / {org_name} / <a href = "{url}">{name} </a></h4>'
            txt = f'<p style="font-size: 12px; color: #808080;">Created {created_time_ago}. Row {self.attrs["last_updated_row"].hid} was last edited {edited_time_ago}'
            try:
                last_edited_by = get_last_edited_user_name(
                    self.attrs["last_updated_row"]
                )
                txt += "  by " + last_edited_by + ".</p>"
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
        except Exception:
            return super()._repr_html_()

    def __repr__(self):
        """method override to customize printing in an interactive session"""

        df_representation = super().__repr__()
        try:
            header = f'{self.attrs["metadata"]["hid"]}\n'
        except Exception:
            header = ""

        return header + df_representation

    @classmethod
    def from_deeporigin(
        cls,
        database_id: str,
        *,
        use_file_names: bool = True,
        reference_format: IDFormat = "human-id",
        client=None,
    ):
        """Create a local Deep Origin DataFrame from a Deep Origin database.

        Args:
            database_id (str): The ID of the Deep Origin database.
            use_file_names (bool, optional): Whether to use the file names in the Deep Origin database. Defaults to True.
            reference_format (IDFormat, optional): The format of the IDs in the Deep Origin database. Defaults to "human-id".

        """

        return api.get_dataframe(
            database_id=database_id,
            use_file_names=use_file_names,
            reference_format=reference_format,
            return_type="dataframe",  # we need this
            client=client,
        )

    def to_deeporigin(self):
        """Write data in dataframe to Deep Origin

        !!! tip "Deep Origin DataFrames can automatically synchronize"
            To automatically save changes to local DataFrames to Deep Origin databases, set the `auto_sync` attribute of the dataframe `True`.
        """

        columns = self._modified_columns.copy()

        for column in columns.keys():
            if column in ["Validation Status", "ID"]:
                continue

            rows = columns[column]

            column_metadata = [
                col for col in self.attrs["metadata"]["cols"] if col["name"] == column
            ]

            if len(column_metadata) == 0:
                # new column to add
                print("Creating new column...")
                # infer column type
                column_type = _infer_column_type(self[column])

                # make a new column
                response = api.add_database_column(
                    database_id=self.attrs["metadata"]["hid"],
                    type=column_type,
                    name=column,
                    key=column,
                )

                # add column metadata to column
                self.attrs["metadata"]["cols"].append(response["data"]["column"])
            else:
                # column already exists
                column_metadata = column_metadata[0]

                if column_metadata["type"] == "file":
                    continue

            if rows == set():
                # we're updating a whole column
                rows = list(self.index)

            api.set_data_in_cells(
                values=self[column],
                row_ids=rows,
                column_id=column,
                database_id=self.attrs["id"],
                columns=self.attrs["metadata"]["cols"],
            )

            print(f"✔︎ Wrote {len(rows)} rows in {column} to Deep Origin database.")

            # remove this from modified columns
            self._modified_columns.pop(column, None)

        # fetch info for the last modified row so we can update what we show
        try:
            self.attrs["last_updated_row"] = api.describe_row(row_id=rows[-1])
        except Exception:
            pass

    @property
    def _constructor(self):
        """this method overrides the _constructor property to return a DataFrame and is required for compatibility with a pandas DataFrame"""

        return DataFrame


@beartype
def _infer_column_type(column: pd.Series):
    """utility function to infer type of data in a pandas column"""

    non_null_values = column.dropna()

    if non_null_values.empty:
        return "string"

    # Helper functions for more complex checks
    def is_bool(val):
        if isinstance(val, str):
            val = val.lower()
            return val in ["true", "false"]
        return isinstance(val, bool)

    def is_int(val):
        try:
            return float(val).is_integer()
        except (ValueError, TypeError):
            return False

    def is_float(val):
        try:
            float(val)
            return True
        except (ValueError, TypeError):
            return False

    def is_date(val):
        try:
            parse(val, fuzzy=False)
            return True
        except (ValueError, TypeError):
            return False

    # Inference checks
    types = {
        "boolean": non_null_values.apply(is_bool).all(),
        "integer": non_null_values.apply(lambda x: is_int(x)).all(),
        "float": non_null_values.apply(lambda x: is_float(x)).all(),
        "date": non_null_values.apply(lambda x: is_date(x)).all(),
        "string": non_null_values.apply(lambda x: isinstance(x, str)).all(),
    }

    # Priority order of types (boolean -> integer -> float -> date -> string)
    for dtype, result in types.items():
        if result:
            return dtype

    # no clear type, fall back to string
    return "string"
