"""
This module defines a class called DataFrame that is a drop-in
replacement for a pandas DataFrame, but also allows for easy
updating of Deep Origin databases.
"""

from datetime import datetime, timezone

from beartype import beartype
from beartype.typing import Optional
from dateutil.parser import parse
import humanize
import pandas as pd

from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils.config import construct_resource_url
from deeporigin.utils.constants import DataType, IDFormat
from deeporigin.utils.network import check_for_updates

check_for_updates()


__NO_NEW_ROWS_MSG__ = "Adding rows to Deep Origin DataFrames is not allowed. "
__NO_NEW_ROWS_FIX__ = (
    "If you want to add rows to the underlying database, use `api.add_database_rows()`."
)


class DataFrame(pd.DataFrame):
    """A subclass of pandas DataFrame that allows for easy updating of a Deep Origin database. This can be used as a drop-in replacement for a pandas DataFrame, and should support all methods a pandas DataFrame supports.

    The primary method of creating an object of this type is to use the [from_deeporigin][src.data_hub.dataframe.DataFrame.from_deeporigin] class method.
    """

    auto_sync: bool = False
    """When `True`, changes made to the dataframe will be automatically synced to the Deep Origin database this dataframe represents."""

    _modified_columns: dict = dict()
    """if data is modified in a dataframe, and auto_sync is False, this list will contain the columns that have been modified so that the Deep Origin database can be updated. If an empty list, the Deep Origin database will not be updated, and the dataframe matches the Deep Origin database at the time of creation."""

    def _track_changes(self, column: str, rows: list):
        """callback that tracks changes made to the DB, and responds appropriately. if auto_sync is true, changes
        are written immediately to DB. if not, then they're tracked in _modified_columns

        Args:
            column (str): the name of the column that was modified
            rows (list): the IDs of rows that were modified
        """

        if self.auto_sync:
            # auto sync enabled, simply write ASAP
            self.to_deeporigin()
        else:
            # auto sync not enabled, so we need to
            # keep track of changes in _modified_columns
            if column not in self._modified_columns.keys():
                # this is the first time we're modifying this column
                self._modified_columns[column] = set(rows)
            else:
                # we've already modified this column before, so update the rows we're touched
                self._modified_columns[column].update(set(rows))

    @property
    def loc(self):
        class _LocIndexer:
            def __init__(self, df):
                self.df = df

            def __getitem__(self, key):
                """this function is called when we slice a DB, so we need to disallow appending rows"""

                # first call the superclass method
                df = super(DataFrame, self.df).loc[key]

                # inherit attributes
                df.attrs = self.df.attrs

                df._modified_columns = self.df._modified_columns

                return df

            def __setitem__(self, key, value):
                """callback for adding a new row or modifying data in existing rows"""

                # disallow making new rows
                if isinstance(key, (list, pd.Index)):
                    rows = key
                    if not all(k in self.df.index for k in key):
                        raise DeepOriginException(
                            message=__NO_NEW_ROWS_MSG__,
                            fix=__NO_NEW_ROWS_FIX__,
                        )
                else:
                    rows = [key]
                    if key not in self.df.index:
                        raise DeepOriginException(
                            message=__NO_NEW_ROWS_MSG__,
                            fix=__NO_NEW_ROWS_FIX__,
                        )

                # first check if the new value is the same as
                # the old value
                old_value = list(self.df.loc[key])

                try:
                    if value == old_value:
                        # noop
                        return
                except Exception:
                    pass

                super(DataFrame, self.df).loc[key] = value

                for col in self.df.columns:
                    self.df._track_changes(col, rows)

        # Return the custom _LocIndexer instance
        return _LocIndexer(self)

    class AtIndexer:
        """this class override is used to intercept calls to at indexer of a pandas dataframe"""

        def __init__(self, df):
            self.df = df

        def __getitem__(self, key):
            """intercept for the set operation"""

            return self.df._get_value(*key)

        def __setitem__(self, key, value) -> None:
            """intercept for the set operation"""

            if isinstance(value, pd.Series) and len(value) > len(self):
                raise DeepOriginException(
                    title="Adding rows to a DataFrame not allowed",
                    message=__NO_NEW_ROWS_MSG__,
                    fix=__NO_NEW_ROWS_FIX__,
                )

            old_value = self.df._get_value(*key)

            # the reason this is in a try block is because
            # this can fail for any number of reasons.
            # for example, types of the two things may be
            # different, or the old or new value may be missing,
            # in which case an equality is meaningless.
            try:
                if value == old_value:
                    # noop
                    return
            except Exception:
                pass

            rows = [key[0]]
            column = key[1]

            # Perform the actual setting operation
            self.df._set_value(*key, value)

            # now update the DB.
            self.df._track_changes(column, rows)

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
            if isinstance(key, str):
                # an empty set means "all "
                self._modified_columns[key] = set()
            elif isinstance(key, list):
                for item in key:
                    self._modified_columns[item] = set()

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
        raise DeepOriginException(
            title="Adding rows to a DataFrame not allowed",
            message=__NO_NEW_ROWS_MSG__,
            fix=__NO_NEW_ROWS_FIX__,
        )

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

            date_str = self.attrs["last_updated_row"].dateUpdated
            date_obj = datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S.%f").replace(
                tzinfo=timezone.utc
            )
            edited_time_ago = humanize.naturaltime(now - date_obj)

            header = f'<h4 style="color: #808080;">Deep Origin / {org_name} / <a href = "{url}">{name} </a></h4>'
            txt = f'<p style="font-size: 12px; color: #808080;">Created {created_time_ago}. Row {self.attrs["last_updated_row"].hid} was last edited {edited_time_ago}'

            txt += ".</p>"

            if self._modified_columns:
                txt += '<p style="color: #808080; font-size: 12px">⚠️ This dataframe contains changes that have not been written back to the Deep Origin database.</p>'
            elif self.auto_sync:
                txt += '<p style="color: #808080; font-size: 12px">🧬 This dataframe will automatically write changes made to it back to Deep Origin.</p>'
            df_html = super()._repr_html_()
            return header + txt + df_html
        except Exception as error:
            print(error)
            return super()._repr_html_()

    def __repr__(self):
        """method override to customize printing in an interactive session"""

        df_representation = super().__repr__()
        try:
            header = f"{self.attrs['metadata']['hid']}\n"
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
        filter: Optional[dict] = None,
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
            filter=filter,
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
                )

                # add column metadata to column
                self.attrs["metadata"]["cols"].append(response.column)
            else:
                # column already exists
                column_metadata = column_metadata[0]

                if column_metadata["type"] == "file":
                    continue

            if rows == set():
                # we're updating a whole column
                rows = list(self.index)

            api.set_data_in_cells(
                values=self[column][list(rows)].to_list(),
                row_ids=list(rows),
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
def _infer_column_type(column: pd.Series) -> DataType:
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
        "text": non_null_values.apply(lambda x: isinstance(x, str)).all(),
    }

    # Priority order of types (boolean -> integer -> float -> date -> string)
    for dtype, result in types.items():
        if result:
            return dtype

    # no clear type, fall back to string
    return "text"
