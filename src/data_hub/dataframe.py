from typing import Optional

import pandas as pd
from deeporigin.data_hub import api


class DataFrame(pd.DataFrame):
    def _repr_html_(self):
        header = "<h3>Deep Origin DataFrame</h3>"
        df_html = super()._repr_html_()
        return header + df_html

    def __repr__(self):
        header = "Deep Origin DataFrame\n"
        df_representation = super().__repr__()
        return header + df_representation

    def sync(self, *, columns: Optional[list] = None):
        print("Syncing...")

        for column in self.columns:
            if column == "Validation Status":
                continue

            print(column)

            api.set_data_in_cells(
                values=self[column],
                row_ids=list(self.index),
                column_id=column,
                database_id=self.attrs["id"],
            )

    @property
    def _constructor(self):
        return DataFrame
