"""this module contains tests for the DataFrame class"""

import numpy as np
from deeporigin.data_hub import api
from deeporigin.data_hub.dataframe import DataFrame

# we use this for testing
DB_NAME = "test-db-client-n23"


# def test_dataframe_float():
#     api.create_database(name=DB_NAME)

#     api.add_database_column(
#         database_id=DB_NAME,
#         name="float",
#         key="float",
#         type="float",
#     )

#     response = api.make_database_rows(database_id=DB_NAME, n_rows=10)
#     row_ids = [row.hid for row in response.rows]

#     values = np.random.random(10)
#     api.set_data_in_cells(values=values, row_ids=row_ids)
