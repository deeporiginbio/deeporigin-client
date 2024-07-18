from deeporigin_data import types
from pydantic import BaseModel

USER_DRN = "drn:identity::user:google-apps|user@deeporigin.com"


class GenericModel(BaseModel):
    """this catch all model is used for some problematic
    models that cannot be converted into dict and cannot be constructed due to errors in the generated code"""

    pass


def list_files_response(assignments=None):
    return [
        types.list_files_response.Data(
            file=types.list_files_response.DataFile(
                id="_file:09fwpdPqHtVdt4jj43ywC",
                contentLength=1343642.0,
                dateCreated="2024-07-03 23:45:44.214",
                name="science.ade2574_sm.pdf",
                status="ready",
                uri="s3://data.deeporigin-com.ijvjf/files/_file:09fwpdPqHtVdt4jj43ywC",
                content_type="application/pdf",
                created_by_user_drn=USER_DRN,
                date_updated="2024-07-03 23:45:44.214",
            ),
            assignments=assignments,
        )
    ]


class MockClient:
    """mock client to respond with static data for testing
    purposes"""

    workspaces = ["ws-placeholder"]
    databases = ["db-placeholder"]
    rows = [f"row-{idx}" for idx in range(10)]
    columns = [f"column-{idx}" for idx in range(5)]
    file = list_files_response()[0].file

    def list_rows(self, filters):
        """callback when we invoke ListRows"""

        if filters == []:
            # list_rows called with no filters. return
            # a workspace, a database, and some rows

            # root workspace
            rows = [
                types.list_rows_response.Data(
                    id=self.workspaces[0],
                    hid=self.workspaces[0],
                    type="workspace",
                    name=self.workspaces[0],
                    parent_id=None,
                )
            ]
            # root database
            rows += [
                types.list_rows_response.Data(
                    hid=self.databases[0],
                    id=self.databases[0],
                    type="database",
                    name=self.workspaces[0],
                    parent_id=self.workspaces[0],
                )
            ]

            # rows in the database
            rows += [
                types.list_rows_response.Data(
                    hid=f"row-{idx}",
                    id=f"row-{idx}",
                    type="row",
                    parentId=self.databases[0],
                ).model_dump()
                for idx in range(10)
            ]
            return rows

        elif (
            "parent" in filters[0].keys()
            and "id" in filters[0]["parent"].keys()
            and filters[0]["parent"]["id"].startswith("db-")
        ):
            # we're asking for rows that belong to a database

            return [
                types.list_rows_response.Data(
                    hid=f"row-{idx}",
                    id=f"row-{idx}",
                    type="row",
                    parentId=filters[0]["parent"]["id"],
                )
                for idx in range(10)
            ]

        elif filters == [dict(parent=dict(is_root=True))] or filters == [
            dict(row_type="workspace")
        ]:
            # return the root workspace
            return [
                types.list_rows_response.Data(
                    hid=self.workspaces[0],
                    id=self.workspaces[0],
                    type="workspace",
                    parentId=None,
                )
            ]
        elif filters == [dict(row_type="database")]:
            # return the example database
            return [
                types.list_rows_response.Data(
                    hid=self.databases[0],
                    id=self.databases[0],
                    type="database",
                    parentId=self.workspaces[0],
                )
            ]
        elif filters == [dict(row_type="row")]:
            # return the example row
            return [
                types.list_rows_response.Data(
                    hid=self.rows[0],
                    id=self.rows[0],
                    type="row",
                    parentId=self.databases[0],
                )
            ]
        else:
            print(filters)
            raise NotImplementedError(
                "This specific request type hasn't been mocked yet"
            )

    def invoke_describe_row(self, data):
        """callback when we invoke DescribeRow"""

        if data["rowId"].startswith("db-"):
            # we are likely asking for a database
            name = self.databases[0]
            return DescribeRowResponseDatabase(
                id=name,
                hid=name,
                name=name,
                type="database",
                parentId="ws_placeholder",
                cols=[
                    dict(
                        ColumnItem(
                            id=f"column-{idx}",
                            key=f"column-{idx}",
                        )
                    )
                    for idx in range(5)
                ],
                hidPrefix="placeholder",
            ).model_dump()

        else:
            # we are asking for a row in a database
            name = data["rowId"]
            fields = [
                FieldItem(
                    columnId=f"column-{idx}",
                    cellId=f"column-{idx}",
                )
                for idx in range(5)
            ]
            row = DescribeRowResponseRow(
                id=name,
                hid=name,
                parentId=self.databases[0],
                fields=fields,
            ).model_dump()

            if not data["fields"]:
                row.pop("fields", None)

        return row

    def invoke_list_database_rows(self, data):
        """callback when we invoke ListDatabaseRows"""

        fields = [
            FieldItem(
                columnId=f"column-{idx}",
                cellId=f"column-{idx}",
            )
            for idx in range(5)
        ]
        return [
            DescribeRowResponseRow(
                id=f"row-{idx}",
                hid=f"row-{idx}",
                parentId=data["databaseRowId"],
                fields=fields,
            ).model_dump()
            for idx in range(5)
        ]

    def list_files(self, filters):
        if filters == [dict(isUnassigned=True)]:
            return list_files_response()
        elif filters == [dict(isUnassigned=False)]:
            list_files_response(
                assignments=[
                    types.list_files_response.DataAssignment(
                        rowId="_row:HmeXSzMWo96G0Or6HpEXK"
                    )
                ],
            )

    def ensure_rows(self, **kwargs):
        """callback when we invoke EnsureRows"""
        return types.ensure_rows_response.Data(
            **{
                "rows": [
                    {
                        "id": "_row:RETrlQNnPAPNSC8q21toF",
                        "parentId": "_database:IpWYEnJE0uIx2TchRVVLQ",
                        "type": "row",
                        "dateCreated": "2024-07-03 13:54:09.888684",
                        "dateUpdated": "2024-07-03 13:54:09.888684",
                        "createdByUserDrn": USER_DRN,
                        "hid": "dna-12",
                        "validationStatus": "invalid",
                    }
                ]
            }
        )

    def create_workspace(self, **kwargs):
        name = kwargs["workspace"]["name"]
        return types.create_workspace_response.Data(
            type="workspace",
            id="_workspace:" + name,
            hid="_workspace:" + name,
            name=name,
            dateCreated="2024-05-08 19:16:09.26078",
            dateUpdated="2024-05-08 19:16:09.26078",
            createdByUserDrn=USER_DRN,
        )

    def create_database(self, **kwargs):
        name = kwargs["database"]["name"]
        parentId = kwargs["database"]["parent_id"]
        return GenericModel(
            id="_database:FzTJKQ11i1VVRfqolWsFg",
            date_created="2024-07-18 13:36:18.017657",
            hid=name,
            hid_prefix=name,
            name=name,
            type="database",
            cols=None,
            created_by_user_drn=USER_DRN,
            creation_block_id=None,
            creation_parent_id=None,
            date_updated="2024-07-18 13:36:18.017657",
            edited_by_user_drn=None,
            editor=None,
            is_template=None,
            parent_id=parentId,
            row_json_schema={"type": "object", "required": [], "properties": {}},
            submission_status=None,
            validation_status=None,
        )

    def delete_database(self, **kwargs):
        return dict()

    def invoke(self, endpoint: str, data: dict):
        """simple returns data without making any network requests"""

        if endpoint == "ListRows":
            return self.invoke_list_rows(data)

        elif endpoint == "EnsureRows":
            return self.invoke_ensure_rows(data)

        elif endpoint == "DescribeRow":
            return self.invoke_describe_row(data)

        elif endpoint == "ConvertIdFormat":
            return [{"id": "_row:W6DjtaCrZ201EGLpmZtGO", "hid": "sample-1"}]

        elif endpoint == "ListDatabaseRows":
            return self.invoke_list_database_rows(data)

        elif endpoint == "DescribeFile":
            return file_description()

        elif endpoint == "DescribeDatabaseStats":
            return {"rowCount": 5}

        elif endpoint == "CreateFileDownloadUrl":
            return {
                "downloadUrl": "https://github.com/deeporiginbio/deeporigin-client/archive/refs/tags/0.0.3.zip"
            }

        elif endpoint == "ListMentions":
            return {
                "mentions": [
                    {
                        "type": "row",
                        "id": "_row:W6DjtaCrZ201EGLpmZtGO",
                        "hid": "sample-1",
                    }
                ]
            }
