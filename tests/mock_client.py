from deeporigin_data import types
from pydantic import BaseModel

USER_DRN = "drn:identity::user:google-apps|user@deeporigin.com"


class GenericModel(BaseModel):
    """this catch all model is used for some problematic
    models that cannot be converted into dict and cannot be constructed due to errors in the generated code"""

    class Config:
        extra = "allow"


def make_fields(n: int = 1):
    return [
        {
            "columnId": f"_column:col-{idx}",
            "cellId": "_cell:ubpEc7OdSkVe6AkVI3vi0",
            "validationStatus": "valid",
            "type": "text",
            "value": "placeholder text",
        }
        for idx in range(n)
    ]


def column_item(idx: int = 1):
    return {
        "id": f"_column:col-{idx}",
        "name": f"col-{idx}",
        "key": f"col-{idx}",
        "parentId": "db-placeholder",
        "type": "text",
        "cardinality": "one",
    }


def list_files_response(assignments=None):
    return [
        types.list_files_response.Data(
            file=types.list_files_response.DataFile(
                id="_file:09fwpdPqHtVdt4jj43ywC",
                contentLength=1343642.0,
                dateCreated="2024-07-03 23:45:44.214",
                name="foo.pdf",
                status="ready",
                uri="s3://data.deeporigin-com.foo/files/_file:09fwpdPqHtVdt4jj43ywC",
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

    folders = ["ws-placeholder"]
    databases = ["db-placeholder"]
    rows = [f"row-{idx}" for idx in range(10)]
    columns = [f"column-{idx}" for idx in range(5)]
    file = list_files_response()[0].file

    def list_rows(self, filters):
        """callback when we invoke ListRows"""

        if filters == []:
            # list_rows called with no filters. return
            # a folder, a database, and some rows

            # root folder
            rows = [
                types.list_rows_response.Data(
                    id=self.folders[0],
                    hid=self.folders[0],
                    type="workspace",
                    name=self.folders[0],
                    parent_id=None,
                )
            ]
            # root database
            rows += [
                types.list_rows_response.Data(
                    hid=self.databases[0],
                    id=self.databases[0],
                    type="database",
                    name=self.folders[0],
                    parent_id=self.folders[0],
                )
            ]

            # rows in the database
            rows += [
                types.list_rows_response.Data(
                    hid=f"row-{idx}",
                    id=f"row-{idx}",
                    type="row",
                    parentId=self.databases[0],
                )
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
            # return the root folder
            return [
                types.list_rows_response.Data(
                    hid=self.folders[0],
                    id=self.folders[0],
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
                    parentId=self.folders[0],
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

    def list_files(self, filters: list = []):
        assignments = None
        if filters == [dict(is_unassigned=False)]:
            # return assigned files, so we need assignments
            assignments = [
                types.list_files_response.DataAssignment(
                    rowId="_row:HmeXSzMWo96G0Or6HpEXK"
                )
            ]
        return list_files_response(assignments=assignments)

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
        parentId = kwargs["database"]["parentId"]
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

    def describe_database_stats(self, **kwargs):
        return types.describe_database_stats_response.Data(rowCount=5)

    def convert_id_format(self, **kwargs):
        return [
            types.convert_id_format_response.Data(
                id="_workspace:57xVgPwGWdhxP5DoQGT30", hid="sandbox"
            )
        ]

    def describe_row(self, row_id, fields: bool = True):
        if row_id.startswith("row-"):
            row_type = "row"
            parent_id = "db-placeholder"
        elif row_id.startswith("db-"):
            row_type = "database"
            parent_id = "ws-placeholder"
        elif row_id.startswith("ws-"):
            row_type = "workspace"
            parent_id = None

        cols = [column_item(idx=idx) for idx in range(5)]

        response = GenericModel(
            id=row_id,
            hid=row_id,
            hid_prefix=None,
            name=row_id,
            type=row_type,
            cols=cols,
            creation_block_id=None,
            creation_parent_id=None,
            editor=None,
            is_template=None,
            row_json_schema=None,
            submission_status=None,
            parent_id=parent_id,
            dateCreated="2024-06-20 19:21:10.81895",
            dateUpdated="2024-07-17 15:33:35.951",
            createdByUserDrn=USER_DRN,
            editedByUserDrn=USER_DRN,
            validationStatus="invalid",
        )

        if fields:
            response.fields = make_fields(5)
        return response

    def delete_database(self, **kwargs):
        return dict()

    def delete_workspace(self, **kwargs):
        return dict()

    def list_mentions(self, **kwargs):
        return types.list_mentions_response.Data(
            mentions=[
                types.list_mentions_response.DataMention(
                    id="_row:k0okQigiru5g6jc8IFfEj",
                    hid="ligand-1",
                    type="row",
                    name="Mol 1",
                    parent_id=None,
                )
            ]
        )

    def list_database_rows(self, database_row_id):
        return [
            types.list_database_rows_response.Data(
                id=database_row_id,
                hid="ligand-4",
                type="row",
                creation_block_id=None,
                creation_parent_id=None,
                fields=[
                    types.list_database_rows_response.DataFieldUnionMember0(
                        cellId="_cell:V0arbWA6UvQ4Q8UOw5Jua",
                        columnId="_column:zRpsD9FirjyIkaszUUTku",
                        type="text",
                        validationStatus="valid",
                        invalid_data=None,
                        system_type="name",
                        value="Mol 4",
                    ),
                    types.list_database_rows_response.DataFieldUnionMember0(
                        cellId="_cell:kJ97kqhgilDGPYJWEzikf",
                        columnId="_column:iooVcZhMFXKYMQCeirREJ",
                        type="text",
                        validationStatus="valid",
                        invalid_data=None,
                        system_type=None,
                        value="CCC",
                    ),
                ],
                is_template=None,
                name="Mol 4",
                parent_id="_database:ZiEtc7k02S8XIVOG5z007",
                submission_status=None,
                parentId="_database:ZiEtc7k02S8XIVOG5z007",
                dateCreated="2024-07-10 19:26:52.597836",
                dateUpdated="2024-07-10 19:29:33.341",
                createdByUserDrn=USER_DRN,
                editedByUserDrn=USER_DRN,
                validationStatus="valid",
            )
        ]

    def create_file_download_url(self, **kwargs):
        return types.create_file_download_url_response.Data(
            downloadUrl="https://github.com/deeporiginbio/deeporigin-client/archive/refs/tags/0.0.3.zip"
        )

    def create_file_upload(
        self,
        *,
        name,
        content_type: str,
        content_length: str,
        **kwargs,
    ):
        file = types.create_file_upload_response.DataFile(
            id="_file:placeholder",
            contentLength=content_length,
            name=name,
            status="ready",
            uri="s3://data.deeporigin-com.ijvjf/files/placeholder",
            contentType=content_type,
            dateCreated="2024-07-25T15:24:09.459Z",
        )

        return types.create_file_upload_response.Data(
            file=file,
            uploadUrl="https://foo.bar/foo",
        )

    def describe_file(self, file_id: str):
        return types.describe_file_response.Data(
            id=file_id,
            content_length=56.0,
            name="local database.csv",
            status="ready",
            uri=f"s3://data.deeporigin-com.ijvjf/files/{file_id}",
            contentLength=56,
            contentType="application/vnd.ms-excel",
            dateCreated="2024-07-03 23:42:51.746",
            dateUpdated="2024-07-03 23:42:51.746",
            createdByUserDrn=USER_DRN,
        )
