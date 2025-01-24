# Draco: a Chemical Data Extractor

Draco is a tool for extracting chemical data (Molecules, Reactions, Procedures, etc.) from documents. This version of Draco (1.0 Closed Alpha) includes the following capabilities:
- Process a single PDF
- Have no page or size limit for the PDF document
- Extract images of molecular fragments as SMILES
- Extract images of full molecules as SMILES
- Receive confidence score for each predicted SMILES
- Receive structural tokens for each predicted SMILES


## Draco Input Requirements
Draco will only accept files in PDF format. There is no file size restriction or page limit. The files have to be uploaded to the Data Hub (see below) for the script to be executed. 

Draco will recognize and analyze only images containing chemical structures - both fragments and full molecules.


## Draco Output Capabilities

![draco_ouput_example](https://s3.us-west-2.amazonaws.com/deeporigin.public/client_doc/draco_image2.png)
Draco produces a file in .xlsx (Excel) format. The file will contain predicted SMILES for each image of chemical structure that PDF document contains.

The report will contain one molecule per row. The following information will be provided with each molecule:
1. **Extracted Image** - an original image as it was reported in the PDF document
2. **Predicted Structure** - 2D rendering of the extracted SMILES string
3. **Confidence** - global confidence score for the entire molecule
4. **Confidence details** - local confidence score for each recognized structural element
5. **SMILES** - extracted SMILES string
6. **Source** - title of the uploaded PDF document
7. **Page** - page number where the extracted chemical structure was found

### Output Files

Draco produces a .xlsx (Excel) file that contains all the extracted chemical structures for a given document. Each row includes the extracted image, the predicted image, the predicted SMILES, the confidence score, and the confidence score associated with each token in the SMILES.

## Running Draco on Deep Origin

To run Draco on Deep Origin, follow these steps:
### 1. Create an account with Deep Origin
See details about how to a Deep Origin account [here](https://docs.deeporigin.io/docs/users)
### 2. Create a database to store input and output files

Create a Column containing the PDF Document and an output column Result, which will store the output files. The type of these columns is “File.”

![draco_database_example](https://github.com/user-attachments/assets/926a4f06-3c27-4b4b-9b47-79fc98e96723)

### 3. Start a tool run on Deep Origin

To start a tool run, use:

```python
from deeporigin.tools import run

job_id = run.draco(
    database_id="<your-db-name>",
    row_id="<row-name>",
    input_file_column_name="<input-file-column-name>"
    output_column_name="<output-column-name>"
)
```

`run.draco` returns the ID of the tool run, which can be used to monitor the status of the run and terminate it if needed. 

`run.draco` prints a message that looks like:

```bash
🧬 Job started with ID: 9f7a3741-e392-45fb-a349-804b7fca07d7
```

To monitor the status of the tool run, use:

```python
from deeporigin.tools.utils import query_run_status
query_run_status("9f7a3741-e392-45fb-a349-804b7fca07d7")
```

To wait for the tool run to finish, use:

```python
from deeporigin.tools.utils import wait_for_job
wait_for_job("9f7a3741-e392-45fb-a349-804b7fca07d7")
```

## Example
Using the database shown above in section "2. Create a database to store input and output files", the code to extract molecules from the uploaded document (patent_test.pdf) is:
```python
from deeporigin.tools import run

job_id = run.draco(
    database_id="draco_use_case",
    row_id="draco-use-case-1",
    input_file_column_name="Document"
    output_column_name="Result"
)
```
