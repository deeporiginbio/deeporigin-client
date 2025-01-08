# Draco: a Chemical Data Extractor

[Draco] is a tool that aims to extract chemical data (Molecules, Reactions, procedures, etc.) from documents. In its first version, Draco focuses at extracting molecules (fragments and complete structures) from PDF documents.

## File Inputs

### 1. PDF File

Upload a PDF file containing chemical structures to extract to the Data Hub.

## Outputs

### Output Files

Draco produces a .xlsx (Excel) file that contains all the extracted chemical structures for a given document. Each row includes the extracted image, the predicted image, the predicted SMILES, the confidence score, and the confidence score associated with each token in the SMILES.

## Running Draco on Deep Origin

To run Draco on Deep Origin, follow these steps:

### 1. Create a database to store input and output files

Create a Column containing the PDF files and an output column, which will store the output files. The type of these columns is File.

![draco_database_example](https://github.com/user-attachments/assets/926a4f06-3c27-4b4b-9b47-79fc98e96723)

### 2. Start a tool run on Deep Origin

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
Using the database shown above in section "1. Create a database to store input and output files", the code to extract molecules from the uploaded document (patent_test.pdf) is:
```python
from deeporigin.tools import run

job_id = run.draco(
    database_id="draco_use_case",
    row_id="draco-use-case-1",
    input_file_column_name="Document"
    output_column_name="Result"
)
```
