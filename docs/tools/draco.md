# Draco: a Chemical Data Extractor

[Draco] is a tool which aims to extract chemical data (Molecules, Reactions, procedures, etc.) from documents. In its first version, Draco focuses at extracting molecules (fragments and complete structures) from PDF documents

## File Inputs

### 1. PDF File

Upload a PDF file containing chemical structures to extract to the Data Hub.

## Parameters

This section describes the parameters for a tool run, that are passed in the `start_run` function.

### IUPAC name (OPTIONAL) - NOT YET AVAILABLE

Use a model to predict the IUPAC name of the predicted structures

## Outputs

### Output Files

Draco produces an .xlsx (Excel) file which contains all the extracted chemical strucutres for a given document. Each row contains the extracted image, the predicted image, the predicted SMILES, the confidence score and the confidence score associated to each token in the SMILES

## Running Draco on Deep Origin

To run Draco on Deep Origin, follow these steps:

### 1. Create a database to store input and output files

### 2. Start a tool run on Deep Origin

For this tool run, we will use the following parameters:

```python
parameters = {
    "input_file": 'your_pdf_file',
}
```

To start a tool run, use:

```python
from deeporigin.tools import draco

job_id = draco.start_run(
    database_id="<your-db-name>",
    row_id="<row-name>",
    input_file_column_name="<input-file-column-name>"
    output_column_name="<output-column-name>"
)
```

`start_run` returns the ID of the tool run, that can be used to monitor the status of the run and terminate it if needed. 

`start_run` prints a message that looks like:

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

## RoadMap

- Add Reaction extraction
- Extraction from Text
- TODO

## TODO

Show a template picture of the data hub