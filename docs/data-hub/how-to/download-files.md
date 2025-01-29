# Download files

This page describes how to download files from Deep Origin to your local computer. 


## Download one or many files from the Data hub

To download file(s) to the Deep Origin data hub, run the following commands:

=== "CLI"

    ```bash
    deeporigin data download-files
    ```

    This will download all files on Deep Origin to the current folder. 

    To download files that have been assigned to a particular row, use:

    ```bash
    deeporigin data download-files --assigned-row-ids <row-id-1>  <row-id-2> ...
    ```

    To download specific files, pass the file IDs using:


    ```bash
    deeporigin data download-files --file-ids <file-1> <file-1> ...
    ```
    

=== "Python"

    ```py
    from deeporigin.data_hub import api
    api.download_files(files)
    ```

    `files` is a list of files to download, and is a list of `ListFilesResponse` objects. To obtain this list, use `api.list_files()`, the output of which can be used as an input to `download_files`. 

    !!! Tip "Download all files"
        To download all files, call `api.list_files()` and pass the output to `download_files`.


