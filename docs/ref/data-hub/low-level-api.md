# API reference: Low-level API

## Background

We use automatically generated code to implement a low-level Python SDK that connects to the Deep Origin data hub API. This page describes that API. In most cases, you will not need to know how this works, nor will you need to use it.

!!! tip "Using the methods in the DeeporiginData class"
    The `deeporigin._api` module provides wrappers that allow you to use methods in the `deeporigin_data.DeeporiginData` class, without needing to instantiate an object and configure it. Each method of this class is exposed as a function with the same name. These functions enable you to easily use the methods of the `deeporigin_data.DeeporiginData` class.

    For example, to use the [list_files](#deeporigin_data._client.DeeporiginData.list_files) method in the `DeeporiginData` class, run:

    ```python
    from deeporigin.data_hub import _api
    _api.list_files()
    ```

    Your IDE should be able to provide information about the arguments and return types of these functions by typing `_api.list_files(` and then pressing tab.    

:::deeporigin_data._client
    options:
      show_root_toc_entry: false
      show_root_heading: false
      show_category_heading: false
      members_order: alphabetical
      show_object_full_path: false
