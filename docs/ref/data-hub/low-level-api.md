# API reference: Low-level API

## Background

We use automatically generated code to implement a low-level Python SDK that connects to the Deep Origin data hub API. This page describes that API. In most cases, you will not need to know how this works, nor will you need to use it.


!!! tip "Using methods in these classes"
    The `deeporigin` package provides a wrapper that allows you to use methods in the `DeepOrigin` class, without needing to instantiate an object or configure it. Every method in this class is wrapped and exposed as a function in the `_api` module, and can be therefore used directly.

    For example, to use the [list_files](#deeporigin_data._client.DeeporiginData.list_files) method in the `DeeporiginData` class, you can:

    ```python
    from deeporigin.data_hub import _api
    _api.list_files()
    ```

    Your IDE should be able to provide contextual help by pressing tab after `_api.list_files(` and give you information about types and parameters of all methods listed here. 

:::deeporigin_data._client
    options:
      show_root_toc_entry: false
      show_root_heading: false
      show_category_heading: false
      members_order: alphabetical
      show_object_full_path: false
