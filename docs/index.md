
# Deep Origin Python client

<div class="grid cards" markdown>

- :fontawesome-solid-handshake-simple: **Use Deep Origin Tools**  <br>Run Docking, molprops, FEP and other molecular modeling tools
- :octicons-unlock-24: **Free and open-source**
<br>Install onto your computer to use your data, variables, secrets, and other resources outside Deep Origin.
- :material-download: **Easy to install**
<br>Just run `pip install deeporigin`.
- :fontawesome-brands-python: **Pure Python**
<br>Lightweight, written in pure Python. Works on any system that can run Python.

</div>

## Example

The Deep Origin CLI and Python client allow you to programmatically
interact with the [Deep Origin OS :octicons-link-external-16:](https://os.deeporigin.io/).

For example:


```python
from deeporigin.drug_discovery import Ligand

ligand = Ligand.from_identifier("serotonin")

ligand.show()
```

<iframe 
    src="./dd/how-to/serotonin.html" 
    width="100%" 
    height="600" 
    style="border:none;"
    title="Visualization of single ligand (serotonin)"
></iframe>