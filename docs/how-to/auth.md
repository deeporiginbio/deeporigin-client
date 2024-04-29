You can sign into the Deep Origin platform using the Python Client.

```python
from deeporigin import sign_into_do_platform
sign_into_do_platform()
```

You will be presented with a prompt that looks like:

```
To connect to the Deep Origin platform, navigate your browser to 

https://<env>auth0.com/activate?user_code=VMPZ-PQFG

and verify the confirmation code is "VMPZ-PQFG", and click the "Confirm" button.
```
Visting that URL, you will see a prompt that looks like:

![](../images/auth-code.png)

On clicking the `Confirm` button, you are shown:

![](../images/auth-confirm.png)

and the python client will now return an access and refresh 
token. These tokens are also cached to disk and will automatically
be used in subsequent interactions. 

!!! info "Authenticating"
    You do not need to authenticate every time you use the client or the CLI. Authenticating once, before first use, should be sufficient.

