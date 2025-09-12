# Sign into Deep Origin

!!! tip "Configure if running locally"
    If you're running this code on your local computer (outside of a Deep Origin Workstation), make sure to [configure](../configure.md#on-your-local-computer) it first.

To use most of the functionality of the Python client, you must first run the following commands to sign into Deep Origin.





```{.python notest}
from deeporigin import auth
auth.authenticate()
```

You will be presented with a prompt similar to below:

```shell
To connect to the Deep Origin OS, navigate your browser to 

https://<env>auth0.com/activate?user_code=VMPZ-PQFG

and verify the confirmation code is "VMPZ-PQFG", and click the "Confirm" button.
```

When you visit that URL, you will see a prompt that looks like:

![](../images/auth-code.png)

After clicking the `Confirm` button, you will see a confirmation similar to below:

![](../images/auth-confirm.png)

After signing in, your access tokens will be cached to disk and then automatically
be used in subsequent interactions with Deep Origin.

!!! info "Authenticating"
    In most cases, you only need to authenticate to the Deep Origin OS once.
    You do not need to authenticate every time you use the client or the CLI.

!!! question "Authenticating on Deep Origin workstations"
    Presently, workstation users must authenticate (once) to the Deep Origin OS. We plan to develop the capability to automatically authenticate workstation users.
