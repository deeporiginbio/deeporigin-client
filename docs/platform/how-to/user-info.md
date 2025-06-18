# Get user information 

## Get info about current user

To get information about the currently logged in user, including the user ID, use:

```{.python notest}
from deeporigin.platform import users_api
users_api.get_account()
```

returns information about the current user. A typical response is:

```json
{
    "attributes": {
        "company": null,
        "email": "curie@deeporigin.com",
        "expertise": null,
        "industries": null,
        "name": "Marie Curie",
        "pendingInvites": [],
        "platform": "OS",
        "title": null
    },
    "id": "google-apps|curie@deeporigin.com",
    "type": "User"
}
```
