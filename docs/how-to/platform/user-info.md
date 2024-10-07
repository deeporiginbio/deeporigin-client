# Get user information 

## Get info about current user

To get information about the currently logged in user, including the user ID, use:

```python
from deeporigin.platform import api
api.whoami()
```

returns information about the current user. A typical response is:

```json
{
    "data": {
        "attributes": {
            "company": null,
            "expertise": null,
            "industries": null,
            "pendingInvites": [],
            "platform": "OS",
            "title": null
        },
        "id": "google-apps|user@deeporigin.com",
        "type": "User"
    },
    "links": {
        "self": "https://os.deeporigin.io/users/me"
    }
}
```

## Get information about a user 

To get information about a user, use:


```python
from deeporigin.platform import api
api.resolve_user("user-id")
```

where `user-id` is the ID of the user, in the format returned by `api.whoami()`. A typical response looks like:


```json
{
    "data": {
        "attributes": {
            "avatar": "https://...",
            "email": "user@deeporigin.com",
            "name": "User Name"
        },
        "id": "918ddd25-ab97-4400-9a14-7a8be1216754",
        "type": "User"
    },
    "links": {
        "self": "https://..."
    }
}
```