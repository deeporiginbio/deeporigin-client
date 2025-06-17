"""small script to refresh tokens"""

from deeporigin import auth

auth.get_tokens(refresh=True)
