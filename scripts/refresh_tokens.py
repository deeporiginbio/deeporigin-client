"""small script to refresh tokens"""


def main():
    from deeporigin import auth

    auth.get_tokens(refresh=True)


if __name__ == "__main__":
    main()
