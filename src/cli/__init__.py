import arguably
from deeporigin import auth


@arguably.command
def variables__install():
    pass


@arguably.command
def variables__auto_install():
    pass


@arguably.command
def variables__uninstall():
    pass


@arguably.command
def data():
    """Commands to interact with the Data Hub"""
    pass


@arguably.command
def data__upload():
    pass


@arguably.command
def data__download():
    pass


@arguably.command
def data__delete():
    pass


@arguably.command
def data__describe():
    pass


@arguably.command
def data__new():
    pass


@arguably.command
def data__show():
    pass


@arguably.command
def data__write():
    pass


@arguably.command
def config__load():
    pass


@arguably.command
def config__save():
    pass


@arguably.command
def config__set():
    pass


@arguably.command
def config__show():
    pass


@arguably.command(alias="auth")
def authenticate():
    """Authenticate with Deep Origin"""
    auth.get_tokens(refresh=False)


@arguably.command
def context():
    pass


def main():
    arguably.run()
