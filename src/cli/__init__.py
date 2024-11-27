import typer
from deeporigin.data_hub import api as data_hub_api

app = typer.Typer()
app.add_typer(data_hub_api.app, name="data")


def main():
    app()
