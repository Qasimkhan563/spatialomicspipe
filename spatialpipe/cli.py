import click
from spatialpipe.pipeline import run_ingest_pipeline


@click.group()
def cli():
    pass


@cli.command()
@click.argument("config")
def ingest(config):
    """Run ingestion pipeline"""
    run_ingest_pipeline(config)
    click.echo("Ingestion complete.")


if __name__ == "__main__":
    cli()
