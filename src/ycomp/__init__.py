from typing import *

import pandas as pd
from click.exceptions import ClickException
from typer_cloup import *

from . import common, ftdna_cmd, snp_cmd, str_cmd, tree_cmd


app = Typer(
    context_settings={
        "help_option_names": ["-h", "--help"],
    },
)
app.add_sub(ftdna_cmd.app)
app.add_sub(tree_cmd.app)
app.add_sub(snp_cmd.app)
app.add_sub(str_cmd.app)

pd.options.mode.chained_assignment = None


@ftdna_cmd.app.callback()
def ftdna() -> None:
    """Work with FTDNA accounts."""

    pass


@tree_cmd.app.callback()
def tree() -> None:
    """Work with haplogroup tree data."""

    pass


@snp_cmd.app.callback()
def snp() -> None:
    """Work with SNP data."""

    pass


@str_cmd.app.callback("str")
def str_() -> None:
    """Work with STR data."""

    pass


def main() -> int:
    import logging
    import os

    import typer_cloup as typer

    try:
        logging.basicConfig(level=os.environ.get("LOGLEVEL", "WARNING").upper())

        command = typer.main.get_command(app)
        result = command(standalone_mode=False)
        return result
    except Abort:
        secho(f"aborted", fg=colors.RED, err=True)
        return 1
    except ClickException as e:
        secho(e, fg=colors.RED, err=True)
        return 1
    except Exception as e:
        if common.debug_mode():
            raise
        else:
            secho(f"UNEXPECTED ERROR: {e}", fg=colors.RED, err=True)
            return 1
