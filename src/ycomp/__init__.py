from typing import *

import pandas as pd
from click.exceptions import ClickException
from typer import *

from . import common, snp_cmd, str_cmd, tree_cmd


app = Typer(
	context_settings = {
		"help_option_names": ["-h", "--help"],
	},
)
app.add_typer(tree_cmd.app)
app.add_typer(snp_cmd.app)
app.add_typer(str_cmd.app)


pd.options.mode.chained_assignment = None


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
	import sys

	import typer

	try:
		command = typer.main.get_command(app)
		result = command(standalone_mode = False)
		return result
	except Abort:
		secho(f"aborted", fg = colors.RED, err = True)
		return 1
	except ClickException as e:
		secho(e, fg = colors.RED, err = True)
		return 1
	except Exception as e:
		if common.debug_mode():
			raise
		else:
			secho(f"UNEXPECTED ERROR: {e}", fg = colors.RED, err = True)
			return 1
