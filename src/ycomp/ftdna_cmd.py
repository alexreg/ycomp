from datetime import datetime
from sre_parse import State
from typing import *

from typer import *

from .common import *
from .ftdna import *


app = Typer()


@app.command()
def refresh() -> None:
	"""Refresh the signed-in FTDNA session."""

	with open_ftdna_login_cache() as shelf:
		if "cookies" in shelf:
			cookies = ftdna_refresh(shelf.get("cookies"), http_timeout = 10)
			if cookies is not None:
				shelf["cookies"] = cookies


@app.command()
def session() -> None:
	"""Show information about signed-in FTDNA session."""

	with open_ftdna_login_cache() as shelf:
		if "cookies" in shelf:
			username: str = shelf["username"]
			dt: datetime = shelf["datetime"]

			expiries = (cookie.get("expires", -1) for cookie in shelf["cookies"])
			min_expiry = min(expiry for expiry in expiries if expiry != -1)
			min_expiry_dt = datetime.fromtimestamp(min_expiry)

			echo(f"Signed in to FTDNA with user `{username}` at {utc_to_local(dt):%Y-%m-%d %H:%M:%S}.")
			if datetime.utcnow() < min_expiry_dt:
				secho(f"Session expires at {utc_to_local(min_expiry_dt):%Y-%m-%d %H:%M:%S}.", fg = colors.GREEN)
			else:
				secho(f"Session expired.", fg = colors.YELLOW)
		else:
			echo(f"Not signed in to FTDNA.")


@app.command()
def signin(
	username: str = Argument(..., help = "The username of the FTDNA account"),
	password: str = Argument(..., help = "The password of the FTDNA account"),
) -> None:
	"""Sign in to FTDNA."""

	with open_ftdna_login_cache() as shelf:
		cookies = ftdna_signin(username, password, http_timeout = 10)
		if cookies is not None:
			shelf["username"] = username
			shelf["datetime"] = datetime.utcnow()
			shelf["cookies"] = cookies


@app.command()
def signout() -> None:
	"""Sign out from FTDNA."""

	with open_ftdna_login_cache() as shelf:
		ftdna_signout(shelf.get("cookies"), http_timeout = 10)

		if "username" in shelf:
			del shelf["username"]

		if "datetime" in shelf:
			del shelf["datetime"]

		if "cookies" in shelf:
			del shelf["cookies"]
