import asyncio
from contextlib import suppress
from datetime import datetime
from sre_parse import State
from typing import *

import pyppeteer
import pyppeteer.errors
from pyppeteer.element_handle import ElementHandle
from typer import *

from .common import *
from .ftdna import *


app = Typer()


def ftdna_signin(username: str, password: str, *, http_timeout: float):
	echo(f"Signing in to FTDNA with user `{username}`...")

	ftdna_signin_url = "https://www.familytreedna.com/sign-in"

	async def run() -> List[Dict]:
		browser = await pyppeteer.launch(headless = True, timeout = http_timeout * 1000)

		page = await browser.newPage()
		await page.goto(ftdna_signin_url)
		userInfoForm: ElementHandle = await page.waitForSelector("form[name = userInfoForm]")
		kitNumberInput: ElementHandle = await userInfoForm.querySelector("input[name = kitNumber]")
		passwordInput: ElementHandle = await userInfoForm.querySelector("input[name = password]")
		submitButton: ElementHandle = await userInfoForm.querySelector("button[type = submit]")

		await kitNumberInput.type(username)
		await passwordInput.type(password)

		wait_for_navigation = asyncio.ensure_future(page.waitForNavigation(waitUntil = "networkidle0"))
		wait_for_error_message = asyncio.ensure_future(page.waitForSelector("div#error-message"))
		(done, pending), _ = await asyncio.gather(
			asyncio.wait(
				[wait_for_navigation, wait_for_error_message],
				return_when = asyncio.FIRST_COMPLETED,
			),
			submitButton.click(),
		)

		for task in pending:
			task.cancel()
			with suppress(asyncio.CancelledError):
				await task

		for task in done:
			await task

		if wait_for_navigation in done:
			cookies = await page.cookies()
		else:
			cookies = None

		await browser.close()

		return cookies

	cookies = asyncio.get_event_loop().run_until_complete(run())

	if cookies:
		echo(f"Successfully signed in to FTDNA.")
	else:
		secho(f"Error logging in to FTDNA.", fg = colors.RED, err = True)

	return cookies


def ftdna_signout(cookies: List[Dict], *, http_timeout: float):
	echo(f"Signing out from FTDNA...")

	ftdna_signout_url = "https://www.familytreedna.com/sign-out"

	async def run():
		browser = await pyppeteer.launch(headless = True, timeout = http_timeout * 1000)

		page = await browser.newPage()

		if cookies is not None:
			await page.setCookie(*cookies)

		await page.goto(ftdna_signout_url, waitUntil = "networkidle0")

		await browser.close()

	asyncio.get_event_loop().run_until_complete(run())

	echo(f"Successfully signed out from FTDNA.")


@app.command()
def session() -> None:
	"""Show information about signed-in FTDNA session."""

	with open_ftdna_login_cache() as shelf:
		if "cookies" in shelf:
			username: str = shelf["username"]
			dt: datetime = shelf["datetime"]
			local_dt = utc_to_local(dt)
			echo(f"Signed in to FTDNA with user `{username}` at {local_dt:%Y-%m-%d %H:%M:%S}.")
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
