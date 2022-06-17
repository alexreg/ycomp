import asyncio
import shelve
from contextlib import suppress
from enum import Enum, auto
from typing import *

import pandas as pd
import pyppeteer
import pyppeteer.errors
from bs4 import NavigableString, Tag
from mechanicalsoup import Form, StatefulBrowser
from platformdirs import PlatformDirs
from pyppeteer.element_handle import ElementHandle

from .common import *


Cookies = List[Dict[str, Union[str, int, bool]]]


class DownloadFtdnaError(Exception, Enum):
	NOT_FOUND = auto()
	RESULTS_UNAVAILABLE = auto()
	RESULTS_HIDDEN = auto()
	UNKNOWN_PAGE_LAYOUT = auto()

	def __str__(self) -> str:
		if self == self.NOT_FOUND:
			return "Group not found"
		elif self == self.RESULTS_UNAVAILABLE:
			return "Results page unavailable"
		elif self == self.RESULTS_HIDDEN:
			return "Results page hidden"
		elif self == self.UNKNOWN_PAGE_LAYOUT:
			return "Unknown page layout"
		else:
			raise ValueError(f"Invalid error kind")


def open_ftdna_login_cache() -> shelve.Shelf:
	return shelve.open(str(cache_dir() / "ftdna-login"))


def ftdna_normalize_columns(kits_df: pd.DataFrame) -> None:
	if "Last Name" in kits_df.columns:
		if "Paternal Ancestor Name" not in kits_df.columns:
			kits_df.rename(columns = {"Last Name": "Paternal Ancestor Name"}, inplace = True)
		else:
			kits_df["Paternal Ancestor Name"] = kits_df["Paternal Ancestor Name"].fillna(kits_df["Last Name"])
			kits_df.drop("Last Name", axis = 1, inplace = True)

	if "Name" in kits_df.columns:
		if "Paternal Ancestor Name" not in kits_df.columns:
			kits_df.rename(columns = {"Name": "Paternal Ancestor Name"}, inplace = True)
		else:
			kits_df["Paternal Ancestor Name"] = kits_df["Paternal Ancestor Name"].fillna(kits_df["Name"])
			kits_df.drop("Name", axis = 1, inplace = True)


def ftdna_fetch_kits(url: str, *, page_size: Optional[int] = None, http_timeout: Optional[float] = None) -> pd.DataFrame:
	def id_to_form_input_name(id: str) -> str:
		return "ctl00$" + id.replace("_", "$")

	kits_df: pd.DataFrame = None

	prelim: bool = True
	gridview_input_name: str
	page_size_input_name: str
	page_df = None

	browser = StatefulBrowser()

	# Configure cookies for requests.
	cookie_jar = browser.get_cookiejar()
	with open_ftdna_login_cache() as shelf:
		cookies = shelf.get("cookies", ())
		for cookie in cookies:
			cookie_jar.set(cookie["name"], cookie["value"], domain = cookie.get("domain"), path = cookie.get("path"))

	echo(f"Fetching kits from <{url}>...")
	response = browser.open(url, timeout = http_timeout)
	if not response.ok or response.url == "https://www.familytreedna.com/":
		raise DownloadFtdnaError.NOT_FOUND

	while True:
		if browser.page.select_one("div#MainContent_pnlNoYResults") is not None:
			raise DownloadFtdnaError.RESULTS_UNAVAILABLE

		if browser.page.select_one("div#MainContent_pnlHiddenYResults") is not None:
			raise DownloadFtdnaError.RESULTS_HIDDEN

		form: Form = browser.select_form("form#form1")
		if form is None:
			raise DownloadFtdnaError.UNKNOWN_PAGE_LAYOUT

		form_tag: Tag = form.form

		gridview_div: Tag = form_tag.select_one("div.AspNet-GridView")
		if gridview_div is None:
			raise DownloadFtdnaError.UNKNOWN_PAGE_LAYOUT

		table = gridview_div.select_one("table")
		if table is None:
			raise DownloadFtdnaError.UNKNOWN_PAGE_LAYOUT

		if prelim:
			page_size_input: Tag = form_tag.select_one("input[id *= 'tbPageSize']")
			page_size_input_name = page_size_input.get("name")

			if page_size is None or int(page_size_input.get("value", 0)) == page_size:
				prelim = False
			else:
				# Submit request to update page size.
				browser[page_size_input_name] = str(page_size)
				echo(f"Updating page size to {page_size}...")
				browser.form.choose_submit_none()
				response = browser.submit_selected(
					data = {
						"__EVENTTARGET": page_size_input_name,
						"__EVENTARGUMENT": "",
					},
					timeout = http_timeout,
				)
				response.raise_for_status()

				continue

		# Extract current and maximum page numbers.
		page = 1
		max_page = 1
		pagination_div: Tag = form_tag.select_one("div.AspNet-GridView-Pagination")
		if pagination_div:
			for child in pagination_div.children:
				if isinstance(child, Tag):
					if child.name.casefold() == "span":
						page = int(child.text.strip())
				elif isinstance(child, NavigableString):
					of_prefix = " of "
					if child.text.startswith(of_prefix):
						max_page = int(child.removeprefix(of_prefix).strip())

		gridview_input_name = id_to_form_input_name(gridview_div.get("id"))

		echo(f"Processing page {page} of {max_page}...")

		prev_page_df = page_df
		page_df = first(pd.read_html(str(table)))
		if page_df is None:
			raise DownloadFtdnaError.UNKNOWN_PAGE_LAYOUT

		# Check if data frame is same as last.
		if prev_page_df is not None and page_df.equals(prev_page_df):
			break

		# Drop header row of table for all but first page.
		if page > 1:
			page_df.drop(page_df.index[0], inplace = True)

		if kits_df is None:
			kits_df = page_df
		else:
			kits_df = pd.concat([kits_df, page_df], axis = 0)

		# Check if there are any more pages remaining to fetch.
		if page == max_page:
			break

		# Submit request for next page.
		next_page = page + 1
		echo(f"Fetching page {next_page}...")
		browser.form.choose_submit_none()
		response = browser.submit_selected(
			data = {
				"__EVENTTARGET": gridview_input_name,
				"__EVENTARGUMENT": f"Page${next_page}",
			},
			timeout = http_timeout,
		)
		response.raise_for_status()

	echo(f"Finished fetching kits.")

	return kits_df


def ftdna_signin(username: str, password: str, *, http_timeout: float) -> Optional[Cookies]:
	echo(f"Signing in to FTDNA with user `{username}`...")

	ftdna_signin_url = "https://www.familytreedna.com/sign-in"

	async def run() -> Optional[Cookies]:
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


def ftdna_signout(cookies: Optional[Cookies], *, http_timeout: float) -> None:
	echo(f"Signing out from FTDNA...")

	ftdna_signout_url = "https://www.familytreedna.com/sign-out"

	async def run() -> None:
		browser = await pyppeteer.launch(headless = True, timeout = http_timeout * 1000)

		page = await browser.newPage()

		if cookies is not None:
			await page.setCookie(*cookies)

		await page.goto(ftdna_signout_url, waitUntil = "networkidle0")

		await browser.close()

	asyncio.get_event_loop().run_until_complete(run())

	echo(f"Successfully signed out from FTDNA.")
