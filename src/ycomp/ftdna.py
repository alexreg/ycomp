import shelve
from enum import Enum, auto
from typing import *

import pandas as pd
from bs4 import NavigableString, Tag
from mechanicalsoup import Form, StatefulBrowser
from platformdirs import PlatformDirs

from .common import *


class DownloadFtdnaError(Exception, Enum):
	NOT_FOUND = auto()
	NOT_GROUP_MEMBER = auto()
	UNKNOWN_PAGE_LAYOUT = auto()

	def __str__(self) -> str:
		if self == self.NOT_FOUND:
			return "Group not found"
		elif self == self.NOT_GROUP_MEMBER:
			return "Not group member or not signed in"
		elif self == self.UNKNOWN_PAGE_LAYOUT:
			return "Unknown page layout"
		else:
			raise ValueError(f"Invalid error kind")


def open_ftdna_login_cache():
	return shelve.open(str(cache_dir() / "ftdna-login"))


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
		if browser.page.select_one("div#MainContent_pnlHiddenYResults") is not None:
			raise DownloadFtdnaError.NOT_GROUP_MEMBER

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
				browser.submit_selected(
					data = {
						"__EVENTTARGET": page_size_input_name,
						"__EVENTARGUMENT": "",
					},
					timeout = http_timeout,
				)

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

		if page > 1:
			# Drop header row of table.
			page_df.drop(page_df.index[0], inplace = True)

		# Check if data frame is same as last.
		if prev_page_df is not None and page_df.equals(prev_page_df):
			break

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
		browser.submit_selected(
			data = {
				"__EVENTTARGET": gridview_input_name,
				"__EVENTARGUMENT": f"Page${next_page}",
			},
			timeout = http_timeout,
		)

	echo(f"Finished fetching kits.")

	return kits_df


def ftdna_normalize_columns(kits_df: pd.DataFrame):
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
