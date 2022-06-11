import importlib.metadata
import os
import shelve
from datetime import datetime, timezone
from enum import Enum, auto
from itertools import takewhile
from os import PathLike
from pathlib import Path
from typing import *

import pandas as pd
from bs4 import NavigableString, Tag
from mechanicalsoup import Form, StatefulBrowser
from platformdirs import PlatformDirs
from typer import *


T = TypeVar("T")
NDFrameT = TypeVar("NDFrameT", bound = pd.core.generic.NDFrame)

package = importlib.metadata.metadata(__package__)
platform_dirs = PlatformDirs(appname = package["Name"], version = package["Version"])

tree_path = "ycomp-yfull-tree.parquet"
snps_path = "ycomp-snps.parquet"
kits_snp_path = "ycomp-kits-snp.parquet"
kits_str_path = "ycomp-kits-str.parquet"


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


def first(seq: Sequence[T]) -> Optional[T]:
	return seq[0] if len(seq) > 0 else None


def utc_to_local(utc_dt: datetime) -> datetime:
	return utc_dt.replace(tzinfo = timezone.utc).astimezone(tz = None)


def float_not(a: NDFrameT) -> NDFrameT:
	return 1 - a


def float_and(a: NDFrameT, b: NDFrameT) -> NDFrameT:
	return a * b


def float_or(a: NDFrameT, b: NDFrameT) -> NDFrameT:
	return float_not(float_and(float_not(a.fillna(0)), float_not(b.fillna(0))))


def float_xor(a: NDFrameT, b: NDFrameT) -> NDFrameT:
	return float_and(float_or(a, b), float_not(float_and(a, b)))


def float_isna(a: NDFrameT) -> NDFrameT:
	return a.isna().astype("float16", copy = False)


def float_notna(a: NDFrameT) -> NDFrameT:
	return a.notna().astype("float16", copy = False)


def debug_mode() -> bool:
	return os.getenv("DEBUG") == "1"


def load_db(path: PathLike) -> pd.DataFrame:
	try:
		return pd.read_parquet(path)
	except (FileNotFoundError, pd.errors.EmptyDataError):
		return None


def store_db(path: PathLike, df: pd.DataFrame) -> None:
	df.to_parquet(path)


def merge_db(path: PathLike, new_df: pd.DataFrame) -> None:
	df = load_db(path)
	if df is None:
		df = new_df
	else:
		# Add new rows, removing old duplicates.
		df = pd.concat([df, new_df], axis = 0)
		df = df[~df.index.duplicated(keep = "last")]

	store_db(path, df)


def delete_db(path: PathLike) -> None:
	from pathlib import Path

	Path(path).unlink(missing_ok = True)


def get_haplogroup_lineage(tree_df: pd.DataFrame, haplogroup: str) -> List[str]:
	lineage = []
	cur_hg = haplogroup
	while pd.notna(cur_hg):
		lineage.append(cur_hg)
		try:
			cur_hg = tree_df.loc[cur_hg, "Parent Haplogroup"]
		except KeyError:
			break

	return lineage


def get_common_lineage(a_lineage: Sequence[str], b_lineage: Sequence[str]) -> Sequence[str]:
	return list(takewhile(lambda hgs: hgs[0] == hgs[1], zip(reversed(a_lineage), reversed(b_lineage))))


def include_kit(tree_df: pd.DataFrame, a_lineage: Sequence[str], b: pd.Series, *, haplogroup: Optional[str], haplogroup_max_diff: Optional[int]) -> bool:
	b_lineage = get_haplogroup_lineage(tree_df, b["Haplogroup"])
	common_lineage = get_common_lineage(a_lineage, b_lineage)
	a_lineage_pos = len(a_lineage) - len(common_lineage)
	b_lineage_pos = len(b_lineage) - len(common_lineage)

	if haplogroup is not None and haplogroup not in b_lineage:
		return False
	if haplogroup_max_diff is not None and abs(a_lineage_pos - b_lineage_pos) > haplogroup_max_diff:
		return False
	return True


def get_cache_dir() -> Path:
	path = Path(platform_dirs.user_cache_dir)
	path.mkdir(parents = True, exist_ok = True)
	return path


def open_ftdna_login_cache():
	return shelve.open(str(get_cache_dir() / "ftdna-login"))


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
