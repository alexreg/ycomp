import os
from enum import Enum, auto
from itertools import takewhile
from os import PathLike
from typing import *

import lxml
import lxml.etree
import lxml.html
import pandas as pd
import requests
from typer import *


T = TypeVar("T")
NDFrameT = TypeVar("NDFrameT", bound = pd.core.generic.NDFrame)

tree_path = "ycomp-yfull-tree.parquet"
snps_path = "ycomp-snps.parquet"
kits_snp_path = "ycomp-kits-snp.parquet"
kits_str_path = "ycomp-kits-str.parquet"


class DownloadFtdnaError(Exception, Enum):
	UNKNOWN_PAGE_LAYOUT = auto()
	NOT_GROUP_MEMBER = auto()

	def __str__(self) -> str:
		if self == self.UNKNOWN_PAGE_LAYOUT:
			return "Unknown page layout"
		elif self == self.NOT_GROUP_MEMBER:
			return "Not group member"


def _match_class(context, arg: str) -> bool:
	context_node = cast(lxml.html.HtmlElement, context.context_node)
	return arg in context_node.classes


def first(seq: Sequence[T]) -> Optional[T]:
	return seq[0] if len(seq) > 0 else None


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


def download_ftdna_kits(url: str, *, page_size: int = 500, http_timeout: float) -> pd.DataFrame:

	def id_to_form_input_name(id: str) -> str:
		return "ctl00$" + id.replace("_", "$")

	kits_df: pd.DataFrame = None

	prelim: bool = True
	pagination: bool = True
	page: int = 1
	data: Dict[str, str] = {}
	gridview_input_name: str
	page_size_input_name: str

	echo(f"Begun fetching kit data from <{url}> (page size = {page_size}).")

	from lxml.html import HtmlElement

	while True:
		if prelim:
			echo(f"Fetching page {page} (preliminary)...")
		else:
			echo(f"Fetching page {page}...")

		response = requests.post(url, data, timeout = http_timeout)

		if response.status_code != requests.codes.ok:
			echo(f"Page not found (HTTP {response.status_code}); finishing...")
			break

		html: HtmlElement = lxml.html.document_fromstring(response.text)

		if first(html.cssselect("div#MainContent_pnlHiddenYResults")) is not None:
			raise DownloadFtdnaError.NOT_GROUP_MEMBER

		form: HtmlElement = first(html.cssselect("form#form1"))
		if form is None:
			raise DownloadFtdnaError.UNKNOWN_PAGE_LAYOUT

		gridview_div: HtmlElement = first(form.cssselect("div.AspNet-GridView"))
		if gridview_div is None:
			raise DownloadFtdnaError.UNKNOWN_PAGE_LAYOUT

		table: HtmlElement = first(gridview_div.cssselect("table"))
		if table is None:
			raise DownloadFtdnaError.UNKNOWN_PAGE_LAYOUT

		if prelim:
			page_size_input: HtmlElement = first(form.cssselect("input[id *= 'tbPageSize']"))
			page_size_input_name = page_size_input.get("name")

			if int(page_size_input.get("value", 0)) == page_size:
				prelim = False
			else:
				# Prepare request to update page size.
				page_size_input.set("value", str(page_size))
				data = dict(form.fields)
				data["__EVENTTARGET"] = page_size_input_name
				data["__EVENTARGUMENT"] = ""

				continue

		# Check if pagination is active.
		pagination = len(form.cssselect("div.AspNet-GridView-Pagination")) > 0

		gridview_input_name = id_to_form_input_name(gridview_div.get("id"))

		echo(f"Processing page {page}...")

		page_df = first(pd.read_html(lxml.html.tostring(table)))

		if page > 1:
			# Drop header row of table.
			page_df.drop(page_df.index[0], inplace = True)

		if kits_df is None:
			kits_df = page_df
		else:
			kits_df = pd.concat([kits_df, page_df], axis = 0)

		if pagination:
			# Prepare request for next page.
			page += 1
			data = dict(form.fields)
			data["__EVENTTARGET"] = gridview_input_name
			data["__EVENTARGUMENT"] = f"Page${page}"
		else:
			break

	echo(f"Finished fetching kits from <{url}>.")

	return kits_df


ns = lxml.etree.FunctionNamespace("http://noldorin.com/xmlns/xpath/css")
ns.prefix = "css"
ns["class"] = _match_class
