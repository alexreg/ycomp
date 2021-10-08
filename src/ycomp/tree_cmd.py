import re
import urllib.parse
from decimal import Decimal
from typing import *

import lxml
import lxml.html
import numpy as np
import pandas as pd
import requests
from typer import *

from .common import *


app = Typer()

yfull_tree_url_template = "https://www.yfull.com/tree/{0}/"

snp_name_pattern = re.compile(r"(?P<name>[A-Z0-9+.=]+)(?P<info>\([A-Z]+\))?")
age_pattern = re.compile(r"formed (?P<age>\d+) ybp, TMRCA (?P<tmrca>\d+) ybp")
age_bounds_pattern = re.compile(r"formed CI (?P<age_cl>\d+)% (?P<age_min>\d+)<->(?P<age_max>\d+) ybp, TMRCA CI (?P<tmrca_cl>\d+)% (?P<tmrca_min>\d+)<->(?P<tmrca_max>\d+) ybp")


@app.command()
def download_yfull(
	haplogroup: str = Option(..., "--haplogroup", "--hg", help = "The haplogroup of the subtree to download."),
) -> None:
	"""Download a subtree of the Y-haplogroup tree on YFull, and store it."""

	try:
		echo(f"Begun downloading YFull tree.")

		url = yfull_tree_url_template.format(urllib.parse.quote(haplogroup))
		response = requests.get(url, timeout = 60)

		if response.status_code == requests.codes.not_found:
			secho(f"Haplogroup {haplogroup} not found in YFull tree.", fg = colors.RED, err = True)
			raise Exit(1)

		response.raise_for_status()

		echo(f"Finished downloading YFull tree.")

		if response.status_code != requests.codes.ok:
			secho(f"Page not found (HTTP {response.status_code}); finishing...", fg = colors.RED, err = True)
			raise Exit(1)
	except (requests.ConnectionError, requests.Timeout) as e:
		secho(f"ERROR: HTTP request failed: {e}", fg = colors.RED, err = True)
		raise Exit(1)

	echo(f"Begun processing YFull tree.")

	from lxml.html import HtmlElement

	html: HtmlElement = lxml.html.document_fromstring(response.text)
	tree_ul: HtmlElement = html.xpath("//ul[@id = 'tree']")[0]

	found_snps: list[Tuple[str, str]] = []

	def snps_to_list(snps: str) -> List[str]:
		def clean_snp(snp: str) -> str:
			match = snp_name_pattern.fullmatch(snp.strip())

			if match:
				return match.group("name")
			else:
				raise ValueError(f"Invalid SNP '{snp}'.")

		if not snps:
			return []

		all_snps_list = []
		multi_snps_list = [snp for snp in snps.split(" * ")]
		for multi_snp in multi_snps_list:
			snps_list = [clean_snp(snp) for snp in multi_snp.split("/")]
			all_snps_list.extend(snps_list)

			for snp in snps_list:
				found_snps.append((snp, multi_snp))

		return all_snps_list

	def traverse_tree(tree_ul: HtmlElement, parent_haplogroup: Optional[str] = None) -> Iterable[pd.Series]:
		item_li_list: list[HtmlElement] = tree_ul.xpath("li")
		for item_li in item_li_list:
			haplogroup_a: HtmlElement = first(item_li.xpath("a[1]"))
			if haplogroup_a is None:
				continue

			snp_span: HtmlElement = first(item_li.xpath("span[css:class('yf-snpforhg')]"))
			plus_snp_span: HtmlElement = first(item_li.xpath("span[css:class('yf-plus-snps')]"))
			age_span: HtmlElement = first(item_li.xpath("span[css:class('yf-age')]"))
			inner_ul: HtmlElement = first(item_li.xpath("ul[1]"))

			haplogroup = haplogroup_a.text_content()

			primary_snps = snps_to_list(snp_span.text_content())

			if plus_snp_span is not None:
				extra_snps = snps_to_list(plus_snp_span.attrib["title"])
			else:
				extra_snps = []

			age = None
			age_cl = None
			age_min = None
			age_max = None
			tmrca = None
			tmrca_cl = None
			tmrca_min = None
			tmrca_max = None
			if age_span is not None:
				age_text = age_span.text_content()
				age_match = age_pattern.fullmatch(age_text)
				if age_match:
					age = int(age_match.group("age"))
					tmrca = int(age_match.group("tmrca"))
				else:
					secho(f"WARNING: unexpected format of text for haplogroup age: '{age_text}'", fg = colors.YELLOW, err = True)

				age_bounds_text = age_span.attrib["title"]
				age_bounds_match = age_bounds_pattern.fullmatch(age_bounds_text)
				if age_bounds_match:
					age_cl = Decimal(int(age_bounds_match.group("age_cl"))) / 100
					age_min = int(age_bounds_match.group("age_min"))
					age_max = int(age_bounds_match.group("age_max"))
					tmrca_cl = Decimal(int(age_bounds_match.group("tmrca_cl"))) / 100
					tmrca_min = int(age_bounds_match.group("tmrca_min"))
					tmrca_max = int(age_bounds_match.group("tmrca_max"))
				else:
					secho(f"WARNING: unexpected format of text for haplogroup age bounds: '{age_bounds_text}'", fg = colors.YELLOW, err = True)

			yield [
				haplogroup,
				parent_haplogroup,
				primary_snps,
				extra_snps,
				age,
				age_cl,
				age_min,
				age_max,
				tmrca,
				tmrca_cl,
				tmrca_min,
				tmrca_max,
			]

			if inner_ul is not None:
				yield from traverse_tree(inner_ul, haplogroup)

	haplogroups = traverse_tree(tree_ul)
	tree_df = pd.DataFrame(
		haplogroups,
		columns = [
			"Haplogroup",
			"Parent Haplogroup",
			"Primary SNPs",
			"Extra SNPs",
			"Age (ybp)",
			"Age (CL%)",
			"Age (min, ybp)",
			"Age (max, ybp)",
			"TMRCA (ybp)",
			"TMRCA (CL%)",
			"TMRCA (min, ybp)",
			"TMRCA (max, ybp)",
		],
	)
	tree_df.set_index("Haplogroup", inplace = True)

	tree_df["Age (ybp)"] = tree_df["Age (ybp)"].astype("Int32")
	tree_df["Age (min, ybp)"] = tree_df["Age (min, ybp)"].astype("Int32")
	tree_df["Age (max, ybp)"] = tree_df["Age (max, ybp)"].astype("Int32")
	tree_df["TMRCA (ybp)"] = tree_df["TMRCA (ybp)"].astype("Int32")
	tree_df["TMRCA (min, ybp)"] = tree_df["TMRCA (min, ybp)"].astype("Int32")
	tree_df["TMRCA (max, ybp)"] = tree_df["TMRCA (max, ybp)"].astype("Int32")

	snps_df = pd.DataFrame(
		found_snps,
		columns = [
			"SNP",
			"Standard Name",
		],
	)
	snps_df.set_index("SNP", inplace = True)

	echo(f"Finished processing YFull tree.")

	merge_db(tree_path, tree_df)
	echo(f"Tree database written to `{tree_path}`.")

	merge_db(snps_path, snps_df)
	echo(f"SNPs database written to `{snps_path}`.")


@app.command()
def prune_tree(
	haplogroups: str = Argument(..., help = "The reg-exp that matches haplogroups to be kept."),
) -> None:
	"""Prune the local tree database."""

	tree_df = load_db(tree_path)
	if tree_df is None:
		secho(f"ERROR: Tree database does not exist.", fg = colors.RED, err = True)
		raise Exit(1)

	echo(f"Found {len(tree_df):,} haplogroups in tree database.")

	haplogroups_pattern = re.compile(haplogroups, re.RegexFlag.IGNORECASE)
	pruned_tree_df = tree_df.filter(regex = haplogroups_pattern, axis = 0)

	retained_fraction = Decimal(len(pruned_tree_df)) / len(tree_df)
	echo(f"Pruned tree has {len(tree_df):,} haplogroups ({retained_fraction:.1%} retained).")

	if not confirm("Are you sure you want to prune the tree?"):
		return

	store_db(tree_path, pruned_tree_df)
	echo(f"Tree database written to `{tree_path}`.")


@app.command()
def delete_tree() -> None:
	"""Delete the local tree database."""

	delete_db(tree_path)
	echo(f"Tree database deleted from `{tree_path}`.")
