import re
import urllib.parse
from itertools import chain
from os import PathLike
from pathlib import Path
from typing import *

import numpy as np
import pandas as pd
from typer import *

from .common import *
from .db import *
from .ftdna import *


app = Typer()

ftdna_url_template = "https://www.familytreedna.com/public/{0}?iframe=ysnp"


def get_relevant_snps(tree_df: pd.DataFrame, max_age = int) -> Set[str]:
	filtered_df = tree_df[tree_df["TMRCA (ybp)"] <= max_age]
	primary_snps = chain.from_iterable(filtered_df["Primary SNPs"])
	extra_snps = chain.from_iterable(filtered_df["Extra SNPs"])
	return set(chain(primary_snps, extra_snps))


def get_irrelevant_snps(tree_df: pd.DataFrame, max_age = int) -> Set[str]:
	filtered_df = tree_df[tree_df["TMRCA (ybp)"] > max_age]
	primary_snps = chain.from_iterable(filtered_df["Primary SNPs"])
	extra_snps = chain.from_iterable(filtered_df["Extra SNPs"])
	return set(chain(primary_snps, extra_snps))


def get_yfull_df(data: Union[IO, PathLike]) -> Tuple[Optional[str], List[str], pd.DataFrame]:
	df: pd.DataFrame = pd.read_csv(data, sep = ";", index_col = 0, header = None, names = [0, 1, 2, 3, 4])
	info_df: pd.DataFrame = df.iloc[: 3, : 1]
	snps_df: pd.DataFrame = df.iloc[3 :]

	haplogroup = cast(str, info_df.loc["Haplogroup", 1])
	terminal_snps = cast(str, info_df.loc["Terminal SNPs", 1]).split(" | ")

	def fix_cols(read_count: Any, rating: Any) -> Tuple[Any, Any]:
		if isinstance(read_count, str) and cast(str, read_count).startswith("*"):
			rating = read_count
			read_count = np.nan
		return read_count, rating

	fixed_cols_df = pd.DataFrame((fix_cols(read_count, rating) for read_count, rating in zip(snps_df[2], snps_df[3])), index = snps_df.index, columns = [2, 3])
	snps_df = pd.concat([snps_df.loc[:, : 1], fixed_cols_df], axis = 1)

	def call_to_bool(call: Any) -> Any:
		if pd.isna(call):
			return np.nan

		if call == "no call":
			return np.nan
		elif call == "ambiguous":
			return np.nan
		elif call == "positive":
			return True
		elif call == "negative":
			return False
		elif call == "false positive":
			return np.nan
		elif call == "false negative":
			return np.nan
		else:
			raise ValueError(f"Invalid SNP call '{call}'.")

	def read_count_to_int(read_count: Any) -> Any:
		if pd.isna(read_count):
			return np.nan

		count, _, read_str = cast(str, read_count).partition(" ")
		if read_str == "read":
			return int(count)
		else:
			raise ValueError(f"Invalid SNP read count '{read_count}'.")

	def rating_to_int(rating: Any) -> Any:
		if pd.isna(rating):
			raise ValueError(f"Invalid SNP rating '{rating}'.")

		rating_s: str = cast(str, rating)
		if rating_s.strip("*") != "":
			raise ValueError(f"Invalid SNP rating '{rating_s}'.")

		return len(rating)

	snps_df[1] = snps_df[1].apply(call_to_bool).astype("boolean")
	snps_df[2] = snps_df[2].apply(read_count_to_int).astype("Int32")
	snps_df[3] = snps_df[3].apply(rating_to_int).astype("int8")

	new_rows: list[pd.Series] = []
	for multi_snp, row in snps_df.iterrows():
		for snp in (snp.strip() for snp in cast(str, multi_snp).split("/")):
			new_rows.append(row.to_frame(snp).T)

	snps_df = pd.concat(new_rows, axis = 0)

	snps_df.rename(
		{
			1: "Call",
			2: "Read Count",
			3: "Rating",
		},
		axis = 1,
		inplace = True,
	)
	snps_df.rename_axis("SNP", axis = 0, inplace = True)

	return haplogroup, terminal_snps, snps_df


@app.command()
def add_yfull(
	kit: Optional[str] = Option(None, "--kit", "-k", help = "The kit number."),
	group: Optional[str] = Option(None, "--group", help = "The group within which the sample clusters."),
	ancestor: Optional[str] = Option(None, "--ancestor", help = "The earliest known patrilineal ancestor."),
	country: Optional[str] = Option(None, "--country", help = "The country from which the ancestor came."),
	haplogroup: Optional[str] = Option(None, "--haplogroup", help = "The haplogroup of the sample."),
	file: Path = Argument(..., exists = True, dir_okay = False, help = "The YFull SNP file for the kit."),
) -> None:
	"""Add a YFull kit to the SNP database."""

	if kit is None:
		match = re.fullmatch(r"SNP_for_(YF\d+)_(\d+)", file.stem)
		if not match:
			raise BadParameter("Could not infer kit name from filename; specify it explicitly.", param = kit)
		kit = match.group(1)

	yfull_haplogroup, _yfull_terminal_snps, yfull_df = get_yfull_df(file)
	info_df = pd.DataFrame(
		{
			"Kit Number": pd.Series(kit, dtype = "str"),
			"Group": pd.Series(group, dtype = "str"),
			"Paternal Ancestor Name": pd.Series(ancestor, dtype = "str"),
			"Country": pd.Series(country, dtype = "str"),
			"Haplogroup": pd.Series(haplogroup or yfull_haplogroup, dtype = "str"),
		}
	)
	info_df.set_index("Kit Number", inplace = True)

	yfull_series = yfull_df["Call"]
	yfull_df = yfull_series.to_frame(kit).T
	yfull_df.rename_axis("Kit Number", axis = 0, inplace = True)

	kit_df = pd.concat([info_df, yfull_df], axis = 1)

	echo(f"Added kit {kit}.")

	merge_db(kits_snp_path, kit_df)
	echo(f"Kits SNP database written to `{kits_snp_path}`.")


@app.command()
def fetch_ftdna(
	ftdna_group: str = Option(..., "--group", "-g", help = "The name of the FTDNA group to fetch kits from."),
	page_size: int = Option(250, "--page-size", "-p", help = "The page size to use when fetching kits."),
) -> None:
	"""Fetch kit SNP data from FTDNA and store it in the database."""

	try:
		url = ftdna_url_template.format(urllib.parse.quote(ftdna_group))
		kits_df = ftdna_fetch_kits(url, page_size = page_size, http_timeout = 15 + 0.2 * page_size)
	except DownloadFtdnaError as e:
		if debug_mode():
			raise

		secho(f"ERROR: {e}", fg = colors.RED, err = True)
		raise Exit(1)

	echo(f"Processing kits from FTDNA...")

	# Clean data.
	kits_df.rename(columns = {"Short Hand": "Haplogroup"}, inplace = True)
	kits_df["Haplogroup"].replace(["-"], None, inplace = True)

	ftdna_normalize_columns(kits_df)

	kits_df.set_index("Kit Number", inplace = True)
	kits_df.index = kits_df.index.astype("str")
	if "Paternal Ancestor Name" in kits_df:
		kits_df["Paternal Ancestor Name"] = kits_df["Paternal Ancestor Name"].astype("str")
	else:
		kits_df["Paternal Ancestor Name"] = pd.Series(np.nan, index = kits_df.index, dtype = "str")
	kits_df["Haplogroup"] = kits_df["Haplogroup"].astype("str")

	def expand_row(row: Any) -> pd.Series:
		def get_snp_value(call: str) -> Optional[Tuple[str, bool]]:
			value: bool
			if call[-1] == "+":
				# Positive call
				value = True
			elif call[-1] == "-":
				# Negative call
				value = False
			elif call[-1] == "*":
				# No/ambiguous call
				return None
			else:
				raise ValueError(f"Invalid SNP call '{call}'.")

			return call[:-1], value

		cols = dict(filter(None, (get_snp_value(call) for call in cast(str, row).split(", ")))) if pd.notna(row) else None
		return pd.Series(cols, dtype = "boolean")

	# Expand out columns for SNPs.
	expanded_rows = kits_df["Confirmed SNPs"].apply(expand_row)

	kits_df = pd.concat(
		[
			pd.Series(np.nan, index = kits_df.index, dtype = "str", name = "Group"),
			kits_df["Paternal Ancestor Name"],
			pd.Series(np.nan, index = kits_df.index, dtype = "str", name = "Country"),
			kits_df["Haplogroup"],
			expanded_rows,
		],
		axis = 1,
	)

	echo(f"Finished processing kits from FTDNA.")

	merge_db(kits_snp_path, kits_df)
	echo(f"Kits SNP database written to `{kits_snp_path}`.")


@app.command()
def merge_str_db() -> None:
	"""Merge in kit information from the STR database."""

	from . import str_cmd

	kits_df = load_db(kits_snp_path)
	if kits_df is None:
		secho(f"ERROR: Kits SNP database does not exist.", fg = colors.RED, err = True)
		raise Exit(1)

	echo(f"Found {len(kits_df):,} kits in SNP database.")

	kits_str_df = load_db(str_cmd.kits_str_path)
	if kits_str_df is None:
		secho(f"ERROR: Kits STR database does not exist.", fg = colors.RED, err = True)
		raise Exit(1)

	echo(f"Found {len(kits_df):,} kits in STR database.")

	str_data_df = kits_str_df[["Group", "Paternal Ancestor Name", "Country", "Haplogroup"]]

	echo("Merging STR metadata into SNP metadata...")

	kits_df.update(str_data_df)

	echo("Finished merging.")

	merge_db(kits_snp_path, kits_df)
	echo(f"Kits SNP database written to `{kits_snp_path}`.")


@app.command()
def analyze(
	self_kit: str = Option(..., "--kit", "-k", help = "The kit to compare against."),
	snp_max_age: int = Option(3500, "--max-age", help = "The maximum age of SNPs to consider."),
	haplogroup: Optional[str] = Option(None, "--haplogroup", help = "The haplogroup clade, to filter by."),
	haplogroup_max_diff: Optional[int] = Option(None, "--haplogroup-max-diff", help = "The maximum difference between generations in the haplogroup tree, to filyer by."),
	output_file: Path = Option("ycomp-analysis-snp.csv", "--output", "-o", dir_okay = False, help = "The (CSV) file to write the analysis to."),
) -> None:
	"""Compare matches in the SNP database."""

	kits_df = load_db(kits_snp_path)
	if kits_df is None:
		secho(f"ERROR: Kits SNP database does not exist.", fg = colors.RED, err = True)
		raise Exit(1)

	echo(f"Found {len(kits_df):,} kits in SNP database.")

	snps_df = load_db(snps_path)
	if snps_df is None:
		secho(f"WARNING: SNPs database does not exist.", fg = colors.YELLOW, err = True)

	if snps_df is not None:
		echo(f"Found {len(snps_df):,} entries in SNPs database.")

	tree_df = load_db(tree_path)
	if tree_df is None:
		secho(f"WARNING: Tree database does not exist.", fg = colors.YELLOW, err = True)

	self_kit_hg = kits_df.loc[self_kit, "Haplogroup"]
	self_kit_lineage = get_haplogroup_lineage(tree_df, self_kit_hg)

	kits_df = kits_df[kits_df.apply(lambda kit: kit_matches_lineage(tree_df, self_kit_lineage, kit, haplogroup = haplogroup, haplogroup_max_diff = haplogroup_max_diff), axis = 1)]
	echo(f"Will compare {len(kits_df.index):,} kits.")

	if tree_df is not None:
		relevant_snps = get_relevant_snps(tree_df, max_age = snp_max_age)
		echo(f"Found {len(relevant_snps):,} relevant SNPs in tree.")
		# irrelevant_snps = get_irrelevant_snps(tree_df, max_age = max_age)
		# echo(f"Found {len(irrelevant_snps):,} irrelevant SNPs in tree.")
	else:
		irrelevant_snps = None # all SNPs

	kits_snps_df = kits_df.iloc[:, 4 :].astype("float16")

	# Drop columns for irrelevant SNPs.
	irrelevant_snps = set(kits_snps_df.columns) - relevant_snps
	kits_snps_df.drop(irrelevant_snps & set(kits_snps_df.columns), axis = 1, inplace = True)

	if snps_df is not None:
		# Merge columns for equivalent SNPs.
		new_cols: OrderedDict[str, str] = OrderedDict()

		def get_standard_snp_name(snp: str) -> str:
			try:
				return snps_df.loc[snp, "Standard Name"]
			except KeyError:
				return snp

		for col in kits_snps_df.columns:
			new_col = get_standard_snp_name(col)
			if new_col in new_cols:
				new_cols[new_col] = float_or(new_cols[new_col], kits_snps_df[col])
			else:
				new_cols[new_col] = kits_snps_df[col]

		kits_snps_df = pd.concat(new_cols, axis = 1)

	echo(f"Will compare {len(kits_snps_df.columns):,} SNPs.")

	try:
		self_kit_s: pd.Series = kits_snps_df.loc[self_kit]
		kits_df.drop(self_kit, axis = 0, inplace = True)
		kits_snps_df = kits_snps_df.drop(self_kit, axis = 0)
	except KeyError:
		secho(f"ERROR: Kit {self_kit} not found.", fg = colors.RED, err = True)
		raise Exit(1)

	echo(f"Starting analysis...")

	def get_snp_list(kit_s: pd.Series) -> str:
		return ", ".join(kit_s[kit_s == 1].index)

	self_kit_isna_s = float_isna(self_kit_s)

	def get_shared_snps(kit_s: pd.Series) -> pd.Series:
		return float_and(self_kit_s, kit_s)

	def get_assumed_shared_snps(kit_s: pd.Series) -> pd.Series:
		return float_or(float_and(self_kit_s, float_isna(kit_s)), float_and(self_kit_isna_s, kit_s))

	shared_snps = kits_snps_df.apply(lambda s: get_shared_snps(s), axis = 1)
	list_shared_snps = shared_snps.apply(get_snp_list, axis = 1)
	num_shared_snps = shared_snps.sum(axis = 1)
	num_comp_shared_snps = shared_snps.count(axis = 1)
	frac_shared_snps = num_shared_snps / num_comp_shared_snps

	assum_shared_snps = kits_snps_df.apply(lambda s: get_assumed_shared_snps(s), axis = 1)
	list_assum_shared_snps = assum_shared_snps.apply(get_snp_list, axis = 1)
	num_assum_shared_snps = assum_shared_snps.sum(axis = 1)
	num_comp_assum_shared_snps = assum_shared_snps.count(axis = 1)
	frac_assum_shared_snps = num_assum_shared_snps / num_comp_assum_shared_snps

	shared_snps = float_or(shared_snps, assum_shared_snps)
	list_all_shared_snps = shared_snps.apply(get_snp_list, axis = 1)
	num_all_shared_snps = num_shared_snps + num_assum_shared_snps
	num_comp_all_shared_snps = num_comp_shared_snps + num_comp_assum_shared_snps
	frac_all_shared_snps = num_all_shared_snps / num_comp_all_shared_snps

	results_df = pd.concat(
		[
			kits_df["Group"],
			kits_df["Paternal Ancestor Name"],
			kits_df["Country"],
			kits_df["Haplogroup"],
			list_shared_snps.rename("Shared SNPs", inplace = True),
			num_shared_snps.astype("int32").rename("Shared SNPs (#)", inplace = True),
			num_comp_shared_snps.astype("int32").rename("Shared SNPs (# compared)", inplace = True),
			frac_shared_snps.astype("float32").rename("Shared SNPs (fraction)", inplace = True),
			list_assum_shared_snps.rename("Assumed Shared SNPs", inplace = True),
			num_assum_shared_snps.astype("int32").rename("Assumed Shared SNPs (#)", inplace = True),
			num_comp_assum_shared_snps.astype("int32").rename("Assumed Shared SNPs (# compared)", inplace = True),
			frac_assum_shared_snps.astype("float32").rename("Assumed Shared SNPs (fraction)", inplace = True),
			list_all_shared_snps.rename("All Shared SNPs", inplace = True),
			num_all_shared_snps.astype("int32").rename("All Shared SNPs (#)", inplace = True),
			num_comp_all_shared_snps.astype("int32").rename("All Shared SNPs (# compared)", inplace = True),
			frac_all_shared_snps.astype("float32").rename("All Shared SNPs (fraction)", inplace = True),
		],
		axis = 1,
	)

	echo(f"Finished analysis.")

	results_df.sort_values("All Shared SNPs (#)", ascending = False).to_csv(output_file)
	echo(f"Analysis written to `{output_file}`.")
