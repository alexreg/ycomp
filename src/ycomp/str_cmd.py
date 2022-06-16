import re
import urllib.parse
from os import PathLike
from pathlib import Path
from typing import *

import numpy as np
import pandas as pd
from more_itertools import distinct_combinations, flatten
from typer import *

from .common import *
from .db import *
from .ftdna import *


app = Typer()

ftdna_url_template = "https://www.familytreedna.com/public/{0}?iframe=yresults"


def get_yfull_df(data: Union[IO, PathLike]) -> pd.DataFrame:
	df: pd.DataFrame = pd.read_csv(data, sep = ";", index_col = 0, header = None)

	# Remove rows specific to YFull and FTDNA.
	df.drop((locus for locus in df.index if "_" in locus), axis = 0, inplace = True)

	def confidence_to_bool(confidence: str) -> bool:
		if confidence == "?":
			# Low confidence
			return True
		elif pd.isna(confidence):
			return False
		else:
			raise ValueError(f"Invalid SNP call confidence '{confidence}'.")

	# Split combined loci values into two columns.
	value_split_df: pd.DataFrame = df[1].str.split(".", 1, expand = True)
	df[1] = value_split_df[0].replace(["?", "n/a"], pd.NA).astype("Int32")
	df[2] = df[2].apply(lambda confidence: confidence_to_bool(confidence)).astype("bool")
	df[3] = value_split_df[1].fillna(value = np.nan)
	df = df.reindex([1, 3, 2], axis = 1)

	df.rename(
		{
			1: "Repeats",
			2: "Low Confidence",
			3: "Variant",
		},
		axis = 1,
		inplace = True,
	)
	df.rename_axis("Locus", axis = 0, inplace = True)

	return df


@app.command()
def add_yfull(
	kit: Optional[str] = Option(None, "--kit", "-k", help = "The kit number."),
	group: Optional[str] = Option(None, "--group", help = "The group within which the sample clusters."),
	ancestor: Optional[str] = Option(None, "--ancestor", help = "The earliest known patrilineal ancestor."),
	country: Optional[str] = Option(None, "--country", help = "The country from which the ancestor came."),
	haplogroup: Optional[str] = Option(None, "--haplogroup", help = "The haplogroup of the sample."),
	file: Path = Argument(..., exists = True, dir_okay = False, help = "The YFull STR file for the kit."),
) -> None:
	"""Add a YFull kit to the STR database."""

	if kit is None:
		match = re.fullmatch(r"STR_for_(YF\d+)_(\d+)", file.stem)
		if not match:
			raise BadParameter("Could not infer kit name from filename; specify it explicitly.", param = kit)
		kit = match.group(1)

	yfull_df = get_yfull_df(file)
	info_df = pd.DataFrame(
		{
			"Kit Number": pd.Series(kit, dtype = "str"),
			"Group": pd.Series(group, dtype = "str"),
			"Paternal Ancestor Name": pd.Series(ancestor, dtype = "str"),
			"Country": pd.Series(country, dtype = "str"),
			"Haplogroup": pd.Series(haplogroup, dtype = "str"),
		}
	)
	info_df.set_index("Kit Number", inplace = True)

	yfull_series = yfull_df["Repeats"]

	# Join together multi-copy loci values.
	new_vals: OrderedDict[str, Optional[list]] = OrderedDict()
	for index in yfull_series.index:
		value = yfull_series[index]
		locus, _, num = index.partition(".")

		if pd.isna(value):
			new_vals[locus] = None
			continue

		if num and locus in new_vals:
			new_col = cast(list, new_vals[locus])
			new_col.append(value)
			new_vals[locus] = new_col
		else:
			new_vals[locus] = [value]

	yfull_series = pd.Series(new_vals).apply(lambda lst: [int(x) for x in lst] if isinstance(lst, list) else None)
	yfull_df = yfull_series.to_frame(kit).T
	yfull_df.rename_axis("Kit Number", axis = 0, inplace = True)

	kit_df = pd.concat([info_df, yfull_df], axis = 1)

	echo(f"Added kit {kit}.")

	merge_db(kits_str_path, kit_df)
	echo(f"Kits STR database written to `{kits_str_path}`.")


@app.command()
def fetch_ftdna(
	ftdna_group: str = Option(..., "--group", "-g", help = "The name of the FTDNA group to fetch kits from."),
	page_size: int = Option(500, "--page-size", "-p", help = "The page size to use when fetching kits."),
) -> None:
	"""Fetch kit STR data from FTDNA and store it in the database."""

	try:
		url = ftdna_url_template.format(urllib.parse.quote(ftdna_group))
		kits_df = ftdna_fetch_kits(url, page_size = page_size, http_timeout = 15 + 0.05 * page_size)
	except DownloadFtdnaError as e:
		if debug_mode():
			raise

		secho(f"ERROR: {e}", fg = colors.RED, err = True)
		raise Exit(1)

	echo(f"Processing kits from FTDNA...")

	# Clean data.
	kits_df["Haplogroup"].replace(["-"], None, inplace = True)

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

	kits_df.set_index("Kit Number", inplace = True)
	kits_df.index = kits_df.index.astype("str")
	kits_df["Paternal Ancestor Name"] = kits_df["Paternal Ancestor Name"].astype("str")
	kits_df["Country"] = kits_df["Country"].astype("str")
	kits_df["Haplogroup"] = kits_df["Haplogroup"].astype("str")

	group = None
	def get_group(kit_s: pd.Series) -> pd.Series:
		nonlocal group

		if kit_s[0] == kit_s[-1]:
			group = kit_s[0]
			return None
		else:
			return group

	# Add group column and assign each kit to corresponding group.
	group_series = kits_df.apply(get_group, axis = 1)
	kits_df.insert(0, "Group", group_series)
	kits_df.drop(kits_df[kits_df["Group"].isna()].index, inplace = True)

	# Parse multi-copy-loci columns.
	new_cols = []
	for col in kits_df.columns[4 :]:
		series = kits_df[col].str.split("-") \
			.apply(lambda lst: [int(x) for x in lst] if isinstance(lst, list) else None)
		series.name = cast(str, series.name).upper()
		new_cols.append(series)

	kits_df = pd.concat(
		[
			kits_df["Group"],
			kits_df["Paternal Ancestor Name"],
			kits_df["Country"],
			kits_df["Haplogroup"],
		] + new_cols,
		axis = 1,
	)

	echo(f"Finished processing kits from FTDNA.")

	merge_db(kits_str_path, kits_df)
	echo(f"Kits STR database written to `{kits_str_path}`.")


@app.command()
def analyze(
	self_kit: str = Option(..., "--kit", "-k", help = "The kit to compare against."),
	haplogroup: Optional[str] = Option(None, "--haplogroup", help = "The haplogroup clade, to filter by."),
	haplogroup_max_diff: Optional[int] = Option(None, "--haplogroup-max-diff", help = "The maximum difference between generations in the haplogroup tree, to filyer by."),
	output_file: Path = Option("ycomp-analysis-str.csv", "--output", "-o", dir_okay = False, help = "The (CSV) file to write the analysis to."),
) -> None:
	"""Compare matches in the STR database."""

	kits_df = load_db(kits_str_path)
	if kits_df is None:
		secho(f"ERROR: Kits STR database does not exist.", fg = colors.RED, err = True)
		raise Exit(1)

	echo(f"Found {len(kits_df):,} kits in STR database.")

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

	kits_loci_df = kits_df.iloc[:, 4 :]

	try:
		self_kit_s: pd.Series = kits_loci_df.loc[self_kit]
		kits_df.drop(self_kit, axis = 0, inplace = True)
		kits_loci_df = kits_loci_df.drop(self_kit, axis = 0)
	except KeyError:
		secho(f"ERROR: Kit {self_kit} not found.", fg = colors.RED, err = True)
		raise Exit(1)

	echo(f"Starting analysis...")

	total_num_loci = kits_df.shape[1]

	def get_abs_diffs(kit_s: pd.Series) -> pd.Series:
		def get_diff(a: List[int], b: List[int]) -> List[int]:
			if a is None or b is None:
				return []

			def get_comb_diff(diff: List[Tuple[int, int, int]]) -> List[int]:
				val_diff = [x for _, _, x in diff]

				for i in set(range(len(a))) - set(i for i, _, _ in diff):
					x = 1000
					if i > 0:
						x = min(x, abs(a[i] - a[i - 1]))
					if i < len(a) - 1:
						x = min(x, abs(a[i] - a[i + 1]))
					val_diff.append(1 + x)

				for j in set(range(len(b))) - set(j for _, j, _ in diff):
					x = 1000
					if j > 0:
						x = min(x, abs(b[j] - b[j - 1]))
					if j < len(b) - 1:
						x = min(x, abs(b[j] - b[j + 1]))
					val_diff.append(1 + x)

				return val_diff

			def get_val_diff(a: int, b: int) -> int:
				if a == 0 or b == 0:
					return 1
				else:
					return abs(a - b)

			min_len = min(len(a), len(b))
			min_diff = get_comb_diff(min(
				(
					[(i, j, get_val_diff(a_val, b_val)) for (i, a_val), (j, b_val) in zip(a_comb, b_comb)]
					for a_comb in distinct_combinations(enumerate(a), min_len)
					for b_comb in distinct_combinations(enumerate(b), min_len)
				),
				key = lambda diff: sum(get_comb_diff(diff)),
			))

			return min_diff

		diffs = (get_diff(a, b) for a, b in zip(self_kit_s, kit_s))
		return pd.Series(flatten(diffs))

	from decimal import Decimal
	from statistics import NormalDist

	def get_confidence_interval(confidence_level: SupportsFloat) -> float:
		return abs(NormalDist().inv_cdf((1.0 - float(confidence_level)) / 2))

	rel_dist_cl = Decimal(95) / 100

	loci_diffs = kits_loci_df.apply(lambda kit_s: get_abs_diffs(kit_s), axis = 1)
	num_comps = loci_diffs.count(axis = 1)
	abs_dists = loci_diffs.sum(axis = 1)
	rel_dists = loci_diffs.mean(axis = 1)
	# Calculate standard errors of relative distances, including finite-population correction.
	rel_dist_std_errs = loci_diffs.std(axis = 1) \
		/ np.sqrt(num_comps) \
		* np.sqrt((total_num_loci - num_comps) / (total_num_loci - 1))
	# Calculate confidence intervals using standard scores.
	rel_dist_cis = rel_dist_std_errs * get_confidence_interval(rel_dist_cl)
	rel_dist_mins = rel_dists - rel_dist_cis
	rel_dist_maxs = rel_dists + rel_dist_cis

	results_df = pd.concat(
		[
			kits_df["Group"],
			kits_df["Paternal Ancestor Name"],
			kits_df["Country"],
			kits_df["Haplogroup"],
			num_comps.astype("int32").rename("Loci Compared", inplace = True),
			abs_dists.astype("int32").rename("Absolute Distance", inplace = True),
			rel_dists.astype("float32").rename("Relative Distance", inplace = True),
			rel_dist_cis.astype("float32").rename(f"Relative Distance (CI, {rel_dist_cl:%})", inplace = True),
			rel_dist_mins.astype("float32").rename("Relative Distance (min)", inplace = True),
			rel_dist_maxs.astype("float32").rename("Relative Distance (max)", inplace = True),
		],
		axis = 1,
	)

	echo(f"Finished analysis.")

	results_df.sort_values("Relative Distance (max)").to_csv(output_file)
	echo(f"Analysis written to `{output_file}`.")
