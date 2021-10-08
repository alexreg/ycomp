from typing import *

import numpy as np
import pandas as pd
from typer import *

from .common import *


app = Typer()


@app.command()
def info(
	use_snp_db: bool = Option(False, "--snp", help = "Use the local SNP database to look up kits. (default)"),
	use_str_db: bool = Option(False, "--str", help = "Use the local STR database to look up kits."),
	self_haplogroup: Optional[str] = Option(None, "--haplogroup", "-hg", help = "The haplogroup to look up."),
	self_kit: Optional[str] = Option(None, "--kit", "-k", help = "The kit whose haplogroup to look up."),
) -> None:
	"""Get basic information about a haplogroup."""

	tree_df = load_db(tree_path)
	if tree_df is None:
		secho(f"ERROR: Tree database does not exist.", fg = colors.RED, err = True)
		raise Exit(1)

	if self_kit is not None:
		if use_str_db:
			kits_db_path = kits_str_path
			kits_db_kind = "STR"
		else:
			kits_db_path = kits_snp_path
			kits_db_kind = "SNP"

		kits_df = load_db(kits_db_path)
		if kits_df is None:
			secho(f"ERROR: Kits {kits_db_kind} database does not exist.", fg = colors.RED, err = True)
			raise Exit(1)
	else:
		kits_df = None

	self_kit_hg = kits_df.loc[self_kit, "Haplogroup"] if self_kit is not None else self_haplogroup
	self_kit_s = tree_df.loc[self_kit_hg]

	self_kit_lineage = get_haplogroup_lineage(tree_df, self_kit_hg)
	lineage_text = " → ".join(reversed(self_kit_lineage))
	primary_snps = ", ".join(self_kit_s["Primary SNPs"])
	extra_snps = ", ".join(self_kit_s["Extra SNPs"])

	echo(f"Haplogroup: {self_kit_hg}")
	echo(f"Lineage: {lineage_text}")
	echo(f"Primary SNPs: {primary_snps}")
	echo(f"Extra SNPs: {extra_snps}")
