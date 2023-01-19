from itertools import takewhile
from pathlib import Path
from typing import *

import pandas as pd
from platformdirs import PlatformDirs

from .common import *


NDFrameT = TypeVar("NDFrameT", bound=pd.core.generic.NDFrame)


tree_path = "ycomp-yfull-tree.parquet"
snps_path = "ycomp-snps.parquet"
kits_snp_path = "ycomp-kits-snp.parquet"
kits_str_path = "ycomp-kits-str.parquet"


def float_not(a: NDFrameT) -> NDFrameT:
    return cast(NDFrameT, 1 - a)


def float_and(a: NDFrameT, b: NDFrameT) -> NDFrameT:
    return a * b


def float_or(a: NDFrameT, b: NDFrameT) -> NDFrameT:
    return float_not(float_and(float_not(a.fillna(0)), float_not(b.fillna(0))))


def float_xor(a: NDFrameT, b: NDFrameT) -> NDFrameT:
    return float_and(float_or(a, b), float_not(float_and(a, b)))


def float_isna(a: NDFrameT) -> NDFrameT:
    return a.isna().astype("float16", copy=False)


def float_notna(a: NDFrameT) -> NDFrameT:
    return a.notna().astype("float16", copy=False)


def load_db(path: GenericPath) -> pd.DataFrame:
    try:
        return pd.read_parquet(path)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return None


def store_db(path: GenericPath, df: pd.DataFrame) -> None:
    df.to_parquet(path)


def merge_db(path: GenericPath, new_df: pd.DataFrame) -> None:
    df = load_db(path)
    if df is None:
        df = new_df
    else:
        # Add new rows, removing old duplicates.
        df = pd.concat([df, new_df], axis=0)
        df = df[~df.index.duplicated(keep="last")]

    store_db(path, df)


def delete_db(path: GenericPath) -> None:
    Path(path).unlink(missing_ok=True)


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


def get_common_lineage(
    a_lineage: Sequence[str], b_lineage: Sequence[str]
) -> Sequence[str]:
    return list(
        hgs[0]
        for hgs in takewhile(
            lambda hgs: hgs[0] == hgs[1], zip(reversed(a_lineage), reversed(b_lineage))
        )
    )


def kit_matches_lineage(
    tree_df: pd.DataFrame,
    a_lineage: Sequence[str],
    b: pd.Series,
    *,
    haplogroup: Optional[str],
    haplogroup_max_diff: Optional[int]
) -> bool:
    b_lineage = get_haplogroup_lineage(tree_df, b["Haplogroup"])
    common_lineage = get_common_lineage(a_lineage, b_lineage)
    a_lineage_pos = len(a_lineage) - len(common_lineage)
    b_lineage_pos = len(b_lineage) - len(common_lineage)

    if haplogroup is not None and haplogroup not in b_lineage:
        return False
    if (
        haplogroup_max_diff is not None
        and abs(a_lineage_pos - b_lineage_pos) > haplogroup_max_diff
    ):
        return False
    return True
