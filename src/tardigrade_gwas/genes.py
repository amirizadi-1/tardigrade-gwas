"""Gene parsing, counting, weighting, and set-similarity helpers."""

from __future__ import annotations

import math
from numbers import Real
from typing import cast

import numpy as np
import pandas as pd

from .constants import MISSING_GENE_TOKENS

def _negative_log10_series(series: pd.Series) -> np.ndarray:
    """Return -log10(values), replacing zero with the smallest positive float."""
    values = series.to_numpy(dtype=np.float64)
    safe_values = np.maximum(values, np.finfo(np.float64).tiny)
    return -np.log10(safe_values)


def _is_missing_scalar(value: object) -> bool:
    """Check scalar missing values without calling pandas.isna on an arbitrary object."""
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, Real):
        return math.isnan(float(value))
    return False


def split_genes(cell: object) -> list[str]:
    """Split a gene field, remove placeholders, and preserve unique order."""
    if _is_missing_scalar(cell):
        return []

    unique_genes: dict[str, None] = {}
    for raw_gene in str(cell).replace(";", ",").split(","):
        gene = raw_gene.strip()
        if gene.upper() in MISSING_GENE_TOKENS:
            continue
        unique_genes.setdefault(gene, None)
    return list(unique_genes)


def extract_gene_set(series: pd.Series) -> set[str]:
    genes: set[str] = set()
    for cell in series.dropna():
        genes.update(split_genes(cell))
    return genes


def extract_gene_counts(series: pd.Series) -> pd.Series:
    all_genes: list[str] = []
    for cell in series.dropna():
        all_genes.extend(split_genes(cell))
    if not all_genes:
        return pd.Series(dtype="int64")
    return pd.Series(all_genes, dtype="string").value_counts()


def extract_gene_weighted_scores(
    dataframe: pd.DataFrame,
    p_column: str,
    gene_column: str,
) -> pd.Series:
    if dataframe.empty:
        return pd.Series(dtype="float64")

    temporary = dataframe.loc[:, [p_column, gene_column]].copy()
    temporary[p_column] = pd.to_numeric(temporary[p_column], errors="coerce")
    temporary = temporary.dropna(subset=[p_column, gene_column])
    temporary = temporary.loc[
        temporary[p_column].between(0.0, 1.0, inclusive="both")
    ].copy()
    if temporary.empty:
        return pd.Series(dtype="float64")

    temporary["_weight"] = _negative_log10_series(temporary[p_column])
    rows: list[tuple[str, float]] = []
    for _, row in temporary.iterrows():
        genes = split_genes(row[gene_column])
        weight = float(row["_weight"])
        rows.extend((gene, weight) for gene in genes)

    if not rows:
        return pd.Series(dtype="float64")

    weighted_rows = pd.DataFrame(
        rows, columns=pd.Index(["gene", "weight"], dtype="object")
    )
    grouped_scores = cast(
        pd.Series, weighted_rows.groupby("gene", sort=False)["weight"].sum()
    )
    return grouped_scores.sort_values(ascending=False)


def jaccard_similarity(set_a: set[str], set_b: set[str]) -> float:
    union = set_a | set_b
    if not union:
        return 0.0
    return len(set_a & set_b) / len(union)


def common_table(
    common_genes: set[str],
    counts_a: pd.Series,
    counts_b: pd.Series,
    name_a: str,
    name_b: str,
) -> pd.DataFrame:
    index = sorted(common_genes)
    dataframe = pd.DataFrame(
        {
            f"{name_a}_count": counts_a.reindex(index).fillna(0).astype(int),
            f"{name_b}_count": counts_b.reindex(index).fillna(0).astype(int),
        },
        index=pd.Index(index, dtype="object"),
    )
    dataframe["total"] = dataframe.sum(axis=1)
    return dataframe.sort_values(
        ["total", f"{name_a}_count", f"{name_b}_count"],
        ascending=False,
    )

