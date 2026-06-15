"""Core multi-GWAS overlap analysis."""

from __future__ import annotations

import itertools
from collections import Counter
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from .constants import DEFAULT_GENE_COL, DEFAULT_P_COL, DEFAULT_P_THRESHOLD
from .genes import (
    common_table,
    extract_gene_counts,
    extract_gene_set,
    extract_gene_weighted_scores,
    jaccard_similarity,
)
from .io import _resolve_column, read_any_table
from .models import AnalysisResult, Dataset, PairKey

def _validate_parameters(p_threshold: float, top_n: int) -> None:
    if not np.isfinite(p_threshold) or not 0.0 < p_threshold <= 1.0:
        raise ValueError(
            "P threshold must be a finite number greater than 0 and at most 1."
        )
    if top_n <= 0:
        raise ValueError("Top N must be greater than zero.")


def _validate_datasets(datasets: Sequence[Dataset]) -> list[Dataset]:
    if not datasets:
        raise ValueError("No datasets were provided.")

    cleaned: list[Dataset] = []
    seen_names: set[str] = set()
    for raw_name, raw_path in datasets:
        name = raw_name.strip()
        path = str(Path(raw_path).expanduser())
        if not name:
            raise ValueError("Dataset names cannot be empty.")
        if name in seen_names:
            raise ValueError(f"Dataset name is duplicated: {name}")
        if not Path(path).is_file():
            raise FileNotFoundError(f"Dataset file not found: {path}")
        seen_names.add(name)
        cleaned.append((name, path))
    return cleaned


def run_multi_gwas(
    datasets: Sequence[Dataset],
    p_column: str = DEFAULT_P_COL,
    gene_column: str = DEFAULT_GENE_COL,
    p_threshold: float = DEFAULT_P_THRESHOLD,
    top_n: int = 20,
) -> AnalysisResult:
    """Analyze significant genes, overlaps, similarities, and consensus scores."""
    _validate_parameters(p_threshold=p_threshold, top_n=top_n)
    cleaned_datasets = _validate_datasets(datasets)
    names = [name for name, _ in cleaned_datasets]

    data: dict[str, pd.DataFrame] = {}
    original_row_counts: dict[str, int] = {}
    valid_p_counts: dict[str, int] = {}
    zero_p_counts: dict[str, int] = {}
    resolved_columns: dict[str, tuple[str, str]] = {}

    for name, path in cleaned_datasets:
        dataframe = read_any_table(path)
        original_row_counts[name] = len(dataframe)
        resolved_p = _resolve_column(
            dataframe,
            p_column,
            aliases=("P", "PVALUE", "P_VALUE", "P-VALUE"),
        )
        resolved_gene = _resolve_column(
            dataframe,
            gene_column,
            aliases=("MAPPED_GENE", "REPORTED GENE(S)", "GENE", "GENES"),
        )

        dataframe = dataframe.copy()
        dataframe[resolved_p] = pd.to_numeric(dataframe[resolved_p], errors="coerce")
        valid_mask = dataframe[resolved_p].between(0.0, 1.0, inclusive="both")
        valid_p_counts[name] = int(valid_mask.sum())
        zero_p_counts[name] = int(dataframe.loc[valid_mask, resolved_p].eq(0.0).sum())
        dataframe = dataframe.loc[valid_mask].copy()

        data[name] = dataframe
        resolved_columns[name] = (resolved_p, resolved_gene)

    filtered: dict[str, pd.DataFrame] = {}
    gene_sets: dict[str, set[str]] = {}
    gene_counts: dict[str, pd.Series] = {}
    gene_weighted_scores: dict[str, pd.Series] = {}

    for name, dataframe in data.items():
        resolved_p, resolved_gene = resolved_columns[name]
        significant = dataframe.loc[dataframe[resolved_p] <= p_threshold].copy()
        filtered[name] = significant
        gene_sets[name] = extract_gene_set(significant[resolved_gene])
        gene_counts[name] = extract_gene_counts(significant[resolved_gene])
        gene_weighted_scores[name] = extract_gene_weighted_scores(
            significant,
            p_column=resolved_p,
            gene_column=resolved_gene,
        )

    pairwise_common_sets: dict[PairKey, set[str]] = {}
    pairwise_top_tables: dict[PairKey, pd.DataFrame] = {}
    pairwise_full_tables: dict[PairKey, pd.DataFrame] = {}

    for name_a, name_b in itertools.combinations(names, 2):
        key = (name_a, name_b)
        common = gene_sets[name_a] & gene_sets[name_b]
        pairwise_common_sets[key] = common
        if common:
            full_table = common_table(
                common,
                gene_counts[name_a],
                gene_counts[name_b],
                name_a,
                name_b,
            )
            pairwise_full_tables[key] = full_table
            pairwise_top_tables[key] = full_table.head(top_n)
        else:
            empty = pd.DataFrame()
            pairwise_full_tables[key] = empty
            pairwise_top_tables[key] = empty

    all_common_set: set[str] = set()
    all_common_table = pd.DataFrame()
    if len(names) >= 2:
        all_common_set = set.intersection(*(gene_sets[name] for name in names))
        if all_common_set:
            index = sorted(all_common_set)
            all_common_table = pd.DataFrame(
                {
                    f"{name}_count": gene_counts[name]
                    .reindex(index)
                    .fillna(0)
                    .astype(int)
                    for name in names
                },
                index=pd.Index(index, dtype="object"),
            )
            all_common_table["total"] = all_common_table.sum(axis=1)
            all_common_table = all_common_table.sort_values(
                ["total", *(f"{name}_count" for name in names)],
                ascending=False,
            )

    similarity_rows = [
        (name_a, name_b, jaccard_similarity(gene_sets[name_a], gene_sets[name_b]))
        for name_a, name_b in itertools.combinations(names, 2)
    ]
    similarity_table = (
        pd.DataFrame(
            similarity_rows,
            columns=pd.Index(["A", "B", "jaccard"], dtype="object"),
        )
        .sort_values("jaccard", ascending=False)
        .reset_index(drop=True)
        if similarity_rows
        else pd.DataFrame(columns=pd.Index(["A", "B", "jaccard"], dtype="object"))
    )

    all_genes_union: set[str] = set()
    for name in names:
        all_genes_union.update(gene_counts[name].index.astype(str))
        all_genes_union.update(gene_weighted_scores[name].index.astype(str))

    consensus_table = pd.DataFrame(
        index=pd.Index(sorted(all_genes_union), dtype="object")
    )
    for name in names:
        consensus_table[f"{name}_count"] = (
            gene_counts[name].reindex(consensus_table.index).fillna(0).astype(int)
        )
        consensus_table[f"{name}_wscore"] = (
            gene_weighted_scores[name]
            .reindex(consensus_table.index)
            .fillna(0.0)
            .astype(float)
        )

    if not consensus_table.empty:
        presence: Counter[str] = Counter()
        for name in names:
            presence.update(gene_sets[name])
        consensus_table["presence"] = (
            pd.Series(presence, dtype="int64")
            .reindex(consensus_table.index)
            .fillna(0)
            .astype(int)
        )
        consensus_table["total_count"] = consensus_table[
            [f"{name}_count" for name in names]
        ].sum(axis=1)
        consensus_table["total_wscore"] = consensus_table[
            [f"{name}_wscore" for name in names]
        ].sum(axis=1)
        consensus_table = consensus_table.sort_values(
            ["presence", "total_wscore", "total_count"],
            ascending=False,
        )

    per_dataset_top_counts: dict[str, pd.DataFrame] = {}
    per_dataset_top_weighted: dict[str, pd.DataFrame] = {}
    for name in names:
        counts = gene_counts[name]
        weighted = gene_weighted_scores[name]
        per_dataset_top_counts[name] = (
            counts.head(top_n).to_frame(name="count")
            if not counts.empty
            else pd.DataFrame(columns=pd.Index(["count"], dtype="object"))
        )
        per_dataset_top_weighted[name] = (
            weighted.head(top_n).to_frame(name="weighted_score")
            if not weighted.empty
            else pd.DataFrame(columns=pd.Index(["weighted_score"], dtype="object"))
        )

    lines: list[str] = [
        "=== Multi-GWAS Gene Overlap Report ===",
        f"Requested P-value column: {p_column}",
        f"Requested gene column:    {gene_column}",
        f"P threshold:              {p_threshold:g}",
        f"Datasets:                 {len(names)}",
        "",
        "---- Per-dataset summary ----",
    ]

    for name in names:
        lines.append(
            f"{name}: rows={original_row_counts[name]} | valid_p_rows={valid_p_counts[name]} "
            f"| zero_p_rows={zero_p_counts[name]} "
            f"| significant_rows={len(filtered[name])} "
            f"| unique_sig_genes={len(gene_sets[name])}"
        )
    lines.append("")

    lines.append("---- Per-dataset TOP genes ----")
    for name in names:
        lines.append(f"[{name}] Top {top_n} genes by COUNT:")
        count_table = per_dataset_top_counts[name]
        lines.append(count_table.to_string() if not count_table.empty else "(none)")
        lines.append("")
        lines.append(f"[{name}] Top {top_n} genes by WEIGHTED SCORE (sum(-log10(p))):")
        weighted_table = per_dataset_top_weighted[name]
        lines.append(
            weighted_table.to_string() if not weighted_table.empty else "(none)"
        )
        lines.append("")

    if len(names) == 1:
        lines.append(
            "Only one dataset was provided; pairwise and all-dataset overlaps do not apply."
        )
    else:
        lines.append("---- Pairwise overlaps ----")
        for key, common in pairwise_common_sets.items():
            name_a, name_b = key
            lines.append(f"{name_a} ∩ {name_b}: common_genes={len(common)}")
            top_table = pairwise_top_tables[key]
            if not top_table.empty:
                lines.append(
                    f"Top {min(top_n, len(top_table))} shared genes by count evidence:"
                )
                lines.append(top_table.to_string())
            lines.append("")

        lines.extend(
            [
                "---- Disease similarity (Jaccard on significant gene sets) ----",
                similarity_table.to_string(index=False)
                if not similarity_table.empty
                else "(not applicable)",
                "",
                "---- Cross-dataset consensus ----",
            ]
        )
        if not consensus_table.empty:
            lines.append("Top 20 genes by presence, weighted score, and count:")
            lines.append(consensus_table.head(20).to_string())
            lines.append("")
            lines.append("Genes present in at least two datasets (Top 50):")
            present_in_two = consensus_table.loc[consensus_table["presence"] >= 2].head(
                50
            )
            lines.append(
                present_in_two.to_string() if not present_in_two.empty else "(none)"
            )
        else:
            lines.append("(none)")
        lines.append("")

        lines.append("---- Intersection of ALL datasets ----")
        lines.append(f"{' ∩ '.join(names)}: common_genes={len(all_common_set)}")
        if not all_common_table.empty:
            lines.append("Top genes shared by all datasets:")
            lines.append(all_common_table.head(20).to_string())
        lines.append("")

    return {
        "report_text": "\n".join(lines),
        "similarity_table": similarity_table,
        "consensus_table": consensus_table,
        "all_common_table": all_common_table,
        "pairwise_full_tables": pairwise_full_tables,
    }

