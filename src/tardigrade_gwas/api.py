"""High-level public API for analysis and plotting."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence, cast

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.projections.polar import PolarAxes

from .analysis import run_multi_gwas
from .constants import (
    DEFAULT_BP_COL,
    DEFAULT_CHR_COL,
    DEFAULT_GENE_COL,
    DEFAULT_MAX_PLOT_POINTS,
    DEFAULT_P_COL,
    DEFAULT_P_THRESHOLD,
)
from .models import Dataset, MultiGWASResult
from .plotting import draw_circular_on_ax, draw_similarity_heatmap


def analyze_datasets(
    datasets: Sequence[Dataset],
    *,
    p_column: str = DEFAULT_P_COL,
    gene_column: str = DEFAULT_GENE_COL,
    p_threshold: float = DEFAULT_P_THRESHOLD,
    top_n: int = 20,
) -> MultiGWASResult:
    """Analyze explicitly named datasets.

    Args:
        datasets: ``[(display_name, file_path), ...]``.
        p_column: Preferred P-value column name. Matching is case-insensitive and
            common aliases are recognized.
        gene_column: Preferred gene column name.
        p_threshold: Significance threshold in the interval ``(0, 1]``.
        top_n: Number of top genes to include in report sections.
    """
    normalized = [(str(name), str(Path(path).expanduser())) for name, path in datasets]
    raw = run_multi_gwas(
        datasets=normalized,
        p_column=p_column,
        gene_column=gene_column,
        p_threshold=p_threshold,
        top_n=top_n,
    )
    return MultiGWASResult.from_mapping(
        tuple(name for name, _ in normalized),
        raw,
    )


def analyze_files(
    files: Sequence[str | Path],
    *,
    names: Sequence[str] | None = None,
    p_column: str = DEFAULT_P_COL,
    gene_column: str = DEFAULT_GENE_COL,
    p_threshold: float = DEFAULT_P_THRESHOLD,
    top_n: int = 20,
) -> MultiGWASResult:
    """Analyze GWAS files, deriving unique dataset names when omitted."""
    paths = [Path(path).expanduser() for path in files]
    if names is not None and len(names) != len(paths):
        raise ValueError("The number of names must match the number of files.")

    display_names = list(names) if names is not None else _derive_unique_names(paths)
    datasets: list[Dataset] = [
        (str(name), str(path)) for name, path in zip(display_names, paths, strict=True)
    ]
    return analyze_datasets(
        datasets,
        p_column=p_column,
        gene_column=gene_column,
        p_threshold=p_threshold,
        top_n=top_n,
    )


def plot_similarity_heatmap(
    result_or_table: MultiGWASResult | pd.DataFrame,
    *,
    names: Sequence[str] | None = None,
    title: str = "Disease similarity (Jaccard)",
    figsize: tuple[float, float] = (8.0, 6.5),
) -> Figure:
    """Create and return a similarity heatmap figure without saving it."""
    if isinstance(result_or_table, MultiGWASResult):
        similarity_table = result_or_table.similarity_table
        resolved_names = tuple(names) if names is not None else result_or_table.dataset_names
    else:
        similarity_table = result_or_table
        if names is None:
            resolved_names = _names_from_similarity_table(similarity_table)
        else:
            resolved_names = tuple(names)

    if not resolved_names:
        raise ValueError("At least one dataset name is required for a heatmap.")

    figure, axis = plt.subplots(figsize=figsize)
    draw_similarity_heatmap(
        axis,
        similarity_table=similarity_table,
        names=resolved_names,
        title=title,
    )
    figure.tight_layout()
    return figure


def plot_circos(
    gwas_path: str | Path,
    *,
    chromosome_column: str = DEFAULT_CHR_COL,
    position_column: str = DEFAULT_BP_COL,
    p_column: str = DEFAULT_P_COL,
    gene_column: str = DEFAULT_GENE_COL,
    p_threshold: float = DEFAULT_P_THRESHOLD,
    maximum_points: int = DEFAULT_MAX_PLOT_POINTS,
    label_top_n: int = 12,
    random_seed: int = 2026,
    figsize: tuple[float, float] = (9.0, 9.0),
) -> Figure:
    """Create and return a Circos-style GWAS figure without saving it."""
    figure, raw_axis = plt.subplots(figsize=figsize, subplot_kw={"projection": "polar"})
    axis = cast(PolarAxes, raw_axis)
    draw_circular_on_ax(
        gwas_path=str(Path(gwas_path).expanduser()),
        ax=axis,
        chromosome_column=chromosome_column,
        position_column=position_column,
        p_column=p_column,
        gene_column=gene_column,
        p_threshold=p_threshold,
        maximum_points=maximum_points,
        label_top_n=label_top_n,
        random_seed=random_seed,
    )
    figure.tight_layout()
    return figure


def _derive_unique_names(paths: Iterable[Path]) -> list[str]:
    names: list[str] = []
    counts: dict[str, int] = {}
    for path in paths:
        name = path.name
        if name.lower().endswith(".gz"):
            name = name[:-3]
        stem = Path(name).stem.strip() or "dataset"
        count = counts.get(stem, 0) + 1
        counts[stem] = count
        names.append(stem if count == 1 else f"{stem}_{count}")
    return names


def _names_from_similarity_table(table: pd.DataFrame) -> tuple[str, ...]:
    if table.empty:
        return ()
    required = {"A", "B"}
    if not required.issubset(table.columns):
        raise ValueError("Similarity table must contain 'A' and 'B' columns.")
    ordered: dict[str, None] = {}
    for value in table["A"].tolist() + table["B"].tolist():
        ordered.setdefault(str(value), None)
    return tuple(ordered)
