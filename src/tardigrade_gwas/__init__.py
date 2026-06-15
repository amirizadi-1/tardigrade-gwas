"""Tardigrade: Multi-GWAS overlap analysis and visualization."""

from ._version import __version__
from .api import analyze_datasets, analyze_files, plot_circos, plot_similarity_heatmap
from .analysis import run_multi_gwas
from .constants import (
    DEFAULT_BP_COL,
    DEFAULT_CHR_COL,
    DEFAULT_GENE_COL,
    DEFAULT_MAX_PLOT_POINTS,
    DEFAULT_P_COL,
    DEFAULT_P_THRESHOLD,
)
from .genes import (
    extract_gene_counts,
    extract_gene_set,
    extract_gene_weighted_scores,
    jaccard_similarity,
    split_genes,
)
from .io import read_any_table
from .models import AnalysisResult, Dataset, MultiGWASResult
from .plotting import build_similarity_matrix, draw_circular_on_ax, draw_similarity_heatmap

__all__ = [
    "__version__",
    "AnalysisResult",
    "Dataset",
    "MultiGWASResult",
    "DEFAULT_P_THRESHOLD",
    "DEFAULT_P_COL",
    "DEFAULT_GENE_COL",
    "DEFAULT_CHR_COL",
    "DEFAULT_BP_COL",
    "DEFAULT_MAX_PLOT_POINTS",
    "analyze_datasets",
    "analyze_files",
    "run_multi_gwas",
    "read_any_table",
    "split_genes",
    "extract_gene_set",
    "extract_gene_counts",
    "extract_gene_weighted_scores",
    "jaccard_similarity",
    "plot_circos",
    "plot_similarity_heatmap",
    "build_similarity_matrix",
    "draw_circular_on_ax",
    "draw_similarity_heatmap",
]
