"""Shared defaults and reference data for Tardigrade GWAS analysis."""

from __future__ import annotations

DEFAULT_P_THRESHOLD = 5e-8
DEFAULT_P_COL = "P-VALUE"
DEFAULT_GENE_COL = "MAPPED_GENE"
DEFAULT_CHR_COL = "CHR_ID"
DEFAULT_BP_COL = "CHR_POS"
DEFAULT_MAX_PLOT_POINTS = 120_000

CHR_SIZES_HG38: dict[str, int] = {
    "1": 248956422,
    "2": 242193529,
    "3": 198295559,
    "4": 190214555,
    "5": 181538259,
    "6": 170805979,
    "7": 159345973,
    "8": 145138636,
    "9": 138394717,
    "10": 133797422,
    "11": 135086622,
    "12": 133275309,
    "13": 114364328,
    "14": 107043718,
    "15": 101991189,
    "16": 90338345,
    "17": 83257441,
    "18": 80373285,
    "19": 58617616,
    "20": 64444167,
    "21": 46709983,
    "22": 50818468,
    "X": 156040895,
    "Y": 57227415,
}

CHR_BASE_COLORS: dict[str, str] = {
    "1": "#1f77b4",
    "2": "#ff7f0e",
    "3": "#2ca02c",
    "4": "#d62728",
    "5": "#9467bd",
    "6": "#8c564b",
    "7": "#e377c2",
    "8": "#7f7f7f",
    "9": "#bcbd22",
    "10": "#17becf",
    "11": "#1f77b4",
    "12": "#ff7f0e",
    "13": "#2ca02c",
    "14": "#d62728",
    "15": "#9467bd",
    "16": "#8c564b",
    "17": "#e377c2",
    "18": "#7f7f7f",
    "19": "#bcbd22",
    "20": "#17becf",
    "21": "#aec7e8",
    "22": "#ffbb78",
    "X": "#2ca02c",
    "Y": "#7f7f7f",
}

MISSING_GENE_TOKENS = {
    "",
    "-",
    ".",
    "NA",
    "N/A",
    "NAN",
    "NONE",
    "NULL",
    "NR",
}
