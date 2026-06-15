"""Input helpers for CSV, TSV, TXT, and gzip-compressed GWAS tables."""

from __future__ import annotations

import csv
import gzip
import os
from pathlib import Path
from typing import Sequence

import pandas as pd

def _resolve_column(
    dataframe: pd.DataFrame,
    requested: str,
    aliases: Sequence[str] = (),
) -> str:
    """Resolve a column exactly first, then case-insensitively."""
    candidates = (requested, *aliases)
    columns = [str(column) for column in dataframe.columns]

    for candidate in candidates:
        if candidate in dataframe.columns:
            return candidate

    normalized = {column.strip().casefold(): column for column in columns}
    for candidate in candidates:
        match = normalized.get(candidate.strip().casefold())
        if match is not None:
            return match

    available = ", ".join(columns[:20])
    suffix = " ..." if len(columns) > 20 else ""
    raise ValueError(
        f"Column '{requested}' was not found. Available columns: {available}{suffix}"
    )


def _detect_delimiter(path: Path) -> str:
    """Detect comma or tab, including gzip-compressed text tables."""
    try:
        opener = gzip.open if path.name.lower().endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8-sig", errors="replace") as handle:
            sample = handle.read(8192)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return dialect.delimiter
    except (OSError, csv.Error, EOFError):
        return "\t"


def read_any_table(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Read a GWAS table while respecting CSV/TSV extensions."""
    file_path = Path(path).expanduser()
    if not file_path.is_file():
        raise FileNotFoundError(f"File not found: {file_path}")

    lower_name = file_path.name.lower()
    if lower_name.endswith((".csv", ".csv.gz")):
        separator = ","
    elif lower_name.endswith((".tsv", ".tsv.gz")):
        separator = "\t"
    else:
        separator = _detect_delimiter(file_path)

    try:
        return pd.read_csv(
            file_path,
            sep=separator,
            low_memory=False,
            encoding="utf-8-sig",
        )
    except UnicodeDecodeError:
        return pd.read_csv(
            file_path,
            sep=separator,
            low_memory=False,
            encoding="latin-1",
        )
    except pd.errors.EmptyDataError as exc:
        raise ValueError(f"The selected file is empty: {file_path}") from exc
    except pd.errors.ParserError as exc:
        raise ValueError(f"Could not parse the table: {file_path}\n{exc}") from exc


def _normalize_chromosomes(series: pd.Series) -> pd.Series:
    normalized = (
        series.astype("string")
        .str.strip()
        .str.replace(r"^chr", "", case=False, regex=True)
        .str.replace(r"\.0$", "", regex=True)
        .str.upper()
    )
    return normalized

