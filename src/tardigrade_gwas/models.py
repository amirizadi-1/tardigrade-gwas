"""Public result models and internal plot specifications."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias, TypedDict

import pandas as pd
from matplotlib.collections import PathCollection
from matplotlib.text import Annotation

Dataset: TypeAlias = tuple[str, str]
PairKey: TypeAlias = tuple[str, str]


class AnalysisResult(TypedDict):
    report_text: str
    similarity_table: pd.DataFrame
    consensus_table: pd.DataFrame
    all_common_table: pd.DataFrame
    pairwise_full_tables: dict[PairKey, pd.DataFrame]


@dataclass(frozen=True, slots=True)
class MultiGWASResult:
    """Structured result returned by the public analysis API."""

    dataset_names: tuple[str, ...]
    report_text: str
    similarity_table: pd.DataFrame
    consensus_table: pd.DataFrame
    all_common_table: pd.DataFrame
    pairwise_full_tables: dict[PairKey, pd.DataFrame]

    @classmethod
    def from_mapping(
        cls,
        dataset_names: tuple[str, ...],
        result: AnalysisResult,
    ) -> "MultiGWASResult":
        return cls(
            dataset_names=dataset_names,
            report_text=result["report_text"],
            similarity_table=result["similarity_table"],
            consensus_table=result["consensus_table"],
            all_common_table=result["all_common_table"],
            pairwise_full_tables=result["pairwise_full_tables"],
        )

    def as_dict(self) -> AnalysisResult:
        """Return the original mapping representation used by the desktop app."""
        return {
            "report_text": self.report_text,
            "similarity_table": self.similarity_table,
            "consensus_table": self.consensus_table,
            "all_common_table": self.all_common_table,
            "pairwise_full_tables": self.pairwise_full_tables,
        }

    def save_report(self, path: str | Path) -> Path:
        destination = Path(path).expanduser()
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(self.report_text, encoding="utf-8")
        return destination

    def save_tables(self, directory: str | Path) -> list[Path]:
        """Save the main result tables as CSV files."""
        output_dir = Path(directory).expanduser()
        output_dir.mkdir(parents=True, exist_ok=True)
        written: list[Path] = []

        tables = {
            "similarity.csv": self.similarity_table,
            "consensus.csv": self.consensus_table,
            "all_common.csv": self.all_common_table,
        }
        for filename, table in tables.items():
            path = output_dir / filename
            table.to_csv(path)
            written.append(path)

        pairwise_dir = output_dir / "pairwise"
        pairwise_dir.mkdir(exist_ok=True)
        for (name_a, name_b), table in self.pairwise_full_tables.items():
            safe_a = _safe_filename(name_a)
            safe_b = _safe_filename(name_b)
            path = pairwise_dir / f"{safe_a}__{safe_b}.csv"
            table.to_csv(path)
            written.append(path)
        return written


def _safe_filename(value: str) -> str:
    cleaned = "".join(character if character.isalnum() or character in "._-" else "_" for character in value)
    return cleaned.strip("._-") or "dataset"


@dataclass(frozen=True, slots=True)
class HeatmapPlotSpec:
    title: str
    result: AnalysisResult
    names: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class CircosPlotSpec:
    title: str
    dataset: Dataset


PlotSpec: TypeAlias = HeatmapPlotSpec | CircosPlotSpec


@dataclass(slots=True)
class GeneLabelArtist:
    collection: PathCollection
    point_index: int
    annotation: Annotation
    gene: str
