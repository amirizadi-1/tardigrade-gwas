"""Matplotlib-based Circos-style and similarity heatmap visualizations."""

from __future__ import annotations

from typing import Sequence, cast

import numpy as np
import pandas as pd
from matplotlib import cm, colors as mcolors
from matplotlib.axes import Axes
from matplotlib.collections import PathCollection
from matplotlib.colors import Colormap
from matplotlib.figure import Figure
from matplotlib.projections.polar import PolarAxes
from matplotlib.text import Annotation

from .analysis import _validate_parameters
from .constants import (
    CHR_BASE_COLORS,
    CHR_SIZES_HG38,
    DEFAULT_BP_COL,
    DEFAULT_CHR_COL,
    DEFAULT_GENE_COL,
    DEFAULT_MAX_PLOT_POINTS,
    DEFAULT_P_COL,
    DEFAULT_P_THRESHOLD,
)
from .genes import _negative_log10_series, split_genes
from .io import _normalize_chromosomes, _resolve_column, read_any_table
from .models import GeneLabelArtist

def _chr_sort_key(chromosome: str) -> int:
    if chromosome == "X":
        return 23
    if chromosome == "Y":
        return 24
    try:
        return int(chromosome)
    except ValueError:
        return 10_000


def _pick_gene(row: pd.Series, gene_columns: Sequence[str]) -> str | None:
    for column in gene_columns:
        value = row.get(column)
        genes = split_genes(value)
        if genes:
            return genes[0]
    return None


def _mix_with_white(
    rgb: tuple[float, float, float], amount: float
) -> tuple[float, float, float]:
    amount = float(np.clip(amount, 0.0, 1.0))
    red, green, blue = rgb
    return (
        red + (1.0 - red) * amount,
        green + (1.0 - green) * amount,
        blue + (1.0 - blue) * amount,
    )


def _mix_with_black(
    rgb: tuple[float, float, float], amount: float
) -> tuple[float, float, float]:
    amount = float(np.clip(amount, 0.0, 1.0))
    red, green, blue = rgb
    return (red * (1.0 - amount), green * (1.0 - amount), blue * (1.0 - amount))


def _make_major_bands(
    chromosome_length: int,
    rng: np.random.Generator,
    target_count: int = 14,
    minimum_bp: int = 2_800_000,
    maximum_bp: int = 9_500_000,
    sigma: float = 0.42,
) -> list[tuple[int, int]]:
    """Create deterministic schematic segments for the chromosome ring."""
    mean_length = max(minimum_bp, int(chromosome_length / max(target_count, 1)))
    bands: list[tuple[int, int]] = []
    position = 0

    while position < chromosome_length:
        length = int(
            np.clip(
                rng.lognormal(mean=np.log(mean_length), sigma=sigma),
                minimum_bp,
                maximum_bp,
            )
        )
        end = min(chromosome_length, position + length)
        bands.append((position, end))
        position = end

    return bands


def _draw_chr_ideogram(
    ax: PolarAxes,
    chromosome: str,
    offsets: dict[str, int],
    genome_length: int,
    ring_inner_radius: float,
    ring_outer_radius: float,
    base_color: str,
    random_seed: int = 2026,
    gap_radians: float = 0.040,
) -> None:
    chromosome_start = offsets[chromosome]
    chromosome_length = CHR_SIZES_HG38[chromosome]
    chromosome_end = chromosome_start + chromosome_length

    theta_start = 2.0 * np.pi * (chromosome_start / genome_length)
    theta_end = 2.0 * np.pi * (chromosome_end / genome_length)
    width = max(0.0, (theta_end - theta_start) - gap_radians)
    center = (theta_start + theta_end) / 2.0
    left_edge = center - width / 2.0
    right_edge = center + width / 2.0

    base_rgb = cast(tuple[float, float, float], mcolors.to_rgb(base_color))
    outline = _mix_with_black(base_rgb, 0.55)
    background = _mix_with_white(base_rgb, 0.22)
    ring_height = ring_outer_radius - ring_inner_radius

    ax.bar(
        x=center,
        width=width,
        height=ring_height,
        bottom=ring_inner_radius,
        color=background,
        edgecolor=outline,
        linewidth=2.1,
        alpha=1.0,
        zorder=2,
    )
    ax.plot(
        [left_edge, left_edge],
        [ring_inner_radius, ring_outer_radius],
        linewidth=2.0,
        alpha=0.9,
        color=outline,
        zorder=5,
    )
    ax.plot(
        [right_edge, right_edge],
        [ring_inner_radius, ring_outer_radius],
        linewidth=2.0,
        alpha=0.9,
        color=outline,
        zorder=5,
    )

    stable_chromosome_seed = random_seed + _chr_sort_key(chromosome) * 104_729
    rng = np.random.default_rng(stable_chromosome_seed)
    bands = _make_major_bands(chromosome_length, rng=rng)

    levels = [
        _mix_with_white(base_rgb, 0.52),
        _mix_with_white(base_rgb, 0.38),
        _mix_with_white(base_rgb, 0.26),
        _mix_with_black(base_rgb, 0.10),
        _mix_with_black(base_rgb, 0.22),
        _mix_with_black(base_rgb, 0.34),
    ]
    probabilities = np.array([1.0, 1.2, 1.5, 1.0, 0.85, 0.65], dtype=float)
    probabilities /= probabilities.sum()
    band_edge = _mix_with_black(base_rgb, 0.64)

    for band_start, band_end in bands:
        genome_start = chromosome_start + band_start
        genome_end = chromosome_start + band_end
        band_theta_start = 2.0 * np.pi * (genome_start / genome_length)
        band_theta_end = 2.0 * np.pi * (genome_end / genome_length)

        band_theta_start = float(np.clip(band_theta_start, left_edge, right_edge))
        band_theta_end = float(np.clip(band_theta_end, left_edge, right_edge))
        segment_width = max(
            0.0, (band_theta_end - band_theta_start) - gap_radians * 0.10
        )
        if segment_width <= 0.0:
            continue

        color = levels[int(rng.choice(len(levels), p=probabilities))]
        ax.bar(
            x=(band_theta_start + band_theta_end) / 2.0,
            width=segment_width,
            height=ring_height,
            bottom=ring_inner_radius,
            color=color,
            edgecolor=band_edge,
            linewidth=0.40,
            alpha=1.0,
            zorder=3,
        )


def _downsample_for_plot(
    dataframe: pd.DataFrame,
    p_column: str,
    maximum_points: int,
    preserve_top: int,
) -> pd.DataFrame:
    """Downsample while always retaining the strongest associations."""
    if maximum_points <= 0:
        raise ValueError("Maximum plot points must be greater than zero.")
    if len(dataframe) <= maximum_points:
        return dataframe.copy()

    top_count = min(max(preserve_top, 500), maximum_points)
    strongest = dataframe.nsmallest(top_count, p_column)
    remaining = dataframe.drop(index=strongest.index.tolist())
    random_count = maximum_points - len(strongest)

    if random_count <= 0:
        return strongest.copy()

    sampled = remaining.sample(
        n=min(random_count, len(remaining)),
        random_state=42,
        replace=False,
    )
    return pd.concat([strongest, sampled], axis=0)


def _get_jet_colormap() -> Colormap:
    """Return Matplotlib's jet colormap without importing `matplotlib.colormaps`."""
    return cast(Colormap, getattr(cm, "jet"))


def draw_circular_on_ax(
    gwas_path: str,
    ax: PolarAxes,
    chromosome_column: str = DEFAULT_CHR_COL,
    position_column: str = DEFAULT_BP_COL,
    p_column: str = DEFAULT_P_COL,
    gene_column: str = DEFAULT_GENE_COL,
    p_threshold: float = DEFAULT_P_THRESHOLD,
    maximum_points: int = DEFAULT_MAX_PLOT_POINTS,
    label_top_n: int = 12,
    random_seed: int = 2026,
) -> list[GeneLabelArtist]:
    """Draw a Circos-like GWAS plot and return clickable gene labels."""
    _validate_parameters(p_threshold=p_threshold, top_n=max(label_top_n, 1))

    raw = read_any_table(gwas_path)
    chromosome_column = _resolve_column(
        raw,
        chromosome_column,
        aliases=("CHR", "CHROM", "CHROMOSOME", "CHR_ID"),
    )
    position_column = _resolve_column(
        raw,
        position_column,
        aliases=("BP", "POS", "POSITION", "CHR_POS", "BASE_PAIR_LOCATION"),
    )
    p_column = _resolve_column(
        raw,
        p_column,
        aliases=("P", "PVALUE", "P_VALUE", "P-VALUE"),
    )

    optional_gene_columns: list[str] = []
    for candidate in (gene_column, "REPORTED GENE(S)", "MAPPED_GENE"):
        try:
            resolved = _resolve_column(raw, candidate)
        except ValueError:
            continue
        if resolved not in optional_gene_columns:
            optional_gene_columns.append(resolved)

    selected_columns = [chromosome_column, position_column, p_column]
    selected_columns.extend(
        column for column in optional_gene_columns if column not in selected_columns
    )
    dataframe = raw.loc[:, selected_columns].copy()

    dataframe[chromosome_column] = _normalize_chromosomes(dataframe[chromosome_column])
    dataframe[position_column] = pd.to_numeric(
        dataframe[position_column].astype("string").str.split(";").str[0],
        errors="coerce",
    )
    dataframe[p_column] = pd.to_numeric(dataframe[p_column], errors="coerce")
    dataframe = dataframe.dropna(
        subset=[chromosome_column, position_column, p_column]
    ).copy()

    valid_p = dataframe[p_column].between(0.0, 1.0, inclusive="both")
    dataframe = dataframe.loc[valid_p].copy()
    dataframe = dataframe.loc[dataframe[chromosome_column].isin(CHR_SIZES_HG38)].copy()

    chromosome_limits = dataframe[chromosome_column].map(CHR_SIZES_HG38)
    valid_position = dataframe[position_column].gt(0) & dataframe[position_column].le(
        chromosome_limits
    )
    dataframe = dataframe.loc[valid_position].copy()

    if dataframe.empty:
        raise ValueError(
            "No valid rows remained after chromosome, position, and P-value filtering."
        )

    dataframe[position_column] = dataframe[position_column].astype(np.int64)
    dataframe["_mlogp"] = _negative_log10_series(dataframe[p_column])

    dataframe = _downsample_for_plot(
        dataframe,
        p_column=p_column,
        maximum_points=maximum_points,
        preserve_top=max(label_top_n * 10, 500),
    )

    chromosomes = sorted(
        dataframe[chromosome_column].astype(str).unique(), key=_chr_sort_key
    )
    offsets: dict[str, int] = {}
    cumulative_position = 0
    for chromosome in chromosomes:
        offsets[chromosome] = cumulative_position
        cumulative_position += CHR_SIZES_HG38[chromosome]
    genome_length = cumulative_position

    dataframe["_genome_bp"] = (
        dataframe[chromosome_column].map(lambda chromosome: offsets[str(chromosome)])
        + dataframe[position_column]
    )
    dataframe["_theta"] = 2.0 * np.pi * dataframe["_genome_bp"] / genome_length

    ax.clear()
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2.0)

    ring_inner_radius = 12.5
    ring_outer_radius = 14.0
    point_inner_radius = 2.0
    point_outer_radius = 11.5

    for chromosome in chromosomes:
        _draw_chr_ideogram(
            ax=ax,
            chromosome=chromosome,
            offsets=offsets,
            genome_length=genome_length,
            ring_inner_radius=ring_inner_radius,
            ring_outer_radius=ring_outer_radius,
            base_color=CHR_BASE_COLORS.get(chromosome, "#1f77b4"),
            random_seed=random_seed,
            gap_radians=0.030,
        )
        start = offsets[chromosome] / genome_length * 2.0 * np.pi
        end = (
            (offsets[chromosome] + CHR_SIZES_HG38[chromosome])
            / genome_length
            * 2.0
            * np.pi
        )
        ax.text(
            (start + end) / 2.0,
            ring_outer_radius + 0.65,
            chromosome,
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=10,
        )

    significance = -np.log10(p_threshold)
    mlogp = dataframe["_mlogp"].to_numpy(dtype=float)
    percentile_cap = float(np.nanpercentile(mlogp, 99.5))
    radial_cap = max(percentile_cap, significance, 1.0)
    clipped_mlogp = np.clip(mlogp, 0.0, radial_cap)
    radial_range = point_outer_radius - point_inner_radius
    point_radii = point_outer_radius - (clipped_mlogp / radial_cap) * radial_range

    color_values = clipped_mlogp / radial_cap
    colors_array = _get_jet_colormap()(color_values)
    ax.scatter(
        dataframe["_theta"].to_numpy(dtype=float),
        point_radii,
        s=12,
        c=colors_array,
        alpha=0.90,
        linewidths=0,
    )

    threshold_radius = (
        point_outer_radius - (min(significance, radial_cap) / radial_cap) * radial_range
    )
    full_circle = np.linspace(0.0, 2.0 * np.pi, 720)
    ax.plot(
        full_circle,
        np.full_like(full_circle, threshold_radius),
        linewidth=1.2,
        alpha=0.90,
        color="#2c3e50",
    )

    polar_spine = ax.spines.get("polar")
    if polar_spine is not None:
        polar_spine.set_visible(False)
    ax.set_frame_on(False)

    label_artists: list[GeneLabelArtist] = []
    if optional_gene_columns and label_top_n > 0:
        top = dataframe.nsmallest(max(label_top_n * 6, label_top_n), p_column).copy()
        top["_gene"] = top.apply(
            lambda row: _pick_gene(row, optional_gene_columns), axis=1
        )
        top = (
            top.dropna(subset=["_gene"])
            .drop_duplicates(subset=["_gene"])
            .head(label_top_n)
        )

        if not top.empty:
            top_mlogp = np.clip(top["_mlogp"].to_numpy(dtype=float), 0.0, radial_cap)
            top_point_radii = (
                point_outer_radius - (top_mlogp / radial_cap) * radial_range
            )
            label_dot_radii = np.full(len(top), ring_outer_radius + 1.8)
            label_collection = cast(
                PathCollection,
                ax.scatter(
                    top["_theta"].to_numpy(dtype=float),
                    label_dot_radii,
                    s=40,
                    c="#000000",
                    picker=5,
                    zorder=11,
                ),
            )

            bbox_style = {
                "facecolor": "#ffffff",
                "edgecolor": "none",
                "boxstyle": "round,pad=0.25",
                "alpha": 0.95,
            }

            for point_index, (_, row) in enumerate(top.iterrows()):
                gene = str(row["_gene"])
                theta = float(row["_theta"])
                point_radius = float(top_point_radii[point_index])
                dot_radius = float(label_dot_radii[point_index])
                ax.plot(
                    [theta, theta],
                    [point_radius, dot_radius],
                    color="gray",
                    linewidth=0.7,
                    alpha=0.6,
                    zorder=9,
                )
                annotation = cast(
                    Annotation,
                    ax.annotate(
                        gene,
                        xy=(theta, dot_radius),
                        xytext=(0, 12),
                        textcoords="offset points",
                        horizontalalignment="center",
                        verticalalignment="bottom",
                        fontsize=9,
                        bbox=bbox_style,
                        visible=False,
                    ),
                )
                label_artists.append(
                    GeneLabelArtist(
                        collection=label_collection,
                        point_index=point_index,
                        annotation=annotation,
                        gene=gene,
                    )
                )

    ax.set_ylim(0.0, ring_outer_radius + 4.0)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(alpha=0.12)
    ax.set_title("GWAS Circos-like plot", pad=22)
    return label_artists


def build_similarity_matrix(
    similarity_table: pd.DataFrame,
    names: Sequence[str],
) -> pd.DataFrame:
    name_index = pd.Index(list(names), dtype="object")
    matrix = pd.DataFrame(
        np.eye(len(name_index), dtype=float),
        index=name_index,
        columns=name_index.copy(),
    )

    if not similarity_table.empty:
        rows = similarity_table.loc[:, ["A", "B", "jaccard"]].itertuples(
            index=False, name=None
        )
        for name_a, name_b, value in rows:
            if name_a in matrix.index and name_b in matrix.columns:
                matrix.loc[name_a, name_b] = float(value)
                matrix.loc[name_b, name_a] = float(value)
    return matrix


def draw_similarity_heatmap(
    ax: Axes,
    similarity_table: pd.DataFrame,
    names: Sequence[str],
    title: str = "Disease similarity (Jaccard)",
) -> None:
    matrix = build_similarity_matrix(similarity_table, names)
    ax.clear()
    ax.imshow(matrix.values, vmin=0.0, vmax=1.0, cmap="Blues")
    ax.set_xticks(range(len(names)))
    ax.set_yticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, horizontalalignment="right", fontsize=9)
    ax.set_yticklabels(names, fontsize=9)

    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            ax.text(
                column_index,
                row_index,
                f"{matrix.iat[row_index, column_index]:.2f}",
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=8,
            )
    ax.set_title(title, fontsize=11)

