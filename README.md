# Tardigrade GWAS

**Tardigrade GWAS** is a lightweight Python library for comparing gene-level signals across multiple GWAS summary-statistics files. It provides significant-gene extraction, pairwise overlap analysis, Jaccard similarity, cross-dataset consensus scoring, report generation, similarity heatmaps, and Circos-style GWAS visualization.

Developed by **[Amir Izadi](https://github.com/amirizadi-1)** and **[Zahra Aghabeygi](https://github.com/zahraaghabeygi)**.

> Tardigrade GWAS is a research utility. It is not a clinical diagnostic system and should not be used as the sole basis for medical decisions.

## Features

- Read CSV, TSV, TXT, `.csv.gz`, and `.tsv.gz` GWAS tables
- Detect common GWAS column aliases case-insensitively
- Filter associations by a configurable P-value threshold
- Parse multi-gene fields separated by commas or semicolons
- Count significant-gene occurrences per dataset
- Calculate gene evidence scores using `sum(-log10(P))`
- Calculate pairwise gene-set intersections
- Calculate pairwise Jaccard similarity
- Build a cross-dataset consensus table
- Find genes shared across all datasets
- Generate text reports and CSV result tables
- Generate Jaccard similarity heatmaps
- Generate Circos-style chromosome plots
- Return standard Pandas DataFrames and Matplotlib Figures for further customization

## Installation

```bash
pip install tardigrade-gwas
```

Upgrade to the latest available version:

```bash
pip install --upgrade tardigrade-gwas
```

Check the installed version:

```python
import tardigrade_gwas as tg

print(tg.__version__)
```

## Quick start

```python
import tardigrade_gwas as tg

result = tg.analyze_files(
    ["disease_a.csv", "disease_b.tsv"],
    names=["Disease A", "Disease B"],
)

print(result.report_text)
print(result.similarity_table)
print(result.consensus_table.head())
```

Save the generated report and all result tables:

```python
result.save_report("outputs/report.txt")
result.save_tables("outputs/tables")
```

## Input data

### Supported file types

- `.csv`
- `.tsv`
- `.txt`
- `.csv.gz`
- `.tsv.gz`

For `.txt` files or files without a recognized CSV/TSV extension, Tardigrade attempts to detect comma or tab delimiters automatically.

### Minimum columns for overlap analysis

The default analysis columns are:

| Purpose | Default column |
|---|---|
| P value | `P-VALUE` |
| Gene annotation | `MAPPED_GENE` |

Example:

```csv
P-VALUE,MAPPED_GENE
1.2e-10,TP53
4.8e-9,BRCA1
2.3e-6,APOE
7.1e-12,"IL6;STAT3"
```

Only rows satisfying `P-value <= p_threshold` are included in the significant-gene analysis. The default threshold is `5e-8`.

### Columns required for Circos-style plots

| Purpose | Default column |
|---|---|
| Chromosome | `CHR_ID` |
| Base-pair position | `CHR_POS` |
| P value | `P-VALUE` |
| Gene annotation | `MAPPED_GENE` |

Example:

```csv
P-VALUE,MAPPED_GENE,CHR_ID,CHR_POS
1.2e-10,TP53,17,7676594
4.8e-9,BRCA1,17,43044295
3.1e-12,APOE,19,44908684
```

### Recognized column aliases

Column matching is case-insensitive.

| Data type | Recognized names |
|---|---|
| P value | `P-VALUE`, `P`, `PVALUE`, `P_VALUE` |
| Gene | `MAPPED_GENE`, `REPORTED GENE(S)`, `GENE`, `GENES` |
| Chromosome | `CHR_ID`, `CHR`, `CHROM`, `CHROMOSOME` |
| Position | `CHR_POS`, `BP`, `POS`, `POSITION`, `BASE_PAIR_LOCATION` |

Custom column names can always be passed explicitly.

## Main analysis API

### `analyze_files`

Use `analyze_files` when you have a list of file paths.

```python
from pathlib import Path
import tardigrade_gwas as tg

result = tg.analyze_files(
    files=[
        Path("data/alzheimer.csv"),
        Path("data/parkinson.csv"),
        Path("data/multiple_sclerosis.tsv"),
    ],
    names=[
        "Alzheimer Disease",
        "Parkinson Disease",
        "Multiple Sclerosis",
    ],
    p_column="P-VALUE",
    gene_column="MAPPED_GENE",
    p_threshold=5e-8,
    top_n=20,
)
```

Parameters:

| Parameter | Type | Default | Description |
|---|---|---:|---|
| `files` | sequence of paths | required | GWAS files to analyze |
| `names` | sequence of strings or `None` | `None` | Display names for the datasets |
| `p_column` | string | `P-VALUE` | Preferred P-value column |
| `gene_column` | string | `MAPPED_GENE` | Preferred gene column |
| `p_threshold` | float | `5e-8` | Significance threshold in `(0, 1]` |
| `top_n` | integer | `20` | Number of top genes shown in report sections |

When `names` is omitted, names are derived from filenames. Duplicate derived names receive suffixes such as `_2`, `_3`, and so on.

### `analyze_datasets`

Use `analyze_datasets` when each dataset already has an explicit display name.

```python
import tardigrade_gwas as tg

result = tg.analyze_datasets(
    [
        ("Alzheimer Disease", "data/ad.csv"),
        ("Parkinson Disease", "data/pd.csv"),
    ],
    p_threshold=5e-8,
)
```

Dataset names must be non-empty and unique.

### `run_multi_gwas`

`run_multi_gwas` is the lower-level analysis function. It returns a dictionary rather than a `MultiGWASResult` object.

```python
raw_result = tg.run_multi_gwas(
    datasets=[
        ("Disease A", "a.csv"),
        ("Disease B", "b.csv"),
    ],
    p_column="P-VALUE",
    gene_column="MAPPED_GENE",
    p_threshold=5e-8,
    top_n=20,
)

print(raw_result["report_text"])
```

For normal use, `analyze_files` or `analyze_datasets` is recommended.

## Understanding the result object

`analyze_files` and `analyze_datasets` return a `MultiGWASResult` object.

```python
result = tg.analyze_files(["a.csv", "b.csv"])
```

The object contains:

```python
result.dataset_names
result.report_text
result.similarity_table
result.consensus_table
result.all_common_table
result.pairwise_full_tables
```

### `dataset_names`

Tuple containing dataset names in input order:

```python
print(result.dataset_names)
```

### `report_text`

A complete plain-text report containing:

- original row counts
- valid P-value counts
- zero P-value counts
- significant row counts
- unique significant-gene counts
- top genes by occurrence count
- top genes by weighted evidence score
- pairwise overlaps
- Jaccard similarity values
- cross-dataset consensus results
- genes shared by all datasets

```python
print(result.report_text)
```

Save it:

```python
path = result.save_report("outputs/multi_gwas_report.txt")
print(path)
```

### `similarity_table`

Pairwise Jaccard similarities between significant-gene sets:

```python
print(result.similarity_table)
```

The Jaccard similarity is:

```text
Jaccard(A, B) = |A intersection B| / |A union B|
```

Its value ranges from `0` to `1`:

- `0`: no shared significant genes
- `1`: identical significant-gene sets

The table has the columns:

```text
A, B, jaccard
```

### `consensus_table`

A gene-level table combining evidence across all input datasets:

```python
print(result.consensus_table.head(20))
```

For every dataset, it contains:

- `<dataset>_count`: number of significant rows containing the gene
- `<dataset>_wscore`: sum of `-log10(P)` for rows containing the gene

It also contains:

- `presence`: number of datasets in which the gene is significant
- `total_count`: total significant-row count across all datasets
- `total_wscore`: total weighted score across all datasets

The table is sorted by `presence`, then `total_wscore`, then `total_count`, all in descending order.

Filter genes present in at least two datasets:

```python
shared = result.consensus_table[
    result.consensus_table["presence"] >= 2
]

print(shared)
```

Get the ten genes with the strongest combined weighted evidence:

```python
top_weighted = result.consensus_table.nlargest(10, "total_wscore")
print(top_weighted)
```

### `all_common_table`

Genes found in every input dataset:

```python
if result.all_common_table.empty:
    print("No gene is shared by all datasets.")
else:
    print(result.all_common_table)
```

This table is only applicable when at least two datasets are analyzed.

### `pairwise_full_tables`

Dictionary of complete pairwise overlap tables:

```python
for (name_a, name_b), table in result.pairwise_full_tables.items():
    print(f"{name_a} vs {name_b}")
    print(table.head())
```

Access one specific comparison:

```python
pair_table = result.pairwise_full_tables[
    ("Disease A", "Disease B")
]
```

If the two datasets have no shared significant genes, the corresponding DataFrame is empty.

### Convert the result to a dictionary

```python
result_dict = result.as_dict()
```

The returned dictionary contains:

```text
report_text
similarity_table
consensus_table
all_common_table
pairwise_full_tables
```

## Save all result tables

```python
written_files = result.save_tables("outputs/tables")

for path in written_files:
    print(path)
```

Output structure:

```text
outputs/tables/
├── similarity.csv
├── consensus.csv
├── all_common.csv
└── pairwise/
    ├── Disease_A__Disease_B.csv
    └── Disease_A__Disease_C.csv
```

## Similarity heatmap

### High-level API

```python
import matplotlib.pyplot as plt
import tardigrade_gwas as tg

result = tg.analyze_files(
    ["disease_a.csv", "disease_b.csv", "disease_c.csv"],
    names=["Disease A", "Disease B", "Disease C"],
)

figure = tg.plot_similarity_heatmap(
    result,
    title="Disease similarity based on significant genes",
    figsize=(9, 7),
)

plt.show()
```

The plot is not saved automatically. Save it only when needed:

```python
figure.savefig(
    "outputs/similarity_heatmap.png",
    dpi=300,
    bbox_inches="tight",
)
```

Other formats supported by Matplotlib include PDF and SVG:

```python
figure.savefig("outputs/heatmap.pdf", bbox_inches="tight")
figure.savefig("outputs/heatmap.svg", bbox_inches="tight")
```

You can pass a custom ordering or subset of names:

```python
figure = tg.plot_similarity_heatmap(
    result,
    names=["Disease C", "Disease A", "Disease B"],
)
```

### Build the similarity matrix without plotting

```python
matrix = tg.build_similarity_matrix(
    result.similarity_table,
    result.dataset_names,
)

print(matrix)
```

### Draw on an existing Matplotlib axis

```python
import matplotlib.pyplot as plt
import tardigrade_gwas as tg

fig, ax = plt.subplots(figsize=(8, 6))

tg.draw_similarity_heatmap(
    ax,
    similarity_table=result.similarity_table,
    names=result.dataset_names,
    title="Disease similarity (Jaccard)",
)

fig.tight_layout()
plt.show()
```

## Circos-style GWAS plot

### High-level API

```python
import matplotlib.pyplot as plt
import tardigrade_gwas as tg

figure = tg.plot_circos(
    "data/disease.csv",
    chromosome_column="CHR_ID",
    position_column="CHR_POS",
    p_column="P-VALUE",
    gene_column="MAPPED_GENE",
    p_threshold=5e-8,
    maximum_points=120_000,
    label_top_n=12,
    random_seed=2026,
    figsize=(10, 10),
)

plt.show()
```

Save the plot:

```python
figure.savefig(
    "outputs/disease_circos.png",
    dpi=300,
    bbox_inches="tight",
)
```

Parameters:

| Parameter | Default | Description |
|---|---:|---|
| `chromosome_column` | `CHR_ID` | Chromosome column |
| `position_column` | `CHR_POS` | Base-pair position column |
| `p_column` | `P-VALUE` | P-value column |
| `gene_column` | `MAPPED_GENE` | Gene column used for labels |
| `p_threshold` | `5e-8` | Significance threshold ring |
| `maximum_points` | `120000` | Maximum plotted points after downsampling |
| `label_top_n` | `12` | Number of top unique genes prepared for annotation |
| `random_seed` | `2026` | Seed controlling deterministic schematic chromosome bands |
| `figsize` | `(9, 9)` | Figure size in inches |

Set `label_top_n=0` to disable gene annotations.

### Draw on an existing polar axis

```python
import matplotlib.pyplot as plt
import tardigrade_gwas as tg

fig, ax = plt.subplots(
    figsize=(10, 10),
    subplot_kw={"projection": "polar"},
)

label_artists = tg.draw_circular_on_ax(
    gwas_path="data/disease.csv",
    ax=ax,
    p_threshold=5e-8,
    maximum_points=120_000,
    label_top_n=12,
)

fig.tight_layout()
plt.show()
```

`draw_circular_on_ax` returns a list of `GeneLabelArtist` objects. Their annotations are initially hidden. To display all generated labels:

```python
for artist in label_artists:
    artist.annotation.set_visible(True)

fig.canvas.draw_idle()
```

### Circos data rules

- Supported chromosomes: `1` through `22`, `X`, and `Y`
- Values such as `chr1`, `chrX`, and `1.0` are normalized
- Positions must be numeric, greater than zero, and within the hg38 chromosome length
- P values must be numeric and between `0` and `1`
- The plot uses hg38 chromosome lengths
- The chromosome band pattern is schematic and is not a true cytoband ideogram
- `maximum_points` reduces the number of plotted points, but the file is still read into memory before downsampling

## Read a GWAS table directly

```python
import tardigrade_gwas as tg

df = tg.read_any_table("data/disease.tsv.gz")

print(df.head())
print(df.columns.tolist())
print(df.shape)
```

The function returns a Pandas DataFrame.

## Gene parsing and utility functions

### Split a gene field

```python
import tardigrade_gwas as tg

genes = tg.split_genes("TP53; BRCA1, APOE")
print(genes)
```

Output:

```python
["TP53", "BRCA1", "APOE"]
```

Commas and semicolons are supported. Duplicate genes within one cell are removed while original order is preserved.

The following tokens are treated as missing values:

```text
"", -, ., NA, N/A, NAN, NONE, NULL, NR
```

### Extract a unique gene set

```python
gene_set = tg.extract_gene_set(df["MAPPED_GENE"])
print(len(gene_set))
```

### Count gene occurrences

```python
counts = tg.extract_gene_counts(df["MAPPED_GENE"])
print(counts.head(20))
```

### Calculate weighted gene scores

```python
scores = tg.extract_gene_weighted_scores(
    df,
    p_column="P-VALUE",
    gene_column="MAPPED_GENE",
)

print(scores.head(20))
```

For each valid row, the weight is:

```text
-log10(P-value)
```

If a row contains multiple genes, the complete row weight is assigned to every listed gene; it is not divided among them.

### Calculate Jaccard similarity directly

```python
set_a = {"TP53", "APOE", "BRCA1"}
set_b = {"TP53", "APOE", "IL6"}

score = tg.jaccard_similarity(set_a, set_b)
print(score)  # 0.5
```

If both sets are empty, the function returns `0.0`.

## Custom column names

Suppose your table uses:

```text
PVALUE, GENE_SYMBOL, CHROMOSOME, BASE_PAIR
```

Analyze it with:

```python
result = tg.analyze_files(
    ["custom.csv"],
    p_column="PVALUE",
    gene_column="GENE_SYMBOL",
)
```

Create a Circos-style plot with:

```python
figure = tg.plot_circos(
    "custom.csv",
    p_column="PVALUE",
    gene_column="GENE_SYMBOL",
    chromosome_column="CHROMOSOME",
    position_column="BASE_PAIR",
)
```

## Error handling

```python
import tardigrade_gwas as tg

try:
    result = tg.analyze_files(["a.csv", "b.csv"])
except FileNotFoundError as error:
    print(f"File error: {error}")
except ValueError as error:
    print(f"Input error: {error}")
```

Common errors include:

- file not found
- empty input file
- malformed CSV/TSV table
- missing required column
- duplicate explicit dataset name
- number of supplied names not matching the number of files
- invalid `p_threshold`
- invalid `top_n`
- no valid chromosome/position/P-value rows for Circos plotting

## Complete workflow example

```python
from pathlib import Path
import matplotlib.pyplot as plt
import tardigrade_gwas as tg

output_dir = Path("outputs")
output_dir.mkdir(parents=True, exist_ok=True)

files = [
    Path("data/alzheimer.csv"),
    Path("data/parkinson.csv"),
    Path("data/multiple_sclerosis.tsv"),
]

names = [
    "Alzheimer Disease",
    "Parkinson Disease",
    "Multiple Sclerosis",
]

result = tg.analyze_files(
    files,
    names=names,
    p_column="P-VALUE",
    gene_column="MAPPED_GENE",
    p_threshold=5e-8,
    top_n=20,
)

print(result.report_text)

result.save_report(output_dir / "multi_gwas_report.txt")
result.save_tables(output_dir / "tables")

heatmap = tg.plot_similarity_heatmap(
    result,
    title="GWAS gene-set similarity",
    figsize=(9, 7),
)
heatmap.savefig(
    output_dir / "similarity_heatmap.png",
    dpi=300,
    bbox_inches="tight",
)

for file_path, dataset_name in zip(files, names, strict=True):
    safe_name = dataset_name.lower().replace(" ", "_")

    circos = tg.plot_circos(
        file_path,
        p_threshold=5e-8,
        maximum_points=120_000,
        label_top_n=12,
        figsize=(10, 10),
    )

    circos.savefig(
        output_dir / f"{safe_name}_circos.png",
        dpi=300,
        bbox_inches="tight",
    )

plt.show()
```

## Public API summary

### Analysis

```python
tg.analyze_files
tg.analyze_datasets
tg.run_multi_gwas
```

### Input

```python
tg.read_any_table
```

### Gene utilities

```python
tg.split_genes
tg.extract_gene_set
tg.extract_gene_counts
tg.extract_gene_weighted_scores
tg.jaccard_similarity
```

### Plotting

```python
tg.plot_similarity_heatmap
tg.plot_circos
tg.build_similarity_matrix
tg.draw_similarity_heatmap
tg.draw_circular_on_ax
```

### Defaults

```python
tg.DEFAULT_P_THRESHOLD
tg.DEFAULT_P_COL
tg.DEFAULT_GENE_COL
tg.DEFAULT_CHR_COL
tg.DEFAULT_BP_COL
tg.DEFAULT_MAX_PLOT_POINTS
```

## Scientific and technical limitations

- Gene symbols are compared as exact strings.
- Gene aliases are not harmonized automatically; for example, `PARK2` and `PRKN` remain distinct.
- Gene-symbol capitalization is not normalized; `TP53` and `tp53` may be treated as different genes.
- Genome-build conversion is not performed.
- Circos positions are validated against hg38 chromosome lengths.
- Circos chromosome bands are schematic, not true cytogenetic bands.
- Zero P values are replaced internally by the smallest positive floating-point value when calculating `-log10(P)`; this can produce very large weighted scores.
- Tables are loaded into memory with Pandas; very large GWAS files may require substantial RAM.
- The package does not currently perform LD clumping, meta-analysis, pathway enrichment, gene ontology analysis, SNP annotation, liftover, QQ plots, or Manhattan plots.
- Outputs should be interpreted in the context of the source datasets and their original study designs.

## Authors and Contributions

- **[Amir Izadi](https://github.com/amirizadi-1)** — Technical development, software engineering, analytical implementation, optimization, Python-library development, testing, debugging, and Windows/Linux deployment.
- **[Zahra Aghabeygi](https://github.com/zahraaghabeygi)** — Visual identity, logo and TARDIGRADE wordmark design, color-system development, user-interface structure, dashboard design, and UX direction.

Project ideation and concept development were shared equally between both contributors. AI-assisted tools were used as supporting resources in the technical and UI/UX workflows; all final decisions, implementation, testing, refinement, and approval were carried out by the human contributors.

For a detailed contribution statement, see [CONTRIBUTORS.md](CONTRIBUTORS.md).

## License

Tardigrade GWAS is distributed under the MIT License.
