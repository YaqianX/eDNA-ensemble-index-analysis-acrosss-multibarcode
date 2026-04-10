# eDNA Multi-Barcode Ensemble Index Pipeline

A reproducible pipeline for computing per-taxon **eDNA ensemble indices** across multiple sampling matrices (biofilm, membrane filter, sediment) and molecular barcodes (16S rRNA, 18S rRNA, COI, 12S rRNA).

> **Reference methodology:**  
> Kelly, R. P., et al. (2019). *Sampling unit matters more than barcode choice in coastal eDNA metabarcoding.* Molecular Ecology Resources, 19(6), 1604–1617.  
> Djurhuus, A., et al. (2020). *Environmental DNA reveals seasonal shifts and potential interactions in a marine community.* Nature Communications, 11(1), 254.

---

## Table of Contents

- [Pipeline Overview](#pipeline-overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Script Reference](#script-reference)
- [Complete Workflow](#complete-workflow)
- [Output Files](#output-files)
- [Running the Integrated Pipeline](#running-the-integrated-pipeline)
- [Citation](#citation)

---

## Pipeline Overview

```
Per barcode (repeat for each of: 16S, 18S, COI, 12S)
─────────────────────────────────────────────────────────────────────
 [otutab.txt]  +  [taxonomy.txt]
        │
        ▼
 merge_otu_taxonomy.py        →  merged OTU table with taxonomy columns
        │
        ▼
 merge_otu_taxon.py           →  taxon-level abundance table (OTUs collapsed)
        │
        ▼
 otutab_normalize.py          →  relative abundance table (0–1 per sample)
        │
        ▼  (one normalized table per barcode)

Cross-barcode integration
─────────────────────────────────────────────────────────────────────
 [16S_normalize.txt]
 [18S_normalize.txt]   ──►  eDNA_ensemble_index_pipeline.py
 [COI_normalize.txt]            │
 [12S_normalize.txt]            │
                         Step 1a: group samples → matrix means
                         Step 1b: divide by max  → eDNA index (0–1)
                         Step 2:  average across barcodes (non-zero only)
                                → ensemble index per taxon per matrix
```

**Key methodological decisions implemented:**

| Decision | Implementation |
|---|---|
| Relative abundance normalization | Per-sample column sum division (0–1) |
| eDNA index (within-barcode) | Taxon value ÷ maximum value across matrices |
| Ensemble index (cross-barcode) | Mean of non-zero eDNA index values; zeros = true non-detection, excluded |
| COI biofilm exclusion | `--exclude-matrix COI:Biofilm`; excluded matrix set to 0, not averaged |
| Taxon identity | Highest resolved level (Species > Genus > … > Kingdom); "Unassigned" dropped |

---

## Requirements

**Python:** 3.9 or higher

**Dependencies:**

```
pandas >= 1.5.0
numpy >= 1.23.0
```

Install with:

```bash
pip install pandas numpy
```

**Input file format:**  
All input files must be tab-separated (`.txt` or `.tsv`).  
The first column must be OTU/ASV identifiers (`#OTU_ID`, `OTUID`, or `OTU_ID`).  
Sample column names should begin with matrix prefix letters (e.g., `B1`, `B2` for Biofilm; `M1`, `M2` for Membrane; `S1`, `S2` for Sediment).

---

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/edna-ensemble-pipeline.git
cd edna-ensemble-pipeline
pip install pandas numpy
```

No additional compilation or installation is required. All scripts are standalone Python modules.

---

## Script Reference

### 1. `merge_otu_taxonomy.py`

Merges a taxonomy annotation file with an OTU abundance table, aligning records by OTU/ASV identifier. Produces a combined table with taxonomy columns prepended.

**Usage:**

```bash
python merge_otu_taxonomy.py \
    -t taxonomy.txt \
    -o otutab.txt \
    -out merged_taxonomy_otutab.txt
```

**Arguments:**

| Flag | Description | Required |
|---|---|---|
| `-t / --taxonomy` | Path to taxonomy annotation file | ✓ |
| `-o / --otutab` | Path to OTU abundance table | ✓ |
| `-out / --output` | Path for merged output file | ✓ |

**Input taxonomy file format:**

| #OTU_ID | Kingdom | Phylum | Class | Order | Family | Genus | Species |
|---|---|---|---|---|---|---|---|
| ASV001 | Bacteria | Cyanobacteria | … | … | … | … | … |

**Output:** Tab-separated file with columns `#OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species, [sample columns]`.

---

### 2. `merge_otu_taxon.py`

Aggregates OTU-level abundances to the finest resolved taxonomic level. Each OTU is assigned to exactly one taxon (Species if annotated, otherwise Genus, Family, …, Kingdom). Unresolved OTUs labeled "Unassigned" at all levels are discarded.

**Usage:**

```bash
python merge_otu_taxon.py \
    -i merged_taxonomy_otutab.txt \
    -o taxon_abundance.txt \
    -l taxon_abundance.log
```

**Arguments:**

| Flag | Description | Required |
|---|---|---|
| `-i / --input` | Merged taxonomy + OTU table (output of Step 1) | ✓ |
| `-o / --output` | Path for taxon-level abundance table | ✓ |
| `-l / --log` | Log file path (default: same name as output + `.log`) | — |

**Output:** Tab-separated file with columns `Classification, Taxon, [sample columns]`.  
`Classification` records the taxonomic level at which the taxon was resolved (e.g., `Genus`, `Family`).

---

### 3. `otutab_normalize.py`

Converts raw read counts to within-sample relative abundances by dividing each value by its column total. Supports both proportional (0–1) and percentage (0–100) output.

**Usage:**

```bash
python otutab_normalize.py \
    -i taxon_abundance.txt \
    -o taxon_normalized.txt \
    -m relative
```

**Arguments:**

| Flag | Description | Default |
|---|---|---|
| `-i / --input` | Input abundance table (output of Step 2) | required |
| `-o / --output` | Output file path | auto-generated |
| `-s / --separator` | Field delimiter (`\t`, `,`, `;`, ` `) | `\t` |
| `-m / --method` | `relative` (0–1) or `percentage` (0–100) | `relative` |
| `-f / --format` | Numeric output format | `%.8f` |
| `-t / --top` | Display top N taxa by mean abundance | `10` |
| `-q / --quiet` | Suppress progress output | `False` |

> **Note:** Always use `-m relative` (the default). The eDNA ensemble index pipeline requires proportional relative abundances (0–1) as input.

---

### 4. `eDNA_ensemble_index_pipeline.py`

The core integration script. Accepts normalized relative abundance tables from one or more barcodes and computes the eDNA ensemble index following a two-step procedure:

- **Step 1a** – Average technical/biological replicates within each sampling matrix per barcode.
- **Step 1b** – Divide by the row maximum across matrices to produce the per-barcode eDNA index (0–1 scale).
- **Step 2** – Average non-zero eDNA index values across barcodes for each taxon–matrix combination. Zero values (true non-detections) are excluded from the cross-barcode mean.

**Mode A — Full pipeline (from normalized tables):**

```bash
python eDNA_ensemble_index_pipeline.py \
    --inputs 18S:18S_normalized.txt COI:COI_normalized.txt \
    --matrix-map "B:Biofilm,M:Membrane,S:Sediment" \
    --exclude-matrix COI:Biofilm \
    --outdir results/
```

**Mode B — Step 2 only (from pre-computed Step 1b files):**

```bash
python eDNA_ensemble_index_pipeline.py \
    --step1-inputs 18S:18S_step1b.txt COI:COI_step1b.txt \
    --exclude-matrix COI:Biofilm \
    --output-merged ensemble_final.txt \
    --outdir results/
```

**Arguments:**

| Flag | Description | Default |
|---|---|---|
| `--inputs` | Full pipeline entry: `BARCODE:path` pairs | — |
| `--step1-inputs` | Step 2-only entry: pre-computed Step 1b files | — |
| `--matrix-map` | Sample prefix → matrix name mapping | `B:Biofilm,M:Membrane,S:Sediment` |
| `--exclude-matrix` | Barcode–matrix pairs to exclude (e.g., `COI:Biofilm`) | none |
| `--outdir / -o` | Output directory (created if absent) | `edna_output/` |
| `--prefix` | Uniform prefix for all output filenames | none |
| `--output-merged` | Custom filename for the final ensemble index file | auto |
| `--out-step1a` | Custom filenames for Step 1a outputs: `BARCODE:filename,...` | auto |
| `--out-step1b` | Custom filenames for Step 1b outputs: `BARCODE:filename,...` | auto |

**Output Status labels:**

| Status value | Meaning |
|---|---|
| `Merged_Average_(Sources:18S+COI)` | Taxon detected by multiple barcodes; ensemble index is mean of non-zero values |
| `Unique_to_18S` | Taxon detected by a single barcode only |
| `Not_detected` | Taxon present in union but zero across all barcodes and matrices |

---

## Complete Workflow

### Example: Four barcodes, three matrices

The following example reflects the study design used in this analysis (Yongding River, 2023–2024), with barcodes 16S rRNA, 18S rRNA, COI, and 12S rRNA and three sampling matrices (biofilm, membrane filter, sediment). COI and 12S rRNA lack biofilm amplification data and are excluded from that matrix.

**Step 1 — Per-barcode preprocessing (repeat for each barcode):**

```bash
# Example shown for 18S rRNA; repeat identically for 16S, COI, 12S
# substituting the appropriate input files.

python merge_otu_taxonomy.py \
    -t data/18S_taxonomy.txt \
    -o data/18S_otutab.txt \
    -out intermediate/18S_merged.txt

python merge_otu_taxon.py \
    -i intermediate/18S_merged.txt \
    -o intermediate/18S_taxon.txt

python otutab_normalize.py \
    -i intermediate/18S_taxon.txt \
    -o intermediate/18S_normalized.txt \
    -m relative
```

**Step 2 — Ensemble index computation:**

```bash
python eDNA_ensemble_index_pipeline.py \
    --inputs \
        16S:intermediate/16S_normalized.txt \
        18S:intermediate/18S_normalized.txt \
        COI:intermediate/COI_normalized.txt \
        12S:intermediate/12S_normalized.txt \
    --matrix-map "B:Biofilm,M:Membrane,S:Sediment" \
    --exclude-matrix COI:Biofilm 12S:Biofilm \
    --prefix Yongding_2024 \
    --outdir results/ensemble/
```

---

## Output Files

After running the complete workflow, the output directory contains the following files for each barcode:

```
results/ensemble/
├── Yongding_2024_16S_step1a_matrix_mean.txt   # Matrix-averaged relative abundance
├── Yongding_2024_16S_step1b_eDNA_index.txt    # Per-barcode eDNA index (0–1)
├── Yongding_2024_18S_step1a_matrix_mean.txt
├── Yongding_2024_18S_step1b_eDNA_index.txt
├── Yongding_2024_COI_step1a_matrix_mean.txt
├── Yongding_2024_COI_step1b_eDNA_index.txt
├── Yongding_2024_12S_step1a_matrix_mean.txt
├── Yongding_2024_12S_step1b_eDNA_index.txt
└── Yongding_2024_ensemble_index_merged.txt    # Final ensemble index table
```

**Final ensemble index table format:**

| Taxon | Status | Biofilm | Membrane | Sediment |
|---|---|---|---|---|
| *Microcystis aeruginosa* | Merged_Average_(Sources:16S+18S) | 0.85 | 0.42 | 0.61 |
| *Gammarus* sp. | Unique_to_COI | 0.00 | 0.73 | 0.95 |
| *Cyprinus carpio* | Unique_to_12S | 0.00 | 0.38 | 0.00 |

All numeric values are written in scientific notation (`%.8e`) to preserve precision for low-abundance taxa.

---

## Running the Integrated Pipeline

For convenience, a single-entry wrapper script `run_pipeline.py` is provided. It executes all four scripts sequentially for a given barcode list and produces the final ensemble index in one command.

```bash
python run_pipeline.py \
    --barcodes 16S 18S COI 12S \
    --otutab-dir data/ \
    --taxonomy-dir data/ \
    --workdir intermediate/ \
    --outdir results/ensemble/ \
    --matrix-map "B:Biofilm,M:Membrane,S:Sediment" \
    --exclude-matrix COI:Biofilm 12S:Biofilm \
    --prefix Yongding_2024
```

Expected input files in `--otutab-dir` and `--taxonomy-dir`:

```
data/
├── 16S_otutab.txt
├── 16S_taxonomy.txt
├── 18S_otutab.txt
├── 18S_taxonomy.txt
├── COI_otutab.txt
├── COI_taxonomy.txt
├── 12S_otutab.txt
└── 12S_taxonomy.txt
```

---

## File Naming Conventions

| Convention | Example |
|---|---|
| Raw OTU table | `{BARCODE}_otutab.txt` |
| Taxonomy annotation | `{BARCODE}_taxonomy.txt` |
| Merged table | `{BARCODE}_merged.txt` |
| Taxon-level table | `{BARCODE}_taxon.txt` |
| Normalized table | `{BARCODE}_normalized.txt` |
| Matrix mean (Step 1a) | `{PREFIX}_{BARCODE}_step1a_matrix_mean.txt` |
| eDNA index (Step 1b) | `{PREFIX}_{BARCODE}_step1b_eDNA_index.txt` |
| Ensemble index (Step 2) | `{PREFIX}_ensemble_index_merged.txt` |

---

## Citation

If you use this pipeline in your research, please cite the following:

> [Your citation here upon publication]

And the underlying methodological references:

> Kelly, R. P., Gallego, R., & Jacobs-Palmer, E. (2019). The effect of tides on nearshore environmental DNA. *PeerJ*, 6, e4521.

> Djurhuus, A., Closek, C. J., Kelly, R. P., Pitz, K. J., Michisaki, R. P., Starks, H. A., … & Chavez, F. P. (2020). Environmental DNA reveals seasonal shifts and potential interactions in a marine community. *Nature Communications*, 11(1), 254.

---

## License

This pipeline is distributed for academic use. Please contact the authors before redistributing or incorporating into commercial software.
