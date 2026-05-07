# Epitope Visualization

Tools for mapping, filtering, and visualizing IEDB epitope data against antigen sequences — with optional PyMOL integration for structural coloring.

## Overview

This repository provides a Python pipeline for:

- **Mapping** IEDB epitopes onto CoV-AbDab antibody entries (CDRH3, CDRL3, VH, VL, and V/J gene-level matching).
- **Parsing** IEDB epitope CSVs into position-level frequency dictionaries.
- **Visualizing** epitope counts per residue position, color-coded by amino acid physicochemical properties.
- **Exporting** PyMOL selections to highlight epitope positions on 3D protein structures.

## Repository Structure

``` Plain Text
epitope_visualization/
├── src/
│   ├── epitope_functions.py      # Core pipeline: parsing, filtering, plotting
│   ├── pymol_functions.py        # PyMOL selection & coloring utilities
│   └── dev_notebook.ipynb        # Development / exploratory notebook
├── CovAbDab_epitopes/            # CSV outputs from CoV-AbDab ↔ IEDB mapping
├── misc_files/                   # Auxiliary data files
├── papers/                       # Related literature
├── aa_physicochemical_properties.csv  # Amino acid property lookup table
└── README.md
```

## Dependencies

Install the required packages:

``` bash
pip install pandas numpy biopython matplotlib
```

For PyMOL functionality, [`pymol`](https://pymol.org/) must be installed and the script must be run from within the PyMOL environment.

## Quick Start

### 1. Epitope Position Analysis

Run the full pipeline from an IEDB CSV to property-colored bar charts:

``` python
from src.epitope_functions import run_epitope_analysis

epitope_csv = "path/to/iedb_export.csv"
antigen_seq = "MFVFLVLLPLVSSQCVNLT..."  # full antigen sequence
properties_table = "aa_physicochemical_properties.csv"

run_epitope_analysis(
    epitope_file=epitope_csv,
    ag_aa_seq=antigen_seq,
    properties_table=properties_table,
    start_pos=300,
    end_pos=600,
)
```

### 2. Individual Pipeline Steps

You can also call the functions individually for finer control:

``` python
from src.epitope_functions import (
    get_pos_dict,
    start_end_pos_dict_filter,
    generate_epi_df,
    aggregate_plot_df,
    plot_all_physical_props,
    plot_max_normalized_property,
)

# Parse and filter
## Individual steps
pos_dict = get_pos_dict("path/to/iedb_export.csv")
filtered = start_end_pos_dict_filter(pos_dict, start=300, end=600)

## One-shot parse & filter function
epi_df = generate_epi_df("path/to/iedb_export.csv", start_pos=300, end_pos=600)

# Aggregate and plot
plot_df = aggregate_plot_df(epi_df)
res_props_df = plot_all_physical_props(plot_df)

prop_cols = [
    "hydrophobicity", "hydrophilicity", "h-bonds",
    "side_chain_vol", "polarity", "polarizability",
    "solv_accessible_surface_area", "net_charge_index", "aa_mass",
] ### default columns of included properties table 

plot_max_normalized_property(res_props_df, prop_cols)
```

### 3. Antigen-Based Property Coloring

Colors by the antigen's own residue identities rather than the epitope residue properties so as to compare the patterns between the two:

``` python
from src.epitope_functions import plot_antigen_property_distribution

plot_antigen_property_distribution(
    plot_df=plot_df,
    ag_aa_seq=antigen_seq,
    prop_cols=prop_cols,
    properties_table="aa_physicochemical_properties.csv",
)
```

## PyMOL Integration

Load `pymol_functions.py` from within a PyMOL session to highlight epitope positions on a 3D structure:

``` pymol
# In the PyMOL command line (after loading your structure):
run src/pymol_functions.py

# Create one selection per epitope row (displayed as cyan cartoon):
create_epitope_selections_from_csv epitope_csv=/path/to/epitopes.csv, sequence="MFVFLVLLP..."

# Create a single combined selection, optionally scoped to a chain:
highlight_all_epitope_positions epitope_csv=/path/to/epitopes.csv, sequence="MFVFLVLLP...", chain="A"
```

### Function Summary

| Function | Purpose |
|---|---|
| `create_epitope_selections_from_csv` | Creates a separate PyMOL selection per CSV row, styled as cyan cartoon |
| `highlight_all_epitope_positions` | Creates a single merged selection from all validated residues across rows |

Both functions validate residues against the provided antigen sequence and convert one-letter amino acid codes to three-letter codes for PyMOL compatibility.

## Physicochemical Properties

The `aa_physicochemical_properties.csv` lookup table contains the following properties for each amino acid:

| Property | Description |
|---|---|
| `hydrophobicity` | Kyte-Doolittle hydrophobicity scale |
| `hydrophilicity` | Hopp-Woods hydrophilicity scale |
| `h-bonds` | Hydrogen bond propensity |
| `side_chain_vol` | Side chain volume |
| `polarity` | Polarity index |
| `polarizability` | Polarizability parameter |
| `solv_accessible_surface_area` | Solvent-accessible surface area |
| `net_charge_index` | Net charge index |
| `aa_mass` | Amino acid molecular weight |

## Data Sources

- **[IEDB](https://www.iedb.org/)** — Immune Epitope Database (epitope CSV exports)
- **[CoV-AbDab](https://opig.stats.ox.ac.uk/webapps/covabdab/)** — Coronavirus Antibody Database

## Author

Isaac Daviet
