## Genomic Landscape & High-Legibility Prevalence Analysis

### Overview of Recent Work
Developed and optimised the analysis pipeline for *Klebsiella pneumoniae* genomic profiles using data exported from PathogenWatch (300 total genomes across high-risk lineages, $n = 50$ per Sequence Type). The processing workflow was restructured to collapse highly granular allele classifications into major, clinically diagnostic resistance families to eliminate visual "barcode noise" and maximize readability for presentation and poster deployment.

### Key Implementation Milestones

1. **Bioinformatics Aggregation Pipeline**
   * Migrated from tracking individual micro-evolutionary point mutations/variants to an overarching structural gene family framework.
   * Utilized string-matching pattern recognition (`tidyverse::str_detect`) to dynamically group raw sub-alleles into core diagnostic lines (e.g., merging `blaNDM-1`, `blaNDM-4`, and `blaNDM-5` into a singular `blaNDM-type` family marker).
   * Preserved multi-drug resistance characteristics by structurally tracking cross-resistance profiles across overlapping functional classes (e.g., maintaining the bifunctional enzyme `aac(6')-Ib-cr` across both Aminoglycosides and Fluoroquinolones).

2. **Figure legibility & Scale Enhancements**
   * **Enlarged Structural Indicators:** Significantly scaled up vertical facet panel labels (Drug Classes) and horizontal text sizing to ensure rapid skimmability from a distance on an A0 landscape layout.
   * **Automated Nomenclature Compliance:** Implemented a dynamic `plotmath` translation engine using expression parsing (`bquote`). This system automatically enforces strict genetic nomenclature, isolating and italicizing the beta-lactamase prefix (`italic("bla")`) while preserving standard roman text formatting for family suffixes (e.g., *bla*KPC-type).
   * **Optical Alignment Rectification:** Reconfigured the canvas anchoring context (`plot.title.position = "plot"`), shifting the entire title/subtitle block from the inner data matrix bound to the absolute geometric center of the figure canvas.
   * **Legend Visibility Amplification:** Physically expanded the dimensions of the viridis prevalence scale color bar and scaled up scale text labels to establish a bold, high-contrast visual anchor.

### Consolidated Mapping Matrix

The visual space was optimized by applying the following logical classification structure:

| Raw Kleborate Feature (Granular) | Grouped Gene Family (Row Label) | Structural Panel Grouping |
| :--- | :--- | :--- |
| `blaNDM-1`, `blaNDM-4`, `blaNDM-5`... | ***bla*NDM-type** | Carbapenems |
| `blaKPC-2`, `blaKPC-3`... | ***bla*KPC-type** | Carbapenems |
| `blaOXA-181`, `blaOXA-232`, `blaOXA-48` | ***bla*OXA-48-like** | Carbapenems |
| `blaCTX-M-15`, `blaCTX-M-14`, `blaCTX-M-27` | ***bla*CTX-M-type** | Cephalosporins (ESBL) |
| `blaOXA-1`, `blaOXA-9` | ***bla*OXA-1/9-like** | Other $\beta$-lactams |
| `strA`, `strB` | ***strA/strB*** | Aminoglycosides |
| `aac(6')-Ib-cr` | ***aac(6')-Ib-cr*** | Aminoglycosides / Fluoroquinolones |
| `sul1`, `sul2` | ***sul1*** / ***sul2*** | Sulphonamides |
| `dfrA12`, `dfrA14` | ***dfrA12*** / ***dfrA14*** | Trimethoprim |

### Output Profile
* **File Generated:** `kp_resistance_heatmap_clean.png`
* **Target Dimensions:** W: 11.0 in, H: 6.8 in (300 DPI high-resolution export)
