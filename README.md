# RNA-Seq Ophiostoma novo-ulmi

This repository contains the analysis associated with the manuscript “Dynamics of Ophiostoma novo-ulmi transcriptome during colonization of resistant and susceptible Ulmus minor hosts”.
It includes the full downstream RNA-Seq analysis workflow, from processed count matrices to statistical analyses and figure generation.

- Combined_Dual_RNA-Seq_Methodology/

This folder documents the bioinformatic preprocessing pipeline used prior to the R-based analyses.
It contains bash commands and scripts integrating several external programs for:
Read mapping against a chimeric reference genome (host + pathogen)
Filtering of reads based on mapping quality and uniqueness
Separation of host and pathogen reads

Generation of gene-level count matrices (FeatureCounts)

- data/

This folder contains all input files required for the RStudio analysis.
All files listed below must be placed inside data/ for the analysis to run correctly.
The analysis relies on relative paths, so R must be launched from the repository root.

- analysis.Rmd (repository root)

Runnable RMarkdown analysis containing the full statistical workflow used to generate all figures and results presented in the manuscript, including:
Data import and formatting
Normalization and library size correction
Differential expression analyses
Functional annotation integration
Figure generation


- install_packages.R (repository root)

Script to install all required R packages used in the analysis.
This script only needs to be run once per R installation.

## Quick start

Place all required input files into a data/ directory at the repository root.
From an R session launched at the repository root, run:
```
source("install_packages.R")  # installs packages (first time only)
```

Render the analysis:

In RStudio: open analysis.md and click Knit.
From R:
```
rmarkdown::render("analysis.Rmd")
````


## Required input files (place under data/) and data description:

-counts_OPhio_EXP.txt
Gene-level count matrix generated using FeatureCounts.

-Seq_name_match_Complete_name.csv
Mapping between sequencing library identifiers and full sample names.

-Genotecas_order.csv
Library ordering and metadata used for sample organization.

-headers_protein_geneBernier_sinOphio.txt
Header file used for parsing reference annotation files (GFF3 and Excel-based annotations).

-genes_proteins.txt
Gene-to-protein ID correspondence used during annotation parsing.

-GenomeAnnotation_OphioH327.csv
Reference genome annotation for O. novo-ulmi strain H327.

-My_Tabla_Colonizacion_virulencia_REF_RStudio.csv
Curated list of genes associated with virulence-related functions.

-123_all_PHI_hits.csv
Gene hits against the PHI-base database.

-GO.csv
Gene Ontology (GO) annotations.

-121_gene_coords.txt
Genomic coordinates of analyzed genes.

-Chr1_8._H327txt.txt
Chromosome lengths for the O. novo-ulmi genome.

-for_R_library_correction.csv
Library size information used for normalization and correction in R.
