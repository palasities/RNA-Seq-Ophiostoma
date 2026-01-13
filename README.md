# RNA-Seq Dynamics — README

This repository contains the analysis for "Dynamics of Ophiostoma novo-ulmi transcriptome during colonization of resistant and susceptible Ulmus minor hosts".

Files:
- analysis.md — runnable RMarkdown analysis.
- install_packages.R — installs required R packages for the analysis.
- README.md  — short instructions.

Quick start
1. Place your required input files into a `data/` directory. 
2. From an R session at the repository root, run:
   - source("install_packages.R")    # installs packages (first time only)
3. Render the analysis:
   - In RStudio: open `analysis.Rmd` and click Knit.
   - Or from R: rmarkdown::render("analysis.Rmd")
4. Figures and intermediate outputs are saved into `results/` by the analysis.

Required input files (place under data/)
- counts_OPhio_EXP.txt
- Seq_name_match_Complete_name.csv
- Genotecas_order.csv
- headers_protein_geneBernier_sinOphio.txt
- genes_proteins.txt
- GenomeAnnotation_OphioH327.csv
- My_Tabla_Colonizacion_virulencia_REF_RStudio.csv
- 123_all_PHI_hits.csv
- GO.csv
- 121_gene_coords.txt
- Chr1_8._H327txt.txt
- 123_all_PHI_hits.csv
- for_R_library_correction.csv
 accordingly.)
