# install_packages.R
# source("install_packages.R")

required <- c(
  "tidyverse","DESeq2","clusterProfiler","enrichplot","ggrepel","ggforce",
  "RColorBrewer","cowplot","EnhancedVolcano","pheatmap","ComplexHeatmap",
  "igraph","ggraph","patchwork","ggVennDiagram","ComplexUpset","UpSetR",
  "multcomp","multcompView","purrr","matrixStats","openxlsx","here","rmarkdown",
  "MASS","emmeans"
)

to_install <- required[!(required %in% installed.packages()[, "Package"])]
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

# Some Bioconductor packages: DESeq2 and clusterProfiler may need BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_needed <- c("DESeq2","clusterProfiler","enrichplot","ComplexHeatmap")
for (pkg in bioc_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

message("All packages installed (or already present). You can now knit analysis.Rmd.")
