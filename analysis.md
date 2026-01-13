# RNA-Seq.Dynamics of Ophiostoma novo ulmi transcriptome
Dynamics of Ophiostoma novo-ulmi transcriptome during colonizacion of resistant and susceptible Ulmus minor hosts

Packages

```{r}
# Core
library(tidyverse)
library(DESeq2)

# Enrichment
library(clusterProfiler)
library(enrichplot)

# Visualization
library(ggrepel)
library(ggforce)
library(RColorBrewer)
library(cowplot)
library(EnhancedVolcano)
library(pheatmap)
library(ComplexHeatmap)
library(igraph)
library(ggraph)
library(patchwork)
library(ggVennDiagram)
library(ComplexUpset)
library(UpSetR)

# Stats
library(multcomp)
library(multcompView)
library(purrr)

# Utilities
library(matrixStats)
library(openxlsx)
```

Annotation and other files required for parsing

```{r}
##protein_ID-gene_ID from reference Annotation Excel.
id_gene_protein_Bernier=read.table("data/headers_protein_geneBernier_sinOphio.txt", col.names =c("protein_ID", "gene_ID")) ##from data folder (the same for the rest)

##protein_ID-gene_ID from reference .gff3. Doesn't match with the reference Annotation .xlsx. 
geneID_proteinID_mia=read.table("data/genes_proteins.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, colClasses = c("character", "numeric"), col.names = c("gene_ID_gff3", "protein_ID_gff3"))

# Reference Annotation Excell:
Anotacion_Bernier <- read.csv("data/GenomeAnnotation_OphioH327.csv",sep = ";", header = FALSE, stringsAsFactors = FALSE)

# 2nd colum as colnames 
colnames(Anotacion_Bernier) <- Anotacion_Bernier[2, ]

# Remove first 2 rows (headers)
Anotacion_Bernier <- Anotacion_Bernier[-c(1,2), ]

##reformat gene names
Anotacion_Bernier$Protein <- gsub(OphioH327gp", "g", Anotacion_Bernier$Protein)

##add id protein:
final_Anotacion=cbind(id_gene_protein_Bernier,geneID_proteinID_mia,Anotacion_Bernier)

##select interested colnames
final_Anotacion_interested=final_Anotacion[,c(1:10,17:30)]

##BroadFunctions related with Colonization / Previously reported in transcriptomic studies:
Colonizacion_Virulencia_REF=read.csv("data/My_Tabla_Colonizacion_virulencia_REF_RStudio.csv", sep=";", header=TRUE)

##add gene ID for left_join:
Colonizacion_Virulencia_REF <- Colonizacion_Virulencia_REF %>%
  mutate(gene_ID = str_extract(Protein, "gp\\d+") %>% 
           str_replace("gp", "g"))
```

Counts matrix from FeatureCounts (and filtered with Combined Dual-RNA Seq methodology) and metadata

```{r}
## Count matrix from featureCounts
counts <- read.table("data/counts_OPhio_EXP.txt", header = TRUE, row.names = 1, check.names = FALSE)

## Rename (extract "Gen-")
colnames(counts) <- sapply(colnames(counts), function(x) {
  if (grepl("Gen-", x)) {
    sub("^.*(Gen-[0-9]+)_.*$", "\\1", x)
  } else {
    x
  }
})

## relation table
relacion_ID <- read.csv2("data/Seq_name_match_Complete_name.csv")

# new names from column 5
a <- colnames(counts)
a_renamed <- a
a_renamed[6:length(a)] <- relacion_ID$Complete_name[match(a[6:length(a)], relacion_ID$Seq_name)]

# Keep the original name if doesn't match
a_renamed[is.na(a_renamed)] <- a[is.na(a_renamed)]
colnames(counts) <- a_renamed
colnames(counts)

## metadata
metadata <- read.csv("data/Genotecas_order.csv", sep=";")

##Local
metadata_local <- metadata %>% filter(Distance == "Local")

## Columns 
counts_local <- counts[, metadata_local$Complete_name]
colnames(counts_local)

##  New vector with good names 
nuevo_orden <- c(
  # ---- VAD2 ----
  "VAD2-29_O_L_6","VAD2-14_O_L_6","VAD2-15_O_L_6","VAD2-44_O_L_6",
  "VAD2-45_O_L_24","VAD2-26_O_L_24","VAD2-28_O_L_24","VAD2-25_O_L_24",
  "VAD2-22_O_L_72","VAD2-17_O_L_72","VAD2-48_O_L_72","VAD2-33_O_L_72",
  "VAD2-47_O_L_144","VAD2-18_O_L_144","VAD2-35_O_L_144","VAD2-23_O_L_144",
    # ---- MDV2.3 ----
    "MDV2.3-21_O_L_6","MDV2.3-18_O_L_6","MDV2.3-33_O_L_6","MDV2.3-49_O_L_6",
  "MDV2.3-17_O_L_24","MDV2.3-45_O_L_24","MDV2.3-13_O_L_24","MDV2.3-30_O_L_24",
  "MDV2.3-44_O_L_72","MDV2.3-28_O_L_72","MDV2.3-32_O_L_72","MDV2.3-15_O_L_72",
  "MDV2.3-26_O_L_144","MDV2.3-19_O_L_144","MDV2.3-50_O_L_144","MDV2.3-46_O_L_144",
  # ---- MDV1 ----
  "MDV1-50_O_L_6","MDV1-15_O_L_6","MDV1-23_O_L_6","MDV1-13_O_L_6",
  "MDV1-44_O_L_24","MDV1-16_O_L_24","MDV1-41_O_L_24","MDV1-33_O_L_24",
  "MDV1-31_O_L_72","MDV1-46_O_L_72","MDV1-18_O_L_72","MDV1-49_O_L_72",
  "MDV1-24_O_L_144","MDV1-9_O_L_144","MDV1-12_O_L_144","MDV1-45_O_L_144"
)

##  Reorder columns from counts_local 
counts_local_ordenado <- counts_local[, nuevo_orden]

## Double check 
stopifnot(identical(colnames(counts_local_ordenado), nuevo_orden))
```

##LOCAL
MODEL DESeq2
Figure 1A. ~Time here together for the PCA Visualization

```{r}
dds_local_time <- DESeqDataSetFromMatrix(
  countData = counts_local, 
  colData = metadata_local,
  design = ~ Time
)

dds_time=DESeq(dds_local_time)
resultsNames(dds_time)

##tranformación logarítima para PCA
vst_local_time <- varianceStabilizingTransformation(dds_time)

##pca
pca_data_time <- plotPCA(vst_local_time, intgroup = c("Time"), returnData = TRUE)
percentVar_Time <- round(100 * attr(pca_data_time, "percentVar"))

##Genotype:
pca_data_time <- pca_data_time %>% 
  mutate(Genotype = sub("-.*", "", name))

##time as factor:
pca_data_time <- pca_data_time %>% 
  mutate(
    Time     = factor(Time, levels = c(6, 24, 72, 144)),   # orden cronológico
    Genotype = factor(Genotype, levels = c("VAD2","MDV2.3","MDV1"))
  )


##PLOT:
 pca_plot_time=ggplot(pca_data_time, aes(PC1, PC2, color =Time, shape = Genotype)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar_Time[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_Time[2], "% variance")) +
  ggtitle("Local PCA Model: y ~ Genotype*Time (vst)") +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm"),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size=5)))

plot(pca_plot_time)

```
Figure 1B.
Under the hypothesis that Time is the main factor separating samples, I will inspect the top PC1 loadings. This may highlight candidate genes associated with the yeast phase (<48 hpi) and genes potentially linked to sporulation (>48 hpi).

```{r}
## 1) VST 
vst_local_time <- varianceStabilizingTransformation(dds_time, blind = TRUE)

## 2) select same genes as plotPCA (ntop=500)
ntop <- 500
mat  <- assay(vst_local_time)                      # genes x muestras
rv   <- rowVars(mat)
sel  <- order(rv, decreasing = TRUE)[seq_len(min(ntop, nrow(mat)))]

## 3) 
pca_same <- prcomp(t(mat[sel, , drop=FALSE]), center = TRUE, scale. = FALSE)

## 4) % var (same as plotPCA)
percentVar_same <- (pca_same$sdev^2) / sum(pca_same$sdev^2) * 100
round(percentVar_same[1:2], 1)

## 5) genes loadings and samples scores
loadings <- as.data.frame(pca_same$rotation) %>%
  tibble::rownames_to_column("gene_id")

scores <- as.data.frame(pca_same$x) %>%
  tibble::rownames_to_column("sample") %>%
  cbind(as.data.frame(colData(dds_time))[.$sample, , drop=FALSE])

##as factor
scores$Time <- factor(scores$Time, levels = c("6", "24", "72", "144")) 

## 6) top PC1 genes
load_PC1 <- loadings %>%
  transmute(gene_id,
            loading_PC1 = PC1,
            contrib_abs = abs(PC1)) %>%
  arrange(desc(contrib_abs))


## 7) TOP "right" (PC1+)  "left" (PC1−)
topN <- 15  
genes_right_PC1 <- load_PC1 %>% filter(loading_PC1 > 0) %>% slice_head(n = topN)
genes_left_PC1  <- load_PC1 %>% filter(loading_PC1 < 0) %>% slice_head(n = topN)

## with Annotation
genes_right_PC1 <- genes_right_PC1 %>% left_join(final_Anotacion_interested, by = c("gene_id" = "gene_ID_gff3"))
genes_left_PC1  <- genes_left_PC1  %>% left_join(final_Anotacion_interested, by = c("gene_id" = "gene_ID_gff3"))

genes_loading_PC1=rbind(genes_right_PC1,genes_left_PC1)
genes_loading_PC1$direction=ifelse(genes_loading_PC1$loading_PC1>0, "+", "-")
genes_top_PC1 <- genes_loading_PC1 %>%
  mutate(gene_ID_loading_graph = paste(gene_ID, `Final Description`, sep = "; "))

genes_top_PC1 <- genes_top_PC1 %>%
  mutate(gene_ID_loading_graph = factor(
    gene_ID_loading_graph,
    levels = genes_top_PC1 %>%
      arrange(desc(direction), desc(abs(loading_PC1))) %>%
      pull(gene_ID_loading_graph)
  ))


##PLOT
loading_PC1_TIME=ggplot(genes_top_PC1, aes(x = gene_ID_loading_graph, y = loading_PC1, fill = direction)) +
  geom_col() +
  coord_flip() +
  labs(title = "Top genes PC1",
       x = "", y = "PC1 loading") +
  theme_minimal()

loading_PC1_TIME
```
Figure 1
-Thais C.Oliveira	Hydrolase	OphioH327gp6791
-Thais C. Oliveria		        OphioH327gp6502
-L.Bernier oral presentation (unpublished data)	Cell-surface-protein	OphioH327gp5966
-Comeau; Martha Nigg; Sbaraini; Thais C.Oliveira; Anna Fijarczyk	SM_BackBone	OphioH327gp5456
-Thais C. Oliveria		OphioH327gp5230
-Comeau; Martha Nigg; Sbaraini; Thais C.Oliveira; Anna Fijarczyk	CYP450	OphioH327gp5008
-L.Bernier oral presentation (unpublished data)	Permease	OphioH327gp1859
-L.Bernier oral presentation (unpublished data)	Permease	OphioH327gp1857

~Time by Genotype checked. Model for the analysis

```{r}
##time genotype as factor
metadata_local <- metadata_local %>% 
  mutate(
    Time     = factor(Time,     levels = c("6", "24", "72", "144")),
    Genotype = factor(Genotype, levels = c("VAD2", "MDV2.3", "MDV1"))
  )

dds_by_genotype <- list()

# genotypes
genotypes <- unique(metadata_local$Genotype)

# Loop
for (g in genotypes) {
  metadata_sub <- metadata_local %>% filter(Genotype == g)
  counts_sub <- counts_local[, metadata_sub$Complete_name]
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = metadata_sub,
    design = ~ Time
  )
  dds <- DESeq(dds)
  dds_by_genotype[[g]] <- dds
}

##vst
vst_by_genotype <- list()
for (g in names(dds_by_genotype)) {
  vst_by_genotype[[g]] <- varianceStabilizingTransformation (dds_by_genotype[[g]], blind = TRUE)
}

##check
metadata_sub_VAD2=metadata_local %>% filter(Genotype == "VAD2")
metadata_sub_VAD2_24=metadata_sub_VAD2 %>% filter(Time == 24)
```
Results ~Time (by Genotype)

```{r}
results_by_genotype <- list()

for (g in genotypes) {
  dds <- dds_by_genotype[[g]]
  
  res_24 <- results(dds, contrast = c("Time", "24", "6"))
  res_72 <- results(dds, contrast = c("Time", "72", "6"))
  res_144 <- results(dds, contrast = c("Time", "144", "6"))
  
  # Save results
  results_by_genotype[[g]] <- list(
    res_24 = res_24,
    res_72 = res_72,
    res_144 = res_144
  )
}
```

Results for the ~Time model (by Genotype). Filtering criteria: |log2FC| > 2 and adjusted p-value < 0.05. 

```{r}
##UP:
sig_genes_up <- list()

for (g in genotypes) {
  res_list <- results_by_genotype[[g]]
  
  sig_24_up <- as.data.frame(res_list$res_24) %>%
    filter(padj < 0.05, log2FoldChange > 2) %>%
    mutate(gene_ID_gff3 = rownames(.))
  
  sig_72_up <- as.data.frame(res_list$res_72) %>%
    filter(padj < 0.05, log2FoldChange > 2) %>%
    mutate(gene_ID_gff3 = rownames(.))
  
  sig_144_up <- as.data.frame(res_list$res_144) %>%
    filter(padj < 0.05, log2FoldChange > 2) %>%
    mutate(gene_ID_gff3 = rownames(.))
  
  sig_genes_up[[g]] <- list(
    sig_24_up = sig_24_up,
    sig_72_up = sig_72_up,
    sig_144_up = sig_144_up
  )
}

##DOWN
sig_genes_down <- list()

for (g in genotypes) {
  res_list <- results_by_genotype[[g]]
  
  sig_24_down <- as.data.frame(res_list$res_24) %>%
    filter(padj < 0.05, log2FoldChange < -2) %>%
    mutate(gene_ID_gff3 = rownames(.))
  
  sig_72_down <- as.data.frame(res_list$res_72) %>%
    filter(padj < 0.05, log2FoldChange < -2) %>%
    mutate(gene_ID_gff3 = rownames(.))
  
  sig_144_down <- as.data.frame(res_list$res_144) %>%
    filter(padj < 0.05, log2FoldChange < -2) %>%
    mutate(gene_ID_gff3 = rownames(.))
  
  sig_genes_down[[g]] <- list(
    sig_24_down = sig_24_down,
    sig_72_down = sig_72_down,
    sig_144_down = sig_144_down
  )
}
```
From list to dataframe

```{r}
##up

##VAD2
Time_sig_up_VAD2_24=sig_genes_up$VAD2$sig_24_up
Time_sig_up_VAD2_72=sig_genes_up$VAD2$sig_72_up
Time_sig_up_VAD2_144=sig_genes_up$VAD2$sig_144_up
##MDV2.3
Time_sig_up_MDV2.3_24=sig_genes_up$MDV2.3$sig_24_up
Time_sig_up_MDV2.3_72=sig_genes_up$MDV2.3$sig_72_up
Time_sig_up_MDV2.3_144=sig_genes_up$MDV2.3$sig_144_up
##MDV1
Time_sig_up_MDV1_24=sig_genes_up$MDV1$sig_24_up
Time_sig_up_MDV1_72=sig_genes_up$MDV1$sig_72_up
Time_sig_up_MDV1_144=sig_genes_up$MDV1$sig_144_up

##down
##VAD2
Time_sig_down_VAD2_24=sig_genes_down$VAD2$sig_24_down
Time_sig_down_VAD2_72=sig_genes_down$VAD2$sig_72_down
Time_sig_down_VAD2_144=sig_genes_down$VAD2$sig_144_down
##MDV2.3
Time_sig_down_MDV2.3_24=sig_genes_down$MDV2.3$sig_24_down
Time_sig_down_MDV2.3_72=sig_genes_down$MDV2.3$sig_72_down
Time_sig_down_MDV2.3_144=sig_genes_down$MDV2.3$sig_144_down
##MDV1
Time_sig_down_MDV1_24=sig_genes_down$MDV1$sig_24_down
Time_sig_down_MDV1_72=sig_genes_down$MDV1$sig_72_down
Time_sig_down_MDV1_144=sig_genes_down$MDV1$sig_144_down


##all with annotation.
Time_sig_up_VAD2_24_annotated <- Time_sig_up_VAD2_24 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")
Time_sig_up_VAD2_72_annotated<- Time_sig_up_VAD2_72 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")
Time_sig_up_VAD2_144_annotated<- Time_sig_up_VAD2_144 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")

Time_sig_up_MDV2.3_24_annotated <- Time_sig_up_MDV2.3_24 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")
Time_sig_up_MDV2.3_72_annotated<- Time_sig_up_MDV2.3_72 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")
Time_sig_up_MDV2.3_144_annotated<- Time_sig_up_MDV2.3_144 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")

Time_sig_up_MDV1_24_annotated <- Time_sig_up_MDV1_24 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")
Time_sig_up_MDV1_72_annotated<- Time_sig_up_MDV1_72 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")
Time_sig_up_MDV1_144_annotated<- Time_sig_up_MDV1_144 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")

##down
Time_sig_down_VAD2_24_annotated <- Time_sig_down_VAD2_24 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")
Time_sig_down_VAD2_72_annotated<- Time_sig_down_VAD2_72 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")
Time_sig_down_VAD2_144_annotated<- Time_sig_down_VAD2_144 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")

Time_sig_down_MDV2.3_24_annotated <- Time_sig_down_MDV2.3_24 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")
Time_sig_down_MDV2.3_72_annotated<- Time_sig_down_MDV2.3_72 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")
Time_sig_down_MDV2.3_144_annotated<- Time_sig_down_MDV2.3_144 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")

Time_sig_down_MDV1_24_annotated <- Time_sig_down_MDV1_24 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")
Time_sig_down_MDV1_72_annotated<- Time_sig_down_MDV1_72 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")
Time_sig_down_MDV1_144_annotated<- Time_sig_down_MDV1_144 %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")
```
BOXPLOTS UP/DOWN 

```{r}
boxplots_list <- list()

for (g in names(results_by_genotype)) {  
  res_list <- results_by_genotype[[g]]
  
  for (t in names(res_list)) {  
    res_df <- as.data.frame(res_list[[t]])
    
    sobreexpresados <- res_df %>% filter(padj < 0.05, log2FoldChange > 2)
    subexpresados <- res_df %>% filter(padj < 0.05, log2FoldChange < -2)
    
    counts_df <- data.frame(
      Categoria = c("Upregulated", "Downregulated"),
      NumGenes = c(nrow(sobreexpresados), nrow(subexpresados))
    )
    
    plot_title <- paste("Nº DEGs -", g, "-", t)
    
    p <- ggplot(counts_df, aes(x = Categoria, y = NumGenes, fill = Categoria)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("#4575B4", "#D73027")) +
      labs(title = plot_title,
           x = "Direction",
           y = "Number of genes") +
      theme_minimal()
    
    # Save plot
    plot_name <- paste0(g, "_", t)
    boxplots_list[[plot_name]] <- p
    
  }
}

##COMBINED BY GENOTYPE:
##VAD2:
my_colors <- c("Upregulated" = "#D73027", "Downregulated" = "#4575B4")
title_theme <- theme_minimal(base_size = 10)

boxplot_24_VAD2=boxplots_list$VAD2_res_24+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(title = "Time 24 vs 6 VAD2", x = NULL, y = "") +
  title_theme +
   scale_y_continuous(limits = c(0, 70)) + 
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    plot.title = element_text(size = 15)
  )

boxplot_72_VAD2=boxplots_list$VAD2_res_72+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(title = "Time 72 vs 6 VAD2", x = NULL, y = "") +
  scale_y_continuous(limits = c(0, 70)) + 
  title_theme +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    plot.title = element_text(size = 15)
  )

boxplot_144_VAD2=boxplots_list$VAD2_res_144+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 70)) + 
  labs(title = "Time 144 vs 6 VAD2", x = "", y = "") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_text(size = 15, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.key.size = unit(0.8, "cm"),
    plot.title = element_text(size = 15)
  )

# legend 
legend <- get_legend(boxplot_144_VAD2)
boxplot_144_VAD2_nolegend <- boxplot_144_VAD2 + theme(legend.position = "none")

combined_plot_VAD2 <- plot_grid(boxplot_24_VAD2, boxplot_72_VAD2, boxplot_144_VAD2_nolegend,
                           labels = c("A", "C", "E"),
                           ncol = 1, align = "v", rel_heights = c(1, 1, 1))
# Add legend
final_plot_VAD2 <- plot_grid(combined_plot_VAD2, legend, ncol = 1, rel_heights = c(1, 0.1))

# Title
final_plot_labeled_VAD2 <- ggdraw(final_plot_VAD2) +
  # draw_label("Expression direction", x = 0.5, y = 0.03, vjust = 0, size = 15) +
draw_label("Number of genes", x = 0, y = 0.5, angle = 90, vjust = 1, size = 15)


##MDV23
  boxplot_24_MDV2.3=boxplots_list$MDV2.3_res_24+
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    labs(title = "Time 24 vs 6 MDV2.3", x = NULL, y = "") +
    title_theme +
     scale_y_continuous(limits = c(0, 70)) + 
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 15, color = "black"),
      axis.ticks.y = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      plot.title = element_text(size = 15)
    )
  boxplot_72_MDV2.3=boxplots_list$MDV2.3_res_72+
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    labs(title = "Time 72 vs 6 MDV2.3", x = NULL, y = "") +
    scale_y_continuous(limits = c(0, 70)) + 
    title_theme +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 15, color = "black"),
      axis.ticks.y = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      plot.title = element_text(size = 15)
    )
  
  boxplot_144_MDV2.3=boxplots_list$MDV2.3_res_144+
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(limits = c(0, 70)) + 
    labs(title = "Time 144 vs 6 MDV2.3", x = "", y = "") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text = element_text(size = 15, color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      legend.key.size = unit(0.8, "cm"),
      plot.title = element_text(size = 15)
    )

  #legend
  legend <- get_legend(boxplot_144_MDV2.3)
  boxplot_144_MDV2.3_nolegend <- boxplot_144_MDV2.3 + theme(legend.position = "none")

  combined_plot_MDV2.3 <- plot_grid(boxplot_24_MDV2.3, boxplot_72_MDV2.3, boxplot_144_MDV2.3_nolegend,
                             labels = c("A", "C", "E"),
                             ncol = 1, align = "v", rel_heights = c(1, 1, 1))
  #Legend
  final_plot_MDV2.3 <- plot_grid(combined_plot_MDV2.3, legend, ncol = 1, rel_heights = c(1, 0.1))
  
  # Title
  final_plot_labeled_MDV2.3 <- ggdraw(final_plot_MDV2.3) +
    # draw_label("Expression direction", x = 0.5, y = 0.03, vjust = 0, size = 15) +
  draw_label("Number of genes", x = 0, y = 0.5, angle = 90, vjust = 1, size = 15)
  
  
##MDV1:
boxplot_24_MDV1=boxplots_list$MDV1_res_24+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(title = "Time 24 vs 6 MDV1", x = NULL, y = "") +
  title_theme +
   scale_y_continuous(limits = c(0, 70)) + 
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    plot.title = element_text(size = 15)
  )
boxplot_72_MDV1=boxplots_list$MDV1_res_72+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(title = "Time 72 vs 6 MDV1", x = NULL, y = "") +
  scale_y_continuous(limits = c(0, 70)) + 
  title_theme +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    plot.title = element_text(size = 15)
  )

boxplot_144_MDV1=boxplots_list$MDV1_res_144+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 70)) + 
  labs(title = "Time 144 vs 6 MDV1", x = "", y = "") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_text(size = 15, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.key.size = unit(0.8, "cm"),
    plot.title = element_text(size = 15)
  )

##Combined plot:
#legend
legend <- get_legend(boxplot_144_MDV1)
boxplot_144_MDV1_nolegend <- boxplot_144_MDV1 + theme(legend.position = "none")

combined_plot_MDV1 <- plot_grid(boxplot_24_MDV1, boxplot_72_MDV1, boxplot_144_MDV1_nolegend,
                           labels = c("A", "C", "E"),
                           ncol = 1, align = "v", rel_heights = c(1, 1, 1))
#Legend
final_plot_MDV1 <- plot_grid(combined_plot_MDV1, legend, ncol = 1, rel_heights = c(1, 0.1))
#Title
final_plot_labeled_MDV1 <- ggdraw(final_plot_MDV1) +
  # draw_label("Expression direction", x = 0.5, y = 0.03, vjust = 0, size = 15) +
draw_label("Number of genes", x = 0, y = 0.5, angle = 90, vjust = 1, size = 15)

## results
final_plot_labeled_VAD2
final_plot_labeled_MDV2.3
final_plot_labeled_MDV1
```

VOLCANO PLOTS

```{r}
# New list
results_annotated_by_genotype <- list()
for (g in genotypes) {
  dds_results <- results_by_genotype[[g]]
  annotated_results <- list()
  for (timepoint in names(dds_results)) {
    res <- as.data.frame(dds_results[[timepoint]]) %>%
      tibble::rownames_to_column(var = "gene_ID_gff3") %>%  #ID
      left_join(final_Anotacion_interested, by = "gene_ID_gff3")  #Annotation
    annotated_results[[timepoint]] <- res
  }

  results_annotated_by_genotype[[g]] <- annotated_results
}


volcanoPlot_list = list()

for (g in names(results_annotated_by_genotype)) {
  res_list = results_annotated_by_genotype[[g]]
  for (t in names(res_list)) {  
    res_df <- as.data.frame(res_list[[t]])
    plot_title = paste("Volcano Plot DEGs -", g, "-", t)
    Vp = EnhancedVolcano(res_df,
                         lab = res_df$gene_ID,
                         x = "log2FoldChange",
                         y = "padj",
                         pCutoff = 0.05,
                         FCcutoff = 2,
                         title = NULL,
                         subtitle = NULL,
                         caption = NULL,
                         drawConnectors = TRUE,
                         labSize = 3.5,
                         axisLabSize = 15,
                         titleLabSize = 0,
                         subtitleLabSize = 0,
                         captionLabSize = 0,
                         legendPosition = "none",
                         xlim = c(-8.5, 8.5),
                        ylim = c(-0.6, 13))
    plot_name = paste0(g, "_", t)
    volcanoPlot_list[[plot_name]] <- Vp
  }
}


##VAD2
combined_VP_VAD2 <- plot_grid(volcanoPlot_list$VAD2_res_24, volcanoPlot_list$VAD2_res_72, volcanoPlot_list$VAD2_res_144,
                           labels = c("B", "D", "F"),
                           ncol = 1, align = "v", rel_heights = c(1, 1, 1))

##MDV23
combined_VP_MDV2.3 <- plot_grid(volcanoPlot_list$MDV2.3_res_24, volcanoPlot_list$MDV2.3_res_72, volcanoPlot_list$MDV2.3_res_144,
                           labels = c("B", "D", "F"),
                           ncol = 1, align = "v", rel_heights = c(1, 1, 1))


##MDV1
combined_VP_MDV1 <- plot_grid(volcanoPlot_list$MDV1_res_24, volcanoPlot_list$MDV1_res_72, volcanoPlot_list$MDV1_res_144,
                           labels = c("B", "D", "F"),
                           ncol = 1, align = "v", rel_heights = c(1, 1, 1))
##plots:
combined_VP_VAD2
combined_VP_MDV2.3
combined_VP_MDV1
```

Figure 2C (up) and 2F (down). VennDiagrams:

```{r}

# parameters
times <- c(24, 72, 144)
genotypes <- c("VAD2", "MDV1", "MDV2.3")

# Main Loop
for (time in times) {
  cat("\nProcesando tiempo:", time, "h\n")
  ## Sets
  set_up <- list()
  set_down <- list()
  for (genotype in genotypes) {
    up_var <- get(paste0("Time_sig_up_", genotype, "_", time, "_annotated"))
    down_var <- get(paste0("Time_sig_down_", genotype, "_", time, "_annotated"))
    set_up[[genotype]] <- up_var$gene_ID_gff3
    set_down[[genotype]] <- down_var$gene_ID_gff3
  }
  
  ## mixed genes (LFC>2 in one genotyoe while LFC<-2 in other genotype)

  mixed_genes <- union(
    union(
      union(intersect(set_up[["VAD2"]], set_down[["MDV1"]]),
            intersect(set_down[["VAD2"]], set_up[["MDV1"]])
      ),
      union(intersect(set_up[["VAD2"]], set_down[["MDV2.3"]]),
            intersect(set_down[["VAD2"]], set_up[["MDV2.3"]])
      )
    ),
    union(
      intersect(set_up[["MDV1"]], set_down[["MDV2.3"]]),
      intersect(set_down[["MDV1"]], set_up[["MDV2.3"]])
    )
  )

  ## remove mixed genes
  for (genotype in genotypes) {
    set_up[[genotype]] <- set_up[[genotype]][!set_up[[genotype]] %in% mixed_genes]
    set_down[[genotype]] <- set_down[[genotype]][!set_down[[genotype]] %in% mixed_genes]
  }
  
  ## List
  gene_sets_up <- setNames(
    list(set_up[[genotypes[1]]], set_up[[genotypes[2]]], set_up[[genotypes[3]]]),
    paste0(genotypes, "_", time, " up")
  )
  
  gene_sets_down <- setNames(
    list(set_down[[genotypes[1]]], set_down[[genotypes[2]]], set_down[[genotypes[3]]]),
    paste0(genotypes, "_", time, " down")
  )
  
  ## Commons and unique genes
  common_VAD2_MDV1_up <- intersect(set_up[["VAD2"]], set_up[["MDV1"]])
  common_VAD2_MDV1_down <- intersect(set_down[["VAD2"]], set_down[["MDV1"]])

  common_VAD2_MDV2.3_up <- intersect(set_up[["VAD2"]], set_up[["MDV2.3"]])
  common_VAD2_MDV2.3_down <- intersect(set_down[["VAD2"]], set_down[["MDV2.3"]])

  common_MDV1_MDV2.3_up <- intersect(set_up[["MDV1"]], set_up[["MDV2.3"]])
  common_MDV1_MDV2.3_down <- intersect(set_down[["MDV1"]], set_down[["MDV2.3"]])

  common_all_up <- Reduce(intersect, list(set_up[["VAD2"]], set_up[["MDV1"]], set_up[["MDV2.3"]]))
  common_all_down <- Reduce(intersect, list(set_down[["VAD2"]], set_down[["MDV1"]], set_down[["MDV2.3"]]))

  all_other_genes_VAD2 <- union(union(set_up[["MDV1"]], set_down[["MDV1"]]),
                                union(set_up[["MDV2.3"]], set_down[["MDV2.3"]]))
  unique_VAD2_up <- setdiff(set_up[["VAD2"]], all_other_genes_VAD2)
  unique_VAD2_down <- setdiff(set_down[["VAD2"]], all_other_genes_VAD2)
  unique_VAD2 <- union(unique_VAD2_up, unique_VAD2_down)

  all_other_genes_MDV1 <- union(union(set_up[["VAD2"]], set_down[["VAD2"]]),
                                union(set_up[["MDV2.3"]], set_down[["MDV2.3"]]))
  unique_MDV1_up <- setdiff(set_up[["MDV1"]], all_other_genes_MDV1)
  unique_MDV1_down <- setdiff(set_down[["MDV1"]], all_other_genes_MDV1)
  unique_MDV1 <- union(unique_MDV1_up, unique_MDV1_down)

  all_other_genes_MDV2.3 <- union(union(set_up[["VAD2"]], set_down[["VAD2"]]),
                                  union(set_up[["MDV1"]], set_down[["MDV1"]]))
  unique_MDV2.3_up <- setdiff(set_up[["MDV2.3"]], all_other_genes_MDV2.3)
  unique_MDV2.3_down <- setdiff(set_down[["MDV2.3"]], all_other_genes_MDV2.3)
  unique_MDV2.3 <- union(unique_MDV2.3_up, unique_MDV2.3_down)

  ## Summary table
  common_counts_labeled <- tibble(
    Comparison = c(
      "VAD2 & MDV1 (up)", "VAD2 & MDV1 (down)",
      "VAD2 & MDV2.3 (up)", "VAD2 & MDV2.3 (down)",
      "MDV1 & MDV2.3 (up)", "MDV1 & MDV2.3 (down)",
      "Common All (up)", "Common All (down)",
      "Unique VAD2", "Unique MDV1", "Unique MDV2.3",
      "Unique VAD2 (up)", "Unique VAD2 (down)",
      "Unique MDV1 (up)", "Unique MDV1 (down)",
      "Unique MDV2.3 (up)", "Unique MDV2.3 (down)"
    ),
    Gene_Count = c(
      length(common_VAD2_MDV1_up), length(common_VAD2_MDV1_down),
      length(common_VAD2_MDV2.3_up), length(common_VAD2_MDV2.3_down),
      length(common_MDV1_MDV2.3_up), length(common_MDV1_MDV2.3_down),
      length(common_all_up), length(common_all_down),
      length(unique_VAD2), length(unique_MDV1), length(unique_MDV2.3),
      length(unique_VAD2_up), length(unique_VAD2_down),
      length(unique_MDV1_up), length(unique_MDV1_down),
      length(unique_MDV2.3_up), length(unique_MDV2.3_down)
    ),
    Gene_List = sapply(list(
      common_VAD2_MDV1_up, common_VAD2_MDV1_down,
      common_VAD2_MDV2.3_up, common_VAD2_MDV2.3_down,
      common_MDV1_MDV2.3_up, common_MDV1_MDV2.3_down,
      common_all_up, common_all_down,
      unique_VAD2, unique_MDV1, unique_MDV2.3,
      unique_VAD2_up, unique_VAD2_down,
      unique_MDV1_up, unique_MDV1_down,
      unique_MDV2.3_up, unique_MDV2.3_down
    ), function(x) paste(x, collapse = "; "))
  )
  
  assign(paste0("common_counts_labeled_", time), common_counts_labeled)
  

  #  plot 
p_up <- ggVennDiagram(gene_sets_up, label = "count", edge_size = 1.2) + 
  scale_fill_gradient(low = "white", high = "tomato") +
  scale_color_manual(values = rep("black", 3)) +
  theme_void() + 
  theme(
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 40)
  ) +
  ggtitle("")

# Aspect
p_up$layers[[which(sapply(p_up$layers, function(l) "GeomLabel" %in% class(l$geom)))]]$aes_params$size <- 10
  
  p_down <- ggVennDiagram(gene_sets_down, label = "count", edge_size = 1.2) + 
  scale_fill_gradient(low = "white", high = "skyblue") +
  scale_color_manual(values = rep("black", 3)) +
  theme_void() + 
  theme(
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    legend.position = "none",
    text = element_text(size = 40)
  ) +
  ggtitle("")
  
  p_down$layers[[which(sapply(p_down$layers, function(l) "GeomLabel" %in% class(l$geom)))]]$aes_params$size <- 10
  
  print(p_up)
  print(p_down)

}
```
Figure 2A UP. Shared and unique genes:

```{r}
# Time_sig_up_VAD2_24_annotated
# Time_sig_up_VAD2_72_annotated
# Time_sig_up_VAD2_144_annotated

sig_tables <- list(

  VAD2_24_up=Time_sig_up_VAD2_24_annotated,
  VAD2_72_up=Time_sig_up_VAD2_72_annotated,
  VAD2_144_up=Time_sig_up_VAD2_144_annotated,
  MDV1_24_up=Time_sig_up_MDV1_24_annotated,
  MDV1_72_up=Time_sig_up_MDV1_72_annotated,
  MDV1_144_up=Time_sig_up_MDV1_144_annotated,
  MDV2.3_24_up=Time_sig_up_MDV2.3_24_annotated,
  MDV2.3_72_up=Time_sig_up_MDV2.3_72_annotated,
  MDV2.3_144_up=Time_sig_up_MDV2.3_144_annotated
)

# 2. unique genes
all_genes_ID <- unique(unlist(lapply(sig_tables, function(x) x$gene_ID)))
# 3. base dataframe
gene_presence_df <- data.frame(Gene = all_genes_ID, stringsAsFactors = FALSE)

# 4. presence
for (name in names(sig_tables)) {
  gene_presence_df[[name]] <- gene_presence_df$Gene %in% sig_tables[[name]]$gene_ID
}

# head
head(gene_presence_df)

#logic columns 0|1
gene_presence_binary <- gene_presence_df
gene_presence_binary[ , -1] <- lapply(gene_presence_binary[ , -1], as.integer)
head(gene_presence_binary)

##rename columns
colnames(gene_presence_binary)= c("genes", "VAD2:24h","VAD2:72h","VAD2:144h", "MDV1:24h", "MDV1:72h", "MDV1:144h", "MDV23:24h", "MDV23:72h", "MDV23:144h")

gene_long <- gene_presence_binary %>%
  pivot_longer(-genes, names_to = "Model", values_to = "Present") %>%
  filter(Present == 1)

# times that a gene appears
gene_counts <- gene_long %>%
  group_by(genes) %>%
  summarise(N = n())

# model or shared
gene_class <- gene_counts %>%
  mutate(Category = case_when(
    N == 1 ~ "Model-specific",
    N > 1 ~ "Shared"
  )) %>%
  left_join(gene_long, by = "genes")


gene_class$Model <- factor(gene_class$Model, levels = c("VAD2:24h","VAD2:72h","VAD2:144h", "MDV23:24h", "MDV23:72h", "MDV23:144h", "MDV1:24h", "MDV1:72h", "MDV1:144h"), labels = c(
    "R-VAD2:24h", "R-VAD2:72h", "R-VAD2:144h",
    "R-MDV2.3:24h", "R-MDV2.3:72h", "R-MDV2.3:144h",
    "S-MDV1:24h", "S-MDV1:72h", "S-MDV1:144h"
    ))


gene_class_pct <- gene_class %>%
  dplyr::count(Model, Category) %>%
  group_by(Model) %>%
  mutate(prop = n / sum(n), label = paste0(round(prop * 100), "%"))

##plot
resumen_shared_unique_models_barplot_UP=ggplot(gene_class_pct, aes(x = Model, y = n, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 5
  ) +
  scale_fill_manual(values = c("Model-specific" = "tomato", "Shared" = "steelblue")) +
  labs(
    title = "Gene Count per Model (with % Unique vs Shared) UP",
    x = "",
    y = "Gene Count",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "top"
  )

resumen_shared_unique_models_barplot_UP
```

Figure 2D DOWN. Shared and unique genes

```{r}
# Time_sig_down_VAD2_24_annotated
# Time_sig_down_VAD2_72_annotated
# Time_sig_down_VAD2_144_annotated

sig_tables_down <- list(

  VAD2_24_down=Time_sig_down_VAD2_24_annotated,
  VAD2_72_down=Time_sig_down_VAD2_72_annotated,
  VAD2_144_down=Time_sig_down_VAD2_144_annotated,
  MDV1_24_down=Time_sig_down_MDV1_24_annotated,
  MDV1_72_down=Time_sig_down_MDV1_72_annotated,
  MDV1_144_down=Time_sig_down_MDV1_144_annotated,
  MDV2.3_24_down=Time_sig_down_MDV2.3_24_annotated,
  MDV2.3_72_down=Time_sig_down_MDV2.3_72_annotated,
  MDV2.3_144_down=Time_sig_down_MDV2.3_144_annotated
)


# unique genes
# all_genes_down <- unique(unlist(lapply(sig_tables_down, rownames)))

# 2. unique genes
all_genes_ID_down <- unique(unlist(lapply(sig_tables_down, function(x) x$gene_ID)))

# 3. base dataframe
gene_presence_df_down <- data.frame(Gene = all_genes_ID_down, stringsAsFactors = FALSE)


# 4. presence
for (name in names(sig_tables_down)) {
  gene_presence_df_down[[name]] <- gene_presence_df_down$Gene %in% sig_tables_down[[name]]$gene_ID
}
head(gene_presence_df_down)


# logic columns 0|1
gene_presence_binary_down <- gene_presence_df_down
gene_presence_binary_down[ , -1] <- lapply(gene_presence_binary_down[ , -1], as.integer)
head(gene_presence_binary_down)

##colnames
colnames(gene_presence_binary_down)= c("genes", "VAD2:24h","VAD2:72h","VAD2:144h", "MDV1:24h", "MDV1:72h", "MDV1:144h", "MDV23:24h", "MDV23:72h", "MDV23:144h")

##summary table
gene_long_down <- gene_presence_binary_down %>%
  pivot_longer(-genes, names_to = "Model", values_to = "Present") %>%
  filter(Present == 1)

# nº times gene appears.
gene_counts_down <- gene_long_down %>%
  group_by(genes) %>%
  summarise(N = n()) 

# shared or specific
gene_class_down <- gene_counts_down %>%
  mutate(Category = case_when(
    N == 1 ~ "Group-specific",
    N > 1 ~ "Shared"
  )) %>%
  left_join(gene_long_down, by = "genes")


gene_class_down$Model <- factor(gene_class_down$Model, levels = c("VAD2:24h","VAD2:72h","VAD2:144h", "MDV23:24h", "MDV23:72h", "MDV23:144h", "MDV1:24h", "MDV1:72h", "MDV1:144h"), labels = c(
    "R-VAD2:24h", "R-VAD2:72h", "R-VAD2:144h",
    "R-MDV2.3:24h", "R-MDV2.3:72h", "R-MDV2.3:144h",
    "S-MDV1:24h", "S-MDV1:72h", "S-MDV1:144h"
    ))


gene_class_pct_down <- gene_class_down %>%
  dplyr::count(Model, Category) %>%
  group_by(Model) %>%
  mutate(prop = n / sum(n), label = paste0(round(prop * 100), "%"))

resumen_shared_unique_models_barplot_down=ggplot(gene_class_pct_down, aes(x = Model, y = n, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 5
  ) +
  scale_fill_manual(values = c("Group-specific" = "tomato", "Shared" = "steelblue")) +
  labs(
    title = "Gene Count per Model (with % Unique vs Shared) DOWN",
    x = "",
    y = "Gene Count",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "top"
  )

resumen_shared_unique_models_barplot_down
```

Figure 2B. UpSet (up)

```{r}
gene_lists <- list(
  `R-VAD2:24h`    = Time_sig_up_VAD2_24_annotated$gene_ID_gff3,
  `R-VAD2:72h`    = Time_sig_up_VAD2_72_annotated$gene_ID_gff3,
  `R-VAD2:144h`   = Time_sig_up_VAD2_144_annotated$gene_ID_gff3,
  
  `R-MDV2.3:24h`  = Time_sig_up_MDV2.3_24_annotated$gene_ID_gff3,
  `R-MDV2.3:72h`  = Time_sig_up_MDV2.3_72_annotated$gene_ID_gff3,
  `R-MDV2.3:144h` = Time_sig_up_MDV2.3_144_annotated$gene_ID_gff3,
  
  `S-MDV1:24h`    = Time_sig_up_MDV1_24_annotated$gene_ID_gff3,
  `S-MDV1:72h`    = Time_sig_up_MDV1_72_annotated$gene_ID_gff3,
  `S-MDV1:144h`   = Time_sig_up_MDV1_144_annotated$gene_ID_gff3
)


# Unique genes
all_genes <- unique(unlist(gene_lists))

# Binary columns
binary_df <- data.frame(gene_ID = all_genes)

for (set_name in names(gene_lists)) {
  binary_df[[set_name]] <- binary_df$gene_ID %in% gene_lists[[set_name]]
}

up_check=ComplexUpset::upset(
  binary_df,
  intersect = names(gene_lists),
  name = "Upregulated genes",
  width_ratio = 0.2,
  sort_sets = FALSE,
  keep_empty_groups = FALSE,
  min_size = 1,  # Ajusta según el solape mínimo que quieras visualizar
  base_annotations = list(
    'Intersection size' = intersection_size(fill = "tomato")
  ),
  set_sizes = upset_set_size(
    geom = geom_bar(fill = "tomato")
  )
) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),  # Quitar nombres en X si no quieres saturar
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = NULL)
```
Figure 2E. Upset (down)

```{r}
gene_lists <- list(
  `R-VAD2:24h`    = Time_sig_down_VAD2_24_annotated$gene_ID_gff3,
  `R-VAD2:72h`    = Time_sig_down_VAD2_72_annotated$gene_ID_gff3,
  `R-VAD2:144h`   = Time_sig_down_VAD2_144_annotated$gene_ID_gff3,
  
  `R-MDV2.3:24h`  = Time_sig_down_MDV2.3_24_annotated$gene_ID_gff3,
  `R-MDV2.3:72h`  = Time_sig_down_MDV2.3_72_annotated$gene_ID_gff3,
  `R-MDV2.3:144h` = Time_sig_down_MDV2.3_144_annotated$gene_ID_gff3,
  
  `S-MDV1:24h`    = Time_sig_down_MDV1_24_annotated$gene_ID_gff3,
  `S-MDV1:72h`    = Time_sig_down_MDV1_72_annotated$gene_ID_gff3,
  `S-MDV1:144h`   = Time_sig_down_MDV1_144_annotated$gene_ID_gff3
)

# Unique genes
all_genes <- unique(unlist(gene_lists))

# Binary columns
binary_df <- data.frame(gene_ID = all_genes)

for (set_name in names(gene_lists)) {
  binary_df[[set_name]] <- binary_df$gene_ID %in% gene_lists[[set_name]]
}

down_check=upset(
  binary_df,
  intersect = names(gene_lists),
  name = "Downregulated genes",
  width_ratio = 0.2,
  sort_sets = FALSE,
  keep_empty_groups = FALSE,
  min_size = 1,  # Ajusta según el solape mínimo que quieras visualizar
  base_annotations = list(
    'Intersection size' = intersection_size(fill = "skyblue")
  ),
  set_sizes = upset_set_size(
    geom = geom_bar(fill = "skyblue")
  )
) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),  # Quitar nombres en X si no quieres saturar
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = NULL)
```

Figure 3. Hierarchical clustering
Prepare the data

```{r}
results_long <- map_dfr(names(results_by_genotype), function(g) {
  inner <- results_by_genotype[[g]]
  
  map_dfr(names(inner), function(comp) {
    x <- inner[[comp]]                      # DESeqResults
    df <- as.data.frame(x)                  # a data.frame base
    df <- tibble::rownames_to_column(df, var = "Gene")
    
    df <- dplyr::select(df, Gene, log2FoldChange)  # <- forzamos dplyr::select
    df$Genotype <- g
    df$Contrast <- comp
    df
  })
})

results_wide <- results_long %>%
  tidyr::unite(Name, Genotype, Contrast, sep = "_") %>%
  tidyr::pivot_wider(names_from = Name, values_from = log2FoldChange)

#rename
results_wide$gene_ID_gff3=results_wide$Gene

##add info:
results_wide_annotated=results_wide %>%
  left_join(final_Anotacion_interested, by = "gene_ID_gff3")

##check:
results_wide_annotated[results_wide_annotated$gene_ID == "g1479", ]
colnames(results_wide_annotated)

## fixed order:
# "VAD2:24h"   "MDV23:24h"  "MDV1:24h"   "VAD2:72h"   "MDV23:72h"  "MDV1:72h"   "VAD2:144h"  "MDV23:144h" "MDV1:144h" 

results_annot_order=results_wide_annotated[,c(13,2,8,5,3,9,6,4,10,7)]

##UP all_genes_ID
##DOWN all_genes_ID_down

up_down_pheatmap <- union(all_genes_ID, all_genes_ID_down)
common_ids <- intersect(all_genes_ID, all_genes_ID_down)

length(up_down_pheatmap)
##sacar solo para los DEGs con ##all_genes_ID <- unique(unlist(lapply(sig_tables, function(x) x$gene_ID)))
results_DEGs <- results_annot_order %>%
  filter(gene_ID %in% up_down_pheatmap)

colnames(results_DEGs)=c("gene_ID","R-VAD2:24h", "R-MDV23:24h", "S-MDV1:24h", 
               "R-VAD2:72h","R-MDV23:72h", "S-MDV1:72h", 
               "R-VAD2:144h","R-MDV23:144h", "S-MDV1:144h")

results_DEGs=as.data.frame(results_DEGs)
rownames(results_DEGs)=results_DEGs$gene_ID
results_DEGs=results_DEGs[,-c(1)]


results_DEGs$gene_ID=rownames(results_DEGs)

##Anotate for starts and circles (Transcriptomic references)
results_DEGs_F3 <- results_DEGs %>%
  left_join(final_Anotacion_interested, by = "gene_ID")

##table REFERENCES
results_DEGs_F3_REF <- results_DEGs_F3 %>%         
  left_join(Colonizacion_Virulencia_REF,        
            by = "gene_ID") %>%                
  distinct(gene_ID, .keep_all = TRUE)           
                                               
##PHI-db
unique_hits_PHIdb=read.csv("123_all_PHI_hits.csv", header=TRUE, sep=";")
PHI_121=unique_hits_PHIdb[,c(2,4)]
colnames(PHI_121)=c("gene_ID", "PHI_description")
PHI_121_virulence <- PHI_121 %>% 
  # remove number (ej. “…27.778”)
  mutate(PHI_description = str_remove(PHI_description, "\\d+(\\.\\d+)?$")) %>% 
  # filter
  filter(PHI_description != "unaffected_pathogenicity")

results_DEGs_F3_REF <- results_DEGs_F3_REF %>%
  left_join(PHI_121_virulence, by = "gene_ID")

results_DEGs_F3_REF <- results_DEGs_F3_REF %>% 
  mutate(stars = paste0(
           if_else(!is.na(BroadFunction), "★", ""),
           if_else(!is.na(PHI_description),       "●", "")),
         label_row = paste0(stars, " ", gene_ID, ";", `Final Description`))

##give that names (label_row) to row names
rownames(results_DEGs)=results_DEGs_F3_REF$label_row
```

Number of clusters by Elbow method

```{r}
X <- results_DEGs %>%
  as.data.frame()
##numeric
X_numeric <- X[, sapply(X, is.numeric)]

# to matrix
X_matrix <- as.matrix(X_numeric)
#prepare for kmeans
X_clean <- X_numeric[complete.cases(X_numeric), , drop = FALSE]
X_clean <- X_clean[apply(X_clean, 1, sd) > 0, , drop = FALSE]

# 2) Scale
X_sc <- t(scale(t(X_clean)))

# 3)  WSS (within-cluster sum of squares) for mutiple k values
set.seed(1)
ks  <- 1:15
wss <- map_dbl(ks, ~ kmeans(X_sc, centers = .x, nstart = 50, iter.max = 100)$tot.withinss)
elbow_df <- data.frame(k = ks, wss = wss)

# 4) plot
elbow_clusters=ggplot(elbow_df, aes(k, wss)) +
  geom_line() + geom_point() +
  labs(x = "Number of clusters (k)", y = "Within-Cluster Sum of Squares (WSS)",
       title = "Método del codo para elegir k")
```

Figure 3. Plot

```{r}
PHMAP=pheatmap::pheatmap(
  results_DEGs[,c(1:9)],
  cutree_rows = 4, ##elbow
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("darkblue","lightblue","white","orange","red"))(256),  # azul=down, blanco=neutro, rojo=up
  na_col = "black",      # color para los NA (no DEG)
  fontsize = 8,
  main = "DEGs per Genotype and Time\n★References, ●PHI-db virulence",
  gaps_col = c(3, 6),     # separar tiempos
  cellwidth = 10,
  cellheight = 18,
  filename = ".png")  ##to pdf less problems with the length gene descriptions
```

Figure 4. Prepare the data for the plot
```{r}
##we extract genes coordinates with the gff3 by gene_ID_gff3 (genes_121_coords)
genes_121_coords=read.table("data/121_gene_coords.txt", header=TRUE)

##reformat
genes_121_coords$gene_ID_gff3=genes_121_coords$Gene_ID
##left_join
genes_121_coords=genes_121_coords%>%
  left_join(final_Anotacion_interested, by="gene_ID_gff3")

##chromosomes and lenghts
chromosomes_Ophio=read.table("Chr1_8._H327txt.txt", header = FALSE, col.names = c("seq_name", "seq_length", "offset", "line_bases", "line_width"))

# Calculate cumulative positions
fai_df_Ophio <- chromosomes_Ophio %>%
  mutate(start = lag(cumsum(seq_length), default = 0) + 1,
         end = cumsum(seq_length))  ##va benne

colnames(fai_df_Ophio)[1:2] <- c("seq_name", "seq_length")
```

Figure 4A. Plot

```{r}
# Levels and labels
chr_levels <- c("OphioH327chr_1", "OphioH327chr_2", "OphioH327chr_3", 
                "OphioH327chr_4", "OphioH327chr_5", "OphioH327chr_6",
                "OphioH327chr_7", "OphioH327chr_8")
chr_labels <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")

# background
rects_df <- fai_df_Ophio %>%
  dplyr::select(`seq_name`, `seq_length`) %>%     # ← nombres exactos
  dplyr::rename(Chromosome = `seq_name`,
                ChrLength  = `seq_length`) %>%
  mutate(xmin = 0, xmax = ChrLength)

# levels and labels
rects_df$Chromosome <- factor(rects_df$Chromosome, levels = chr_levels, labels = chr_labels)
genes_121_coords$Chromosome <- factor(genes_121_coords$Chromosome, levels = chr_levels, labels = chr_labels)
table(genes_121_coords$Chromosome)

# plot
p_where <- ggplot() +
  geom_rect(data = rects_df, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey95", color = NA) +
  geom_point(data = genes_121_coords,
             aes(x = Start.x, y = 0.1, color = Chromosome),
             position = position_jitter(height = 0.1), size = 2.5, alpha = 0.85) +
  geom_text_repel(data = genes_121_coords,
                  aes(x = Start.x, y = 0.1, label = gene_ID, color = Chromosome),
                  size = 4.5, max.overlaps = 100,
                  segment.size = 0.2, segment.color = "grey50") +
  facet_grid(rows = vars(Chromosome), scales = "fixed", switch = "y") +
  theme_bw(base_size = 12) +
  theme(
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    strip.placement    = "outside",
    strip.background   = element_rect(fill = NA, color = NA),
    strip.text.y.left  = element_text(angle = 0, hjust = 1, size = 10, face = "bold", margin = margin(0,0,0,0)),
    panel.spacing.y    = unit(0.1, "lines"),
    theme(legend.position = "none"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank()
  ) +
  labs(
    # title = "Differentially Expressed Unique Genes by Chromosome",
    x = "Genomic Position (bp)",
    y = NULL,
    color = "Comparison"
  ) +
  coord_cartesian(ylim = c(-1, 1))

print(p_where)
```

Figure 4B. Prepare the data with GO db from reference GO annotation.

```{r}
#global:
GO_db_2017=read.csv("data/GO.csv", sep=";")

##split:
GO_molecular_function=subset(GO_db_2017, gotermType == "molecular_function")
GO_molecular_function$protein_ID_gff3=GO_molecular_function$X.proteinId
GO_biolical_proces=subset(GO_db_2017, gotermType == "biological_process")
GO_celular_component=subset(GO_db_2017, gotermType == "cellular_component")
##Annotation
final_Anotacion_interested_MF=final_Anotacion_interested%>%
  left_join(GO_molecular_function, by="protein_ID_gff3")
```

Figure 4B. Enrichment plot. MF

```{r}
GO_molecular_function$protein_ID_gff3=GO_molecular_function$protein_ID

# 1. TERM2GENE 
TERM2GENE_mf <- GO_molecular_function %>%
  dplyr::select(goName, protein_ID_gff3) %>%
  dplyr::distinct()

annotated_list <- list(
  VAD2_24   = Time_sig_up_VAD2_24_annotated,
  VAD2_72   = Time_sig_up_VAD2_72_annotated,
  VAD2_144  = Time_sig_up_VAD2_144_annotated,
  MDV2.3_24 = Time_sig_up_MDV2.3_24_annotated,
  MDV2.3_72 = Time_sig_up_MDV2.3_72_annotated,
  MDV2.3_144= Time_sig_up_MDV2.3_144_annotated,
  MDV1_24   = Time_sig_up_MDV1_24_annotated,
  MDV1_72   = Time_sig_up_MDV1_72_annotated,
  MDV1_144  = Time_sig_up_MDV1_144_annotated
)

# 3. Run enricher on all groups
enrichment_results <- lapply(names(annotated_list), function(name) {
  
  genes_interest <- annotated_list[[name]]$protein_ID_gff3
  genes_interest <- intersect(genes_interest, TERM2GENE_mf$protein_ID_gff3)
  
  ego <- enricher(
    gene = genes_interest,
    TERM2GENE = TERM2GENE_mf,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  
  return(ego)
})

# 4. Name results
names(enrichment_results) <- names(annotated_list)
```

Figure 4B. Dotplot by group

```{r}
for (name in names(enrichment_results)) {
  ego <- enrichment_results[[name]]
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {  # Solo si hay resultados
    print(
      dotplot(ego, showCategory = 20) +  
        ggtitle(paste("Dotplot:", name))
    )
  } else {
    message(paste("No enrichment for:", name))
  }
}
```
Figure 4B. Combined dotplot

```{r}
# 1. From ego to dataframe

df_list <- lapply(enrichment_results, function(x) {
  if (!is.null(x) && nrow(as.data.frame(x)) > 0) {
    as.data.frame(x)
  } else {
    NULL
  }
})

# 2. by time
# Asignamos el tiempo a cada conjunto
df_list$VAD2_24$Genotype <- "VAD2"; df_list$VAD2_24$Time <- "24h"
df_list$MDV2.3_24$Genotype <- "MDV2.3"; df_list$MDV2.3_24$Time <- "24h"
df_list$MDV1_24$Genotype <- "MDV1"; df_list$MDV1_24$Time <- "24h"

df_list$VAD2_72$Genotype <- "VAD2"; df_list$VAD2_72$Time <- "72h"
df_list$MDV2.3_72$Genotype <- "MDV2.3"; df_list$MDV2.3_72$Time <- "72h"
df_list$MDV1_72$Genotype <- "MDV1"; df_list$MDV1_72$Time <- "72h"

df_list$VAD2_144$Genotype <- "VAD2"; df_list$VAD2_144$Time <- "144h"
df_list$MDV2.3_144$Genotype <- "MDV2.3"; df_list$MDV2.3_144$Time <- "144h"
df_list$MDV1_144$Genotype <- "MDV1"; df_list$MDV1_144$Time <- "144h"

# 3. Combined by time
grap_24_enrich_combine <- bind_rows(df_list$VAD2_24, df_list$MDV2.3_24, df_list$MDV1_24)
grap_72_enrich_combine <- bind_rows(df_list$VAD2_72, df_list$MDV2.3_72, df_list$MDV1_72)
grap_144_enrich_combine <- bind_rows(df_list$VAD2_144, df_list$MDV2.3_144, df_list$MDV1_144)

# 4.parse ratios
parse_ratio <- function(ratio_str) {
  sapply(strsplit(as.character(ratio_str), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}

# 5. Calculate FoldEnrichment
grap_24_enrich_combine$CountRatio <- parse_ratio(grap_24_enrich_combine$GeneRatio)
grap_24_enrich_combine$BgRatio    <- parse_ratio(grap_24_enrich_combine$BgRatio)
grap_24_enrich_combine$FoldEnrichment <- grap_24_enrich_combine$CountRatio / grap_24_enrich_combine$BgRatio

grap_72_enrich_combine$CountRatio <- parse_ratio(grap_72_enrich_combine$GeneRatio)
grap_72_enrich_combine$BgRatio    <- parse_ratio(grap_72_enrich_combine$BgRatio)
grap_72_enrich_combine$FoldEnrichment <- grap_72_enrich_combine$CountRatio / grap_72_enrich_combine$BgRatio

grap_144_enrich_combine$CountRatio <- parse_ratio(grap_144_enrich_combine$GeneRatio)
grap_144_enrich_combine$BgRatio    <- parse_ratio(grap_144_enrich_combine$BgRatio)
grap_144_enrich_combine$FoldEnrichment <- grap_144_enrich_combine$CountRatio / grap_144_enrich_combine$BgRatio

# 6. Combine
combined_enrichment <- bind_rows(grap_24_enrich_combine,
                                 grap_72_enrich_combine,
                                 grap_144_enrich_combine)

# 7. As factor
combined_enrichment$Time <- factor(combined_enrichment$Time, 
                                   levels = c("24h", "72h", "144h"))

combined_enrichment$Description <- factor(combined_enrichment$Description, 
                                          levels = rev(unique(combined_enrichment$Description)))

combined_enrichment$GO_wrapped <- str_wrap(combined_enrichment$Description, width = 30)

combined_enrichment$log10_padj <- -log10(combined_enrichment$p.adjust)

combined_enrichment <- combined_enrichment %>%
  filter(!is.na(Description))


# Ordered
combined_enrichment$Genotype <- factor(
  combined_enrichment$Genotype,
  levels = c("VAD2", "MDV2.3", "MDV1")
)

combined_enrichment <- combined_enrichment %>%
  mutate(Genotype = factor(Genotype,
                           levels = c("VAD2", "MDV2.3", "MDV1"),
                           ordered = TRUE))

# Plot a dot plot faceted by Time (with Genotype on the x-axis).
dotplot_combined_MF <- ggplot(combined_enrichment, aes(x = Genotype, y = GO_wrapped)) +
  geom_point(
    aes(size = FoldEnrichment, fill = log10_padj),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  facet_wrap(~ Time, ncol = 3) +  # <-- Ahora facet por Time
  scale_fill_gradient(low = "orange", high = "red", name = "-log10(adj. p-value)") +
  scale_size(range = c(5, 15)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 18, face = "bold")   (24vs6, 72vs6, 144vs6)
  ) +
  labs(
    title = "GO Molecular Function Enrichment Comparison by Time",
    x = "",
    y = "",
    size = "FoldEnrichment"
  ) +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )

print(dotplot_combined_MF)
```

Figure 4C. Enrichment plot. BP

```{r}
GO_biolical_proces$protein_ID_gff3=GO_biolical_proces$X.proteinId

# 1. TERM2GENE
TERM2GENE_bp <- GO_biolical_proces %>%
  dplyr::select(goName, protein_ID_gff3) %>%
  dplyr::distinct()

# 2. 
annotated_list <- list(
  VAD2_24   = Time_sig_up_VAD2_24_annotated,
  VAD2_72   = Time_sig_up_VAD2_72_annotated,
  VAD2_144  = Time_sig_up_VAD2_144_annotated,
  MDV2.3_24 = Time_sig_up_MDV2.3_24_annotated,
  MDV2.3_72 = Time_sig_up_MDV2.3_72_annotated,
  MDV2.3_144= Time_sig_up_MDV2.3_144_annotated,
  MDV1_24   = Time_sig_up_MDV1_24_annotated,
  MDV1_72   = Time_sig_up_MDV1_72_annotated,
  MDV1_144  = Time_sig_up_MDV1_144_annotated
)

# 3. Run enricher on all groups
enrichment_results_bp <- lapply(names(annotated_list), function(name) {
  
genes_interest_bp <- annotated_list[[name]]$protein_ID_gff3
genes_interest_bp <- intersect(genes_interest_bp, TERM2GENE_bp$protein_ID_gff3)

  ego_bp <- enricher(
    gene = genes_interest_bp,
    TERM2GENE = TERM2GENE_bp,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  
  return(ego_bp)
})

# 4. Name
names(enrichment_results_bp) <- names(annotated_list)
```

Figure 4C. Dotplot by group

```{r}
for (name in names(enrichment_results_bp)) {
  ego <- enrichment_results_bp[[name]]
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {  
    print(
      dotplot(ego, showCategory = 20) +  
        ggtitle(paste("Dotplot:", name))
    )
    
  } else {
    message(paste("No enrichment for:", name))
  }
}
```

Figure 4C. Combined dotplot

```{r}
# 1.ego to dataframe
df_list_bp <- lapply(enrichment_results_bp, function(x) {
  if (!is.null(x) && nrow(as.data.frame(x)) > 0) {
    as.data.frame(x)
  } else {
    NULL
  }
})

# 2. by time
# Asignamos el tiempo a cada conjunto
df_list_bp$VAD2_24$Genotype <- "VAD2"; df_list_bp$VAD2_24$Time <- "24h"
df_list_bp$MDV2.3_24$Genotype <- "MDV2.3"; df_list_bp$MDV2.3_24$Time <- "24h"
df_list_bp$MDV1_24$Genotype <- "MDV1"; df_list_bp$MDV1_24$Time <- "24h"

df_list_bp$VAD2_72$Genotype <- "VAD2"; df_list_bp$VAD2_72$Time <- "72h"
df_list_bp$MDV2.3_72$Genotype <- "MDV2.3"; df_list_bp$MDV2.3_72$Time <- "72h"
df_list_bp$MDV1_72$Genotype <- "MDV1"; df_list_bp$MDV1_72$Time <- "72h"

df_list_bp$VAD2_144$Genotype <- "VAD2"; df_list_bp$VAD2_144$Time <- "144h"
df_list_bp$MDV2.3_144$Genotype <- "MDV2.3"; df_list_bp$MDV2.3_144$Time <- "144h"
df_list_bp$MDV1_144$Genotype <- "MDV1"; df_list_bp$MDV1_144$Time <- "144h"

# 3. Combine 
grap_24_enrich_combine <- bind_rows(df_list_bp$VAD2_24, df_list_bp$MDV2.3_24, df_list_bp$MDV1_24)
grap_72_enrich_combine <- bind_rows(df_list_bp$VAD2_72, df_list_bp$MDV2.3_72, df_list_bp$MDV1_72)
grap_144_enrich_combine <- bind_rows(df_list_bp$VAD2_144, df_list_bp$MDV2.3_144, df_list_bp$MDV1_144)

# 4. 
parse_ratio <- function(ratio_str) {
  sapply(strsplit(as.character(ratio_str), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}

# 5. FoldEnrichment
grap_24_enrich_combine$CountRatio <- parse_ratio(grap_24_enrich_combine$GeneRatio)
grap_24_enrich_combine$BgRatio    <- parse_ratio(grap_24_enrich_combine$BgRatio)
grap_24_enrich_combine$FoldEnrichment <- grap_24_enrich_combine$CountRatio / grap_24_enrich_combine$BgRatio

grap_72_enrich_combine$CountRatio <- parse_ratio(grap_72_enrich_combine$GeneRatio)
grap_72_enrich_combine$BgRatio    <- parse_ratio(grap_72_enrich_combine$BgRatio)
grap_72_enrich_combine$FoldEnrichment <- grap_72_enrich_combine$CountRatio / grap_72_enrich_combine$BgRatio

grap_144_enrich_combine$CountRatio <- parse_ratio(grap_144_enrich_combine$GeneRatio)
grap_144_enrich_combine$BgRatio    <- parse_ratio(grap_144_enrich_combine$BgRatio)
grap_144_enrich_combine$FoldEnrichment <- grap_144_enrich_combine$CountRatio / grap_144_enrich_combine$BgRatio

# 6. Combined
combined_enrichment_bp <- bind_rows(grap_24_enrich_combine,
                                 grap_72_enrich_combine,
                                 grap_144_enrich_combine)

# 7. factor
combined_enrichment_bp$Time <- factor(combined_enrichment_bp$Time, 
                                   levels = c("24h", "72h", "144h"))

combined_enrichment_bp$Description <- factor(combined_enrichment_bp$Description, 
                                          levels = rev(unique(combined_enrichment_bp$Description)))

combined_enrichment_bp$GO_wrapped <- str_wrap(combined_enrichment_bp$Description, width = 30)

combined_enrichment_bp$log10_padj <- -log10(combined_enrichment_bp$p.adjust)

combined_enrichment_bp <- combined_enrichment_bp %>%
  filter(!is.na(Description))

combined_enrichment_bp <- combined_enrichment_bp %>%
  mutate(Genotype = factor(Genotype,
                           levels = c("VAD2", "MDV2.3", "MDV1"),
                           ordered = TRUE))


# Plot a dot plot faceted by Time (with Genotype on the x-axis)
dotplot_combined_bp <- ggplot(combined_enrichment_bp, aes(x = Genotype, y = GO_wrapped)) +
  geom_point(
    aes(size = FoldEnrichment, fill = log10_padj),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  facet_wrap(~ Time, ncol = 3) +  # <-- Ahora facet por Time
  scale_fill_gradient(low = "orange", high = "red", name = "-log10(adj. p-value)") +
  scale_size(range = c(5, 15)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    # plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    # strip.text = element_text(size = 18, face = "bold")  # tamaño del título del facet (24vs6, 72vs6, 144vs6)
  ) +
  labs(
    title = "GO Molecular Function Enrichment Comparison by Time",
    x = "",
    y = "",
    size = "FoldEnrichment"
  ) +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "cm")
  )

print(dotplot_combined_bp)

```

FIGURE 4BC togheter

```{r}

levels(combined_enrichment$Genotype)
levels(combined_enrichment_bp$Genotype)


##MF:
dotplot_combined <- ggplot(combined_enrichment, aes(x = Genotype, y = GO_wrapped)) +
  geom_point(
    aes(size = FoldEnrichment, fill = log10_padj),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  facet_wrap(~ Time, ncol = 3) +  # <-- Ahora facet por Time
  scale_fill_gradient(low = "orange", high = "red", name = "-log10(adj. p-value)") +
  scale_size(range = c(1, 30)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 18, face = "bold")  (24vs6, 72vs6, 144vs6)
  ) +
  labs(
    title = "",
    x = "",
    y = "",
    size = "FoldEnrichment"
  ) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.7, "cm"),
    # ⬇️ ELIMINAR nombres eje X (VAD2, etc.)
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
  )

print(dotplot_combined)


##BP:
dotplot_combined_bp <- ggplot(combined_enrichment_bp, aes(x = Genotype, y = GO_wrapped)) +
  geom_point(
    aes(size = FoldEnrichment, fill = log10_padj),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  facet_wrap(~ Time, ncol = 3) +  # <-- Ahora facet por Time
  scale_fill_gradient(low = "orange", high = "red", name = "-log10(adj. p-value)") +
  scale_size(range = c(1, 30)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, face="bold"),
    # plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    # strip.text = element_text(size = 18, face = "bold")  facet (24vs6, 72vs6, 144vs6)
  ) +
  labs(
    # title = "GO Molecular Function Enrichment Comparison by Time",
    x = "",
    y = "",
    size = "FoldEnrichment"
  ) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.7, "cm"),
    strip.text.x = element_blank(),        # ningún texto
    strip.background = element_blank()     # sin fondo/gris
  )


dotplot_combined_bp

# Figure 4BC
figure_combined_AB <- (dotplot_combined / dotplot_combined_bp) +
  plot_layout(guides = "collect") +         # << une las leyendas
  plot_annotation(tag_levels = "A")         # << añade etiquetas A/B

print(figure_combined_AB)


##Fine-tune

size_vec <- c(combined_enrichment$FoldEnrichment,
              combined_enrichment_bp$FoldEnrichment)
fill_vec <- c(combined_enrichment$log10_padj,
              combined_enrichment_bp$log10_padj)

size_lim   <- range(size_vec, na.rm = TRUE)
fill_lim   <- range(fill_vec, na.rm = TRUE)

# breaks 
size_breaks <- pretty(size_lim, 4)
fill_breaks <- pretty(fill_lim, 5)

size_lim    <- range(c(combined_enrichment$FoldEnrichment,
                       combined_enrichment_bp$FoldEnrichment), na.rm = TRUE)
size_breaks <- pretty(size_lim, 4)

fill_lim    <- range(c(combined_enrichment$log10_padj,
                       combined_enrichment_bp$log10_padj), na.rm = TRUE)
fill_breaks <- pretty(fill_lim, 5)
# Common scale
size_scale <- scale_size_continuous(
  limits = range(c(combined_enrichment$FoldEnrichment,
                   combined_enrichment_bp$FoldEnrichment), na.rm = TRUE),
  range = c(4, 20),  # Tamaño de puntos
  name = "FoldEnrichment"
)

fill_scale <- scale_fill_gradient(
  limits = range(c(combined_enrichment$log10_padj,
                   combined_enrichment_bp$log10_padj), na.rm = TRUE),
  low = "orange", high = "red",
  name = expression(-log[10]~adj.~italic(p))
)

# Panel MF 
dotplot_MF <- dotplot_combined +
  size_scale + fill_scale +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(10, 10, 0, 10)
  )

# Panel BP 
dotplot_BP <- dotplot_combined_bp +
  size_scale + fill_scale +
  theme(
    legend.position = "right",  
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 15),
    legend.key.size = unit(0.6, "cm"),
    plot.margin = margin(0, 10, 10, 10)
  )

# 3. Combined: MF 60 %, BP 40 %

figure_AB <- (dotplot_MF / dotplot_BP) +
  plot_layout(heights = c(4, 2), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20))

print(figure_AB)
```

~Genotype model based on DEGs
Time 24

```{r}
##DEGs up and down.
ids_unicos_24_all=unique(c(Time_sig_up_VAD2_24_annotated$gene_ID_gff3,Time_sig_up_MDV2.3_24_annotated$gene_ID_gff3,Time_sig_up_MDV1_24_annotated$gene_ID_gff3, Time_sig_down_VAD2_24_annotated$gene_ID_gff3, Time_sig_down_MDV2.3_24_annotated$gene_ID_gff3, Time_sig_down_MDV1_24_annotated$gene_ID_gff3))

##counts:
counts_local_filtered_24 <- counts_local[rownames(counts_local) %in% ids_unicos_24_all, ]
##filtering
counts_local_filtered_24_cols_24 <- grep("24$", colnames(counts_local_filtered_24), value = TRUE)

##
counts_local_filtered_24 <- counts_local_filtered_24[, counts_local_filtered_24_cols_24]


##metadata:
# Filter the metadata to keep the same samples
metadata_24 <- metadata_local[metadata_local$Complete_name %in% counts_local_filtered_24_cols_24, ]

# Order
metadata_24 <- metadata_24[match(colnames(counts_local_filtered_24), metadata_24$Complete_name), ]

# DESeq2
dds_24 <- DESeqDataSetFromMatrix(
  countData = counts_local_filtered_24,
  colData   = metadata_24,
  design    = ~ Genotype
)

dds_24 <- DESeq(dds_24)
vsd_24 <- varianceStabilizingTransformation(dds_24, blind = TRUE)

##PCA
mat_24 <- assay(vsd_24)
pca_24=prcomp(t(mat_24))  ##mat_24 tiene forma genes × muestras. Pero la PCA quiere observar cómo se agrupan las muestras, no los genes.

# % var
percentVar_24 <- (pca_24$sdev^2) / sum(pca_24$sdev^2) * 100
round(percentVar_24[1], 1)
round(percentVar_24[2], 1)

##loadings
loadings_24=as.data.frame(pca_24$rotation)
loadings_24$Gene <- rownames(loadings_24)

##variables (sample replicates)
scores_24 <- as.data.frame(pca_24$x)
##column Genotype
scores_24$Genotype <- colData(dds_24)$Genotype

## top genes loadings:
top_PC1_24=head(order(abs(loadings_24[,1]), decreasing = TRUE),15)
top_PC2_24=head(order(abs(loadings_24[,2]), decreasing=TRUE),15)

##position
genes_PC1_24 <- rownames(loadings_24)[top_PC1_24]
genes_PC2_24 <- rownames(loadings_24)[top_PC2_24]

# top of PC1 y PC2
top_genes_24 <- unique(c(genes_PC1_24, genes_PC2_24))

# # Scale
scale_factor_24 <- 10  
loadings_scaled_24 <- loadings_24
loadings_scaled_24$PC1 <- loadings_24$PC1 * scale_factor_24
loadings_scaled_24$PC2 <- loadings_24$PC2 * scale_factor_24

loadings_scaled_top_24 <- loadings_scaled_24[loadings_scaled_24$Gene %in% top_genes_24, ]

##rename
loadings_scaled_top_24$gene_ID_gff3=loadings_scaled_top_24$Gene
loadings_scaled_top_24_annotated=loadings_scaled_top_24%>%
  left_join(final_Anotacion_interested, by="gene_ID_gff3")

# biplot:
loadings_PCA_24=ggplot(scores_24, aes(x = PC1, y = PC2, color = Genotype)) +
  geom_point(size = 3) +
  geom_segment(data = loadings_scaled_top_24_annotated,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2,"cm")), color = "ivory3") +
   geom_text_repel(data = loadings_scaled_top_24_annotated,
                  aes(x = PC1, y = PC2, label = gene_ID),
                  color = "ivory4", size = 3.5,
                  max.overlaps = Inf,
                  box.padding = 0.3,  
                  point.padding = 0.2) +
  theme_classic() +
labs(
    x = paste0("PC1: ", round(percentVar_24[1], 1), "% variance"),
    y = paste0("PC2: ", round(percentVar_24[2], 1), "% variance"),
    title = "PCA biplot with top 25 contributing genes in 24 hpi"
  )
loadings_PCA_24
```

~Genotype model based on DEGs
Time 72

```{r}
##DEGs up and down.
ids_unicos_72_ALL=unique(c(Time_sig_up_VAD2_72_annotated$gene_ID_gff3,Time_sig_up_MDV2.3_72_annotated$gene_ID_gff3,Time_sig_up_MDV1_72_annotated$gene_ID_gff3, Time_sig_down_VAD2_72_annotated$gene_ID_gff3, Time_sig_down_MDV2.3_72_annotated$gene_ID_gff3, Time_sig_down_MDV1_72_annotated$gene_ID_gff3))
s

##counts:
counts_local_filtered_72 <- counts_local[rownames(counts_local) %in% ids_unicos_72_ALL, ]

##filter
counts_local_filtered_72_cols_72 <- grep("72$", colnames(counts_local_filtered_72), value = TRUE)
counts_local_filtered_72 <- counts_local_filtered_72[, counts_local_filtered_72_cols_72]


##metadata:
metadata_72 <- metadata_local[metadata_local$Complete_name %in% counts_local_filtered_72_cols_72, ]
metadata_72 <- metadata_72[match(colnames(counts_local_filtered_72), metadata_72$Complete_name), ]

#DESeq2
dds_72 <- DESeqDataSetFromMatrix(
  countData = counts_local_filtered_72,
  colData   = metadata_72,
  design    = ~ Genotype
)

dds_72 <- DESeq(dds_72)
vsd_72 <- varianceStabilizingTransformation(dds_72, blind = TRUE)
mat_72 <- assay(vsd_72)
pca_72=prcomp(t(mat_72))  

# % var
percentVar_72 <- (pca_72$sdev^2) / sum(pca_72$sdev^2) * 100
round(percentVar_72[1], 1)
round(percentVar_72[2], 1)

##loadings
loadings_72=as.data.frame(pca_72$rotation)
loadings_72$Gene <- rownames(loadings_72)

##scores
scores_72 <- as.data.frame(pca_72$x)
##añado columna Genotype
scores_72$Genotype <- colData(dds_72)$Genotype

##top genes loadings:
top_PC1_72=head(order(abs(loadings_72[,1]), decreasing = TRUE),15)
top_PC2_72=head(order(abs(loadings_72[,2]), decreasing=TRUE),15)

##positions:
genes_PC1_72 <- rownames(loadings_72)[top_PC1_72]
genes_PC2_72 <- rownames(loadings_72)[top_PC2_72]

# genes top PC1 y PC2
top_genes_72 <- unique(c(genes_PC1_72, genes_PC2_72))

# # Scale
scale_factor_72 <- 30 
loadings_scaled_72 <- loadings_72
loadings_scaled_72$PC1 <- loadings_72$PC1 * scale_factor_72
loadings_scaled_72$PC2 <- loadings_72$PC2 * scale_factor_72

loadings_scaled_top_72 <- loadings_scaled_72[loadings_scaled_72$Gene %in% top_genes_72, ]

##rename
loadings_scaled_top_72$gene_ID_gff3=loadings_scaled_top_72$Gene
loadings_scaled_top_72_annotated=loadings_scaled_top_72%>%
  left_join(final_Anotacion_interested, by="gene_ID_gff3")

# biplot
loadings_PCA_72=ggplot(scores_72, aes(x = PC1, y = PC2, color = Genotype)) +
  geom_point(size = 3) +
  geom_segment(data = loadings_scaled_top_72_annotated,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2,"cm")), color = "ivory3") +
   geom_text_repel(data = loadings_scaled_top_72_annotated,
                  aes(x = PC1, y = PC2, label = gene_ID),
                  color = "ivory4", size = 3.5,
                  max.overlaps = Inf,
                  box.padding = 0.3,  
                  point.padding = 0.2) +
  theme_classic() +
labs(
    x = paste0("PC1: ", round(percentVar_72[1], 1), "% variance"),
    y = paste0("PC2: ", round(percentVar_72[2], 1), "% variance"),
    title = "PCA biplot with top 29 contributing genes in 72 hpi"
  )

loadings_PCA_72

```

~Genotype model based on DEGs
Time 144

```{r}
##DEGs up and down.
ids_unicos_144_ALL=unique(c(Time_sig_up_VAD2_144_annotated$gene_ID_gff3,Time_sig_up_MDV2.3_144_annotated$gene_ID_gff3,Time_sig_up_MDV1_144_annotated$gene_ID_gff3, Time_sig_down_VAD2_144_annotated$gene_ID_gff3, Time_sig_down_MDV2.3_144_annotated$gene_ID_gff3, Time_sig_down_MDV1_144_annotated$gene_ID_gff3))

## counts:
counts_local_filtered_144 <- counts_local[rownames(counts_local) %in% ids_unicos_144_ALL, ]

## filter
counts_local_filtered_144_cols_144 <- grep("144$", colnames(counts_local_filtered_144), value = TRUE)
counts_local_filtered_144 <- counts_local_filtered_144[, counts_local_filtered_144_cols_144]

## metadata:
metadata_144 <- metadata_local[metadata_local$Complete_name %in% counts_local_filtered_144_cols_144, ]

# Check
metadata_144 <- metadata_144[match(colnames(counts_local_filtered_144), metadata_144$Complete_name), ]

# DESeq2
dds_144 <- DESeqDataSetFromMatrix(
  countData = counts_local_filtered_144,
  colData   = metadata_144,
  design    = ~ Genotype
)

dds_144 <- DESeq(dds_144)
vsd_144 <- varianceStabilizingTransformation(dds_144, blind = TRUE)

mat_144 <- assay(vsd_144)
pca_144=prcomp(t(mat_144))  

# var
percentVar_144 <- (pca_144$sdev^2) / sum(pca_144$sdev^2) * 100
round(percentVar_144[1], 1)
round(percentVar_144[2], 1)

## loadings
loadings_144=as.data.frame(pca_144$rotation)
loadings_144$Gene <- rownames(loadings_144)

## scores
scores_144 <- as.data.frame(pca_144$x)
##añado columna Genotype
scores_144$Genotype <- colData(dds_144)$Genotype

## top genes loadings:
top_PC1_144=head(order(abs(loadings_144[,1]), decreasing = TRUE),15)
top_PC2_144=head(order(abs(loadings_144[,2]), decreasing=TRUE),15)

## positions:
genes_PC1_144 <- rownames(loadings_144)[top_PC1_144]
genes_PC2_144 <- rownames(loadings_144)[top_PC2_144]

# genes top PC1 y PC2
top_genes_144 <- unique(c(genes_PC1_144, genes_PC2_144))

## Scale
scale_factor_144 <- 40  # ajusta este valor hasta que las flechas se vean proporcionadas
loadings_scaled_144 <- loadings_144
loadings_scaled_144$PC1 <- loadings_144$PC1 * scale_factor_144
loadings_scaled_144$PC2 <- loadings_144$PC2 * scale_factor_144

loadings_scaled_top_144 <- loadings_scaled_144[loadings_scaled_144$Gene %in% top_genes_144, ]

## rename
loadings_scaled_top_144$gene_ID_gff3=loadings_scaled_top_144$Gene
loadings_scaled_top_144_annotated=loadings_scaled_top_144%>%
  left_join(final_Anotacion_interested, by="gene_ID_gff3")


# biplot
loadings_PCA_144=ggplot(scores_144, aes(x = PC1, y = PC2, color = Genotype)) +
  geom_point(size = 3) +
  geom_segment(data = loadings_scaled_top_144_annotated,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2,"cm")), color = "ivory3") +
   geom_text_repel(data = loadings_scaled_top_144_annotated,
                  aes(x = PC1, y = PC2, label = gene_ID),
                  color = "ivory4", size = 3.5,
                  max.overlaps = Inf, # permite más solapamientos controlados
                  box.padding = 0.3,  # separación alrededor de cada etiqueta
                  point.padding = 0.2) +
  theme_classic() +
labs(
    x = paste0("PC1: ", round(percentVar_144[1], 1), "% variance"),
    y = paste0("PC2: ", round(percentVar_144[2], 1), "% variance"),
    title = "PCA biplot with top 28 contributing genes in 144 hpi"
  )

loadings_PCA_144

```

##DISTAL

Figure 5.

Metdatada and why no DESeq2

```{r}
##Distal
metadata_distal <- metadata %>% filter(Distance == "Distal")
counts_distal <- counts[, metadata_distal$Complete_name]
metadata_distal$Genotype=as.factor(metadata_distal$Genotype)
metadata_distal$Time=as.factor(metadata_distal$Time)
```
~Time:

```{r}
metadata_distal <- metadata_distal %>% 
  mutate(
    Time     = factor(Time,     levels = c("6", "24", "72", "144")),
    Genotype = factor(Genotype, levels = c("VAD2", "MDV2.3", "MDV1"))
  )

dds_distal_time <- DESeqDataSetFromMatrix(
  countData = counts_distal,
  colData = metadata_distal,
  design = ~ Time
)

head(rowSums(counts(dds_distal_time)))  
summary(rowSums(counts(dds_distal_time)))
colnames(dds_distal_time)

sample_info <- data.frame(sample = colnames(dds_distal_time)) %>%
  separate(sample, into = c("Genotipo", "Muestra", "O", "D", "Tiempo"), sep = "[-_]", remove = FALSE)

# factor
sample_info$Tiempo <- factor(sample_info$Tiempo, levels = c("6", "24", "72", "144"))
sample_info$Genotipo=factor(sample_info$Genotipo, levels=c("VAD2","MDV2.3","MDV1"))


# Count matrix
count_matrix <- counts(dds_distal_time)

# dframe before pivot_longer()
str(data.frame(sample_info, t(count_matrix)))

# Create a data frame ensuring that categorical columns are not mixed with the count data
df_combined <- bind_cols(sample_info, as.data.frame(t(count_matrix)))

# Counts summary by genotype and time
count_summary <- df_combined %>%
  pivot_longer(cols = starts_with("gene_"), names_to = "Gene", values_to = "Counts") %>% 
  group_by(Genotipo, Tiempo) %>%
  summarise(Ophiostoma_presence = sum(Counts, na.rm = TRUE), .groups = "drop")

##english
count_summary$Genotype=count_summary$Genotipo

# Ver resultados
print(count_summary)

df_counts <- df_combined %>%
  pivot_longer(starts_with("gene_"),
               names_to  = "Gene",
               values_to = "Counts") %>%
  group_by(sample, Genotype = Genotipo, Tiempo) %>%   # renombramos aquí
  summarise(total_counts = sum(Counts), .groups = "drop")

##add mean and se
df_summary <- df_counts %>%
  group_by(Genotype, Tiempo) %>%
  summarise(mean_counts = mean(total_counts),
            se_counts   = sd(total_counts) / sqrt(n()),
            n           = n(),
            .groups     = "drop")
```

 ~Time by Genotype

```{r}
dds_by_genotype_distal <- list()
genotypes <- unique(metadata_distal$Genotype)
for (g in genotypes) {
  metadata_sub <- metadata_distal %>% filter(Genotype == g)
  counts_sub <- counts_distal[, metadata_sub$Complete_name]  # filter
  keep <- rowSums(counts_sub > 0) >= 2
  counts_sub <- counts_sub[keep, ]
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = metadata_sub,
    design = ~ Time
  )

  # DESeq2
  dds <- DESeq(dds, sfType = "poscounts")

  dds_by_genotype_distal[[g]] <- dds
}

table(keep)

##vst
vst_by_genotype_distal <- list()
for (g in names(dds_by_genotype_distal)) {
  vst_by_genotype_distal[[g]] <- varianceStabilizingTransformation (dds_by_genotype_distal[[g]], blind = TRUE)
}

##DESeq2 does not make sense here
metadata_sub_VAD2=metadata_distal %>% filter(Genotype == "VAD2")
metadata_sub_VAD2_24=metadata_sub_VAD2 %>% filter(Time == 24)
```

Normalize for library size.

```{r}
Tamaño_libreria=read.csv("data/for_R_library_correction.csv", sep=";", dec = ",")

Length_library_distal=Tamaño_libreria%>%
  left_join(metadata_distal, by="Seq_name")

library_distal <- Length_library_distal %>%
  filter(!is.na(Genotype))

#rename
library_distal <- library_distal %>%
  rename(sample = Complete_name)

#left join
df_counts_distal_library=df_counts %>%
  left_join(library_distal, by="sample")

#select columns
df_counts_distal_library=df_counts_distal_library[,c(1,2,3,4,5,6,8,9)]
head(df_counts_distal_library)


df_summary_library <- df_counts_distal_library %>%
  group_by(Genotype.x, Tiempo) %>%
  summarise(mean_counts = mean(total_counts),
            se_counts   = sd(total_counts) / sqrt(n()),
            n           = n(),
            Mseqs       = sum(Mseqs),
            .groups     = "drop")



##########HERE HERE HERE #########################


##¿OFFSET OR COVARIATE?:
library(MASS)
library(emmeans)

m_offset <- glm.nb(total_counts ~ Genotype.x * Tiempo + offset(log(Mseqs)), data = df_counts_distal_library)
m_cov <- glm.nb(total_counts ~ Genotype.x * Tiempo + log(Mseqs), data = df_counts_distal_library)

# wich better?
AIC(m_offset, m_cov)  ## p = 0.0015 < 0.05 the improvement is statistically significant.
anova(m_offset, m_cov, test = "Chisq") ##χ² = 21.28, df = 1, p < 0.00001

##¿Interaction?
##without interacción:
model_final <- glm.nb(total_counts ~ Genotype.x + Tiempo + log(Mseqs), data = df_counts_distal_library)
##con interacción:
model_final2 <- glm.nb(total_counts ~ Genotype.x * Tiempo + log(Mseqs), data = df_counts_distal_library) ##la interacción funciona mejor.
AIC(model_final, model_final2)
anova(model_final, model_final2, test = "Chisq")

##with interaction and log(Mseqs)

##signitfivativo en qué tiempo:

##Asterisks indicate significant differences relative to the resistant genotype (R-VAD2) at the same time point according to a negative binomial generalized ##linear model (GLM) adjusted for library size. 

summary(model_final2) ##here are the results for the figure 5.
```

Figure 5 GRAPH:

```{r}

# Asegurar factores y referencias
df_DISTAL_model <- df_counts_distal_library %>%
  mutate(
    Genotype.x = factor(Genotype.x, levels = c("MDV1","MDV2.3","VAD2")),
    Tiempo     = factor(Tiempo,     levels = c("6","24","72","144")),
    # Evita ceros o negativos en Mseqs si vas a usar offset:
    Mseqs      = as.numeric(Mseqs)
  )

pos=position_dodge(0.9)

# Gráfico con asteriscos
plot_ophiostoma <- ggplot(df_summary,
                          aes(x = factor(Tiempo, levels = sort(unique(Tiempo))),
                              y = mean_counts,
                              fill = Genotype)) +
  geom_col(position = pos, colour = "black") +
  geom_errorbar(aes(ymin = mean_counts - se_counts,
                    ymax = mean_counts + se_counts),
                width = 0.2, position = pos) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "",
    x = "hpi",
    y = "Mean read counts ± SE"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y  = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.key.size = unit(1.1, "cm")
  )

print(plot_ophiostoma)
```

END.
