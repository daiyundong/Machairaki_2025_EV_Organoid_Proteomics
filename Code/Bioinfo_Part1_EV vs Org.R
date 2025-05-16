# -------------------------------------------------------------------
# Organoid vs EV Combined Proteomics Pipeline (R)
# -------------------------------------------------------------------
# This script imports batch‑1 and batch‑2 proteomics exports, attaches
# clinical metadata, creates unified sample labels, merges the two
# batches, filters proteins
# by missing‑value threshold, performs log2 + median normalisation and a
# differential‑expression test contrasting organoid (Org) versus
# extracellular‑vesicle (EV) samples. Downstream visualisations include a
# Volcano plot and annotation‑rich heatmaps.
# -------------------------------------------------------------------

# ----------------------------
# SECTION 0 – SET‑UP
# ----------------------------
setwd("./Data_used_in_bioinfo") # point to folder containing raw xlsx files
batch1 <- list()   # container for batch‑1 data (cz031)
batch2 <- list()   # container for batch‑2 data (cz039)

# ----------------------------
# SECTION 1 – LIBRARIES
# ----------------------------
library(readxl)
library(dplyr)
library(stringr)
library(tibble)
library(HarmonizR)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(pheatmap)
library(grid)
library(gridExtra)

# ----------------------------
# SECTION 2 – RAW & CLINICAL TABLES
# ----------------------------
# Clinical sample sheets per batch
batch1[["clinical"]] <- read_excel("20240501_cz050_samples list for tymora_updated.xlsx", sheet = "batch 1 (cz031)")
batch2[["clinical"]] <- read_excel("20240501_cz050_samples list for tymora_updated.xlsx", sheet = "batch 2 (cz039)")

# Raw quantitative protein matrices per batch
batch1[["raw"]] <- read_excel("JHMI-Vasso-batch1-organoid_and_EV-quant.xlsx", sheet = "Proteins")
batch2[["raw"]] <- read_excel("JHMI-Vasso-batch2-organoid-and-EV-quant_proteins.xlsx", sheet = "Proteins")

# Additional sample‑level metadata (age, sex, etc.)
sample.clinical <- read_excel("sample_metadata.xlsx")

# ----------------------------
# SECTION 3 – SAMPLE LABELS & METADATA JOIN
# ----------------------------
sample.clinical$CellLineID                     <- tolower(sample.clinical$CellLineID)
batch1[["clinical"]]$`original cell line name` <- tolower(batch1[["clinical"]]$`original cell line name`)
batch2[["clinical"]]$`original cell line name` <- tolower(batch2[["clinical"]]$`original cell line name`)

# Helper: create unique labels (diagnosis × matrix × id × treatment)
create_label <- function(clinical.data){
  clinical.data %>%
    mutate(label = case_when(
      diagnosis == "n/a"     ~ "n/a",
      diagnosis == "AD"      ~ paste0("AD.",     ifelse(`sample type` == "organoids", "Org.", "EV."), id, ifelse(treatment == "+EO", ".eo", ".cl")),
      diagnosis == "AD+NPS"  ~ paste0("ADNPS.",  ifelse(`sample type` == "organoids", "Org.", "EV."), id, ifelse(treatment == "+EO", ".eo", ".cl")),
      diagnosis == "healthy" ~ paste0("NCI.",    ifelse(`sample type` == "organoids", "Org.", "EV."), id, ifelse(treatment == "+EO", ".eo", ".cl"))
    ))
}

batch1[["clinical"]] <- create_label(batch1[["clinical"]])
batch2[["clinical"]] <- create_label(batch2[["clinical"]])

# Merge ancillary clinical variables
colnames(sample.clinical)[2] <- "original cell line name"
batch1[["clinical"]] <- merge(batch1[["clinical"]], sample.clinical[, c(2, 4:7)], by = "original cell line name", all.x = TRUE)
batch2[["clinical"]] <- merge(batch2[["clinical"]], sample.clinical[, c(2, 4:7)], by = "original cell line name", all.x = TRUE)

# ----------------------------
# SECTION 4 – COLUMN RENAMING & BASIC CLEAN‑UP
# ----------------------------
batch1[["clean"]] <- batch1[["raw"]][, c(1, 2, 13:ncol(batch1[["raw"]]))]
batch2[["clean"]] <- batch2[["raw"]][, c(1, 2, 13:ncol(batch2[["raw"]]))]

# Replace Thermo Fisher column headers with sample labels
replace_colnames <- function(preclean.data, clinical.data){
  sample_numbers <- str_extract(colnames(preclean.data), "(?<=Abundance: F)\\d+(?=: Sample,)")
  name_map <- clinical.data %>% filter(`sample #` %in% as.numeric(sample_numbers)) %>% select(`sample #`, label) %>% deframe()
  new_cols <- colnames(preclean.data)
  for (i in seq_along(sample_numbers)) {
    s <- sample_numbers[i]
    if (!is.na(s) && s %in% names(name_map)) new_cols[i] <- name_map[[s]]
  }
  new_cols
}
colnames(batch1[["clean"]]) <- replace_colnames(batch1[["clean"]], batch1[["clinical"]])
colnames(batch2[["clean"]]) <- replace_colnames(batch2[["clean"]], batch2[["clinical"]])

colnames(batch1[["clean"]])[2] <- "Gene"
colnames(batch2[["clean"]])[2] <- "Gene"

# Re‑order to keep Gene column before first abundance column
batch1[["clean"]] <- batch1[["clean"]][, c(1, 2, 75, 3:74)]
batch2[["clean"]] <- batch2[["clean"]][, c(1, 2, 99, 3:98)]

# ----------------------------
# SECTION 5 – MERGE BATCHES (NO ADJUSTMENT)
# ----------------------------
merge_noadj <- list()
merge_noadj[["clean"]] <- merge(batch1[["clean"]], batch2[["clean"]][, c(1, 4:ncol(batch2[["clean"]]))], by = "Accession", all = FALSE)

# ----------------------------
# SECTION 6 – BATCH‑EFFECT CORRECTION (HarmonizR)
# ----------------------------
batch_effect_adjust <- function(clean.data){
  merge_df <- clean.data[, 4:ncol(clean.data)]
  SampleNum <- as.numeric(str_extract(colnames(merge_df), "(?<=\\.(EV|Org)\\.)\\d+"))
  batchType <- ifelse(SampleNum %in% c(1, 2, 3), 1, 2)
  Description <- data.frame(ID = colnames(merge_df), sample = seq_len(ncol(merge_df)), batch = batchType)
  merge_df <- cbind(clean.data[, 1], merge_df)
  write.table(merge_df, "test.tsv", row.names = FALSE, sep = "\t")
  write.csv(Description, "Description.csv", row.names = FALSE)
  tmp <- HarmonizR::harmonizR("test.tsv", "Description.csv")
  tmp$Accession <- rownames(tmp)
  adjust_df <- merge(clean.data[, 1:3], tmp, by = "Accession")
  adjust_df
}

merge_adj <- list()
merge_adj[["clean"]] <- batch_effect_adjust(merge_noadj[["clean"]])
file.remove("./test.tsv", "./Description.csv", "./cured_data.tsv")

# ----------------------------
# SECTION 7 – MISSING‑VALUE FILTER
# ----------------------------
filter_na_rows_improve <- function(my.data, na_threshold = 0.8){
  ad.cl  <- grep('(AD|ADNPS)\\.[A-Za-z]+\\.\\d+\\.cl', colnames(my.data), value = TRUE)
  nci.cl <- grep('NCI\\.[A-Za-z]+\\.\\d+\\.cl', colnames(my.data), value = TRUE)
  ad.eo  <- grep('(AD|ADNPS)\\.[A-Za-z]+\\.\\d+\\.eo', colnames(my.data), value = TRUE)
  nci.eo <- grep('NCI\\.[A-Za-z]+\\.\\d+\\.eo', colnames(my.data), value = TRUE)
  frac <- function(cols) apply(my.data[, cols, drop = FALSE], 1, function(x) mean(x != 0 & !is.na(x)))
  keep <- frac(ad.cl) > na_threshold | frac(nci.cl) > na_threshold | frac(ad.eo) > na_threshold | frac(nci.eo) > na_threshold
  my.data[keep, , drop = FALSE]
}
merge_noadj[["filtermiss"]] <- filter_na_rows_improve(merge_noadj[["clean"]])
merge_adj[["filtermiss"]]   <- filter_na_rows_improve(merge_adj[["clean"]])

# ----------------------------
# SECTION 8 – LOG2 + MEDIAN NORMALISATION
# ----------------------------
norm_transform <- function(x){
  meta <- x[, 1:3]
  q    <- x[, -c(1:3)]
  log2q <- log2(q)                 # raw intensities already offset‑free
  normq <- sweep(log2q, 2, apply(log2q, 2, median, na.rm = TRUE), "-")
  cbind(meta, normq)
}
merge_noadj[["norm"]] <- norm_transform(merge_noadj[["filtermiss"]])
merge_adj[["norm"]]   <- norm_transform(merge_adj[["filtermiss"]])

# ----------------------------
# SECTION 9 – DIFFERENTIAL EXPRESSION (Org vs EV)
# ----------------------------
perform_t_test <- function(my.data){
  expr <- my.data[, sapply(my.data, is.numeric), drop = FALSE]
  diag_vec <- sapply(strsplit(colnames(expr), "\\."), `[`, 2)  # EV vs Org indicator
  g1 <- unique(diag_vec)[1]; g2 <- unique(diag_vec)[2]
  res <- data.frame(Gene = my.data$Gene, Accession = my.data$Accession, Modifications = my.data$Modifications, p_value = NA, log2FC = NA)
  for (i in seq_len(nrow(expr))){
    v1 <- as.numeric(expr[i, diag_vec == g1]); v2 <- as.numeric(expr[i, diag_vec == g2])
    if (sum(!is.na(v1)) > 1 && sum(!is.na(v2)) > 1){
      t <- t.test(v2, v1)  # log2FC = g2 – g1
      res$p_value[i] <- t$p.value
      res$log2FC[i]  <- mean(v2, na.rm = TRUE) - mean(v1, na.rm = TRUE)
    }
  }
  res$padj <- p.adjust(res$p_value, "BH")
  cbind(res, expr)
}

merge_noadj[["DEAnalysis"]] <- perform_t_test(merge_noadj[["norm"]])
merge_adj[["DEAnalysis"]]   <- perform_t_test(merge_adj[["norm"]])

# Retain DEPs at padj < 0.05 & |log2FC| > 1
merge_noadj[["DEGp005"]] <- merge_noadj[["DEAnalysis"]] %>% filter(padj < 0.05, abs(log2FC) > 1)
merge_adj[["DEGp005"]]   <- merge_adj[["DEAnalysis"]]   %>% filter(padj < 0.05, abs(log2FC) > 1)

# ----------------------------
# SECTION 10 – VISUALISATIONS
# ----------------------------
## 10A  Volcano plot

volcano_plot <- function(deg_data, titlename, p_thre = 1e-100, log2fc_thre = 1){
  temp <- deg_data[,1:6]
  temp$select <- ifelse(temp$padj<p_thre & abs(temp$log2FC)>log2fc_thre,temp$Gene,NA)
  temp$group <- ifelse(temp$log2FC > 1 & temp$padj < 0.05, 
                       "Upregulated", 
                       ifelse(temp$log2FC < -1 & temp$padj < 0.05,
                              "Downregulated",
                              "No Difference"))
  temp <- temp[!is.na(temp$padj),]
  
  ggplot(
    temp, aes(x = log2FC, y = -log10(padj), color = factor(group))
  ) + 
    geom_point(size=1.5) +
    geom_label_repel(
      color = "white",
      arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"),
      point.padding = NA,
      box.padding = 0.1,
      aes(label = select, fill = group), 
      size = 2,
      fontface = "bold",
      max.overlaps=30
    ) +
    scale_fill_manual(
      values = c("Upregulated" = "tomato", "Downregulated" = "skyblue", "No Difference" = "grey"),labels=NULL
    ) +
    scale_color_manual(
      values = c("Upregulated" = "tomato", "Downregulated" = "skyblue", "No Difference" = "grey")
    ) +
    theme_bw() +
    labs(x="log2(fold change)", y="-log10 (padj)", title = titlename) +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      legend.title = element_blank() 
    ) +
    geom_hline(yintercept = -log10(0.05),linetype=2,cex=0.5,color = "grey")+
    geom_vline(xintercept = c(-1,1),linetype=2,cex=0.5,color = "grey")
}
volcano_plot(merge_adj[["DEAnalysis"]], titlename = "Org vs EV")
ggsave("Figure2F_Volcano_Org vs EV.pdf",width = 6, height = 6)

## 10B  Heatmap
### Version 1
heatmap_plot <- function(my.data,
                         show_protein_label = FALSE,
                         mytitle,
                         fontsize_col = 6,  
                         angle_col    = 45) { 
  numeric_data <- my.data[, sapply(my.data, is.numeric), drop = FALSE]
  sample_names <- colnames(numeric_data)
  splitted     <- strsplit(sample_names, "\\.")
  diag_vec <- sapply(splitted, function(x) x[2])
  group <- diag_vec
  group_factor <- factor(group)

  annotation_df <- data.frame(Group = group_factor)
  rownames(annotation_df) <- sample_names
  unique_groups <- levels(group_factor)
  group_colors  <- c("EV"="#E9C46A","Org"="#299D8F")
  annotation_colors <- list(Group = c("EV" = "#E9C46A", "Org" = "#299D8F"))
  
  
  pheatmap(
    numeric_data[,4:ncol(numeric_data)],
    scale             = "row",
    show_rownames     = show_protein_label,
    show_colnames     = TRUE,
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    annotation_col    = annotation_df,
    annotation_colors = annotation_colors,
    main              = mytitle,
    fontsize_col      = fontsize_col,
    angle_col         = angle_col
  )
}
dev.off()
pdf("Figure2E_Heatmap1_Org vs EV.pdf", width = 25, height = 6)
heatmap_plot(merge_adj[["DEGp005"]][,c(1:3,7:ncol(merge_adj[["DEGp005"]]))],mytitle="Heatmap - After Adjustment")
dev.off()

### Version 2
heatmap_plot2 <- function(my.data,
                          mytitle,
                          show_protein_label = FALSE,
                          annotation_colors  = NULL) {
  numeric_data <- my.data[, sapply(my.data, is.numeric), drop = FALSE]
  sample_names <- colnames(numeric_data)
  splitted     <- strsplit(sample_names, "\\.")
  
  diag_vec <- sapply(splitted, function(x) x[1])
  type_vec <- sapply(splitted, function(x) x[2])
  treat_vec <- sapply(splitted, function(x) if (length(x) >= 4) x[4] else NA)
  batch_vec <- sapply(splitted, function(x) {
    if (length(x) < 3) {
      "batch2"
    } else {
      num_part <- suppressWarnings(as.integer(x[3]))
      if (!is.na(num_part) && num_part %in% c(1, 2, 3)) {
        "batch1"
      } else {
        "batch2"
      }
    }
  })

  annotation_df <- data.frame(type = factor(type_vec),
                              batch = factor(batch_vec),
                              treat = factor(treat_vec), 
                              diag = factor(diag_vec),
                              stringsAsFactors = FALSE)
  rownames(annotation_df) <- sample_names
  
  annotation_colors <- list(diag = c("AD"= "#FF9999","ADNPS" = "#8B0000","NCI"= "skyblue"),
                            treat = c("cl" = "#FCCDE5", "eo"="#B3DE69"),
                            batch = c("batch1" = "#66c2a6", "batch2" = "#8ea0cb"),
                            type = c("EV"="#E9C46A","Org"="#299D8F"))

  pheatmap(
    mat               = numeric_data,
    scale             = "row",
    show_rownames     = F,
    show_colnames     = F,
    cluster_rows      = TRUE,
    cluster_cols      = T,
    annotation_col    = annotation_df,
    annotation_colors = annotation_colors,
    main              = mytitle
  )
}
p <- heatmap_plot2(merge_adj[["DEGp005"]][,c(1:3,7:ncol(merge_adj[["DEGp005"]]))],mytitle="Heatmap - After Adjustment - with EV vs Org DEPs")
g <- p$gtable
g <- gtable::gtable_add_cols(g, unit(2, "cm"))
pdf("Figure2E_Heatmap2_Org vs EV.pdf", width = 6, height = 6)
grid::grid.draw(g)
dev.off()

