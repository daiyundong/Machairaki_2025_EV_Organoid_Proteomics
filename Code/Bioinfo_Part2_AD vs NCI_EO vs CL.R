# -------------------------------------------------------------------
# Extracellular Vesicle & Organoid Proteomics Pipeline (R)
# -------------------------------------------------------------------
# This script ingests raw quantitative proteomics data, merges them with
# sample‑level metadata, cleans and normalises the expression matrix,
# performs batch‑effect correction, exploratory PCA, differential
# expression testing (disease‑ and treatment‑based contrasts), and a
# suite of downstream visualisations (boxplots, volcano plots, violin
# plots, PCA overlays and annotated heatmaps).
# -------------------------------------------------------------------

# ----------------------------
# SECTION 0 – SET‑UP
# ----------------------------
setwd("./Data_used_in_bioinfo")                       # point to folder containing raw xlsx files
EV  <- list()                         # container for extracellular‑vesicle data
Org <- list()                         # container for organoid data

# ----------------------------
# SECTION 1 – LIBRARIES
# ----------------------------
library(dplyr)
library(readxl)
library(stringr)
library(tibble)
library(HarmonizR)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(ggfortify)
library(pheatmap)
library(purrr)

# ----------------------------
# SECTION 2 – CLINICAL METADATA
# ----------------------------
# Read sample sheets for the two experimental batches
batch1.clinical <- read_excel("20240501_cz050_samples list for tymora_updated.xlsx", sheet = "batch 1 (cz031)")
batch2.clinical <- read_excel("20240501_cz050_samples list for tymora_updated.xlsx", sheet = "batch 2 (cz039)")

# Split into EV / organoid and bind batches together
EV[["clinical"]]  <- rbind(batch1.clinical[batch1.clinical$`sample type` == "EVs",],
                           batch2.clinical[batch2.clinical$`sample type` == "EVs",])
Org[["clinical"]] <- rbind(batch1.clinical[batch1.clinical$`sample type` == "organoids",],
                           batch2.clinical[batch2.clinical$`sample type` == "organoids",])

# ----------------------------
# SECTION 3 – RAW QUANTITATIVE MATRICES
# ----------------------------
EV[["raw"]]  <- read_excel("JHMI-Vasso-batch1_and_2-media-EV-quant.xlsx",  sheet = "Proteins")
Org[["raw"]] <- read_excel("032424-JHMI-Vasso-organoid-batch_1_2-quant.xlsx", sheet = "Proteins")

# ----------------------------
# SECTION 4 – SAMPLE‑LEVEL METADATA & LABELS
# ----------------------------
sample.clinical <- read_excel("sample_metadata.xlsx")

sample.clinical$CellLineID                     <- tolower(sample.clinical$CellLineID)
EV[["clinical"]]$`original cell line name`  <- tolower(EV[["clinical"]]$`original cell line name`)
Org[["clinical"]]$`original cell line name` <- tolower(Org[["clinical"]]$`original cell line name`)

# Function: build unique sample labels incorporating disease, matrix
# (EV/organoid), numeric ID and treatment status
create_label <- function(clinical.data) {
  clinical.data %>%
    mutate(
      label = case_when(
        diagnosis == "n/a"    ~ "n/a",
        diagnosis == "AD"     ~ paste0("AD.",     ifelse(`sample type` == "organoids", "Org.", "EV."), id, ifelse(treatment == "+EO", ".eo", ".cl")),
        diagnosis == "AD+NPS" ~ paste0("ADNPS.",  ifelse(`sample type` == "organoids", "Org.", "EV."), id, ifelse(treatment == "+EO", ".eo", ".cl")),
        diagnosis == "healthy"~ paste0("NCI.",    ifelse(`sample type` == "organoids", "Org.", "EV."), id, ifelse(treatment == "+EO", ".eo", ".cl"))
      )
    )
}
EV[["clinical"]]  <- create_label(EV[["clinical"]])
Org[["clinical"]] <- create_label(Org[["clinical"]])

# Merge extra clinical attributes (sex, age, etc.)
colnames(sample.clinical)[2] <- "original cell line name"
EV[["clinical"]]  <- merge(EV[["clinical"]],  sample.clinical[, c(2, 4:7)], by = "original cell line name", all.x = TRUE)
Org[["clinical"]] <- merge(Org[["clinical"]], sample.clinical[, c(2, 4:7)], by = "original cell line name", all.x = TRUE)

# ----------------------------
# SECTION 5 – COLUMN RENAMING & BASIC CLEAN‑UP
# ----------------------------
EV[["clean"]]  <- EV[["raw"]][,  c(1, 2, 13:ncol(EV[["raw"]]))]
Org[["clean"]] <- Org[["raw"]][, c(1, 2, 13:ncol(Org[["raw"]]))]

# replace ‘Abundances (Normalised): <sample #>’ headers with
# sample labels constructed above, appending ".1", ".2"… where the
# same label appears multiple times (e.g. technical replicates)
replace_colnames <- function(preclean.data, clinical.data) {
  sample_numbers <- str_extract(colnames(preclean.data), "(?<=Abundances \\(Normalized\\): )\\d+")
  name_map <- clinical.data %>%
    filter(`sample #` %in% as.numeric(sample_numbers)) %>%
    select(`sample #`, label) %>%
    deframe()
  new_colnames <- colnames(preclean.data)
  for (i in seq_along(sample_numbers)) {
    sample_num <- sample_numbers[i]
    if (!is.na(sample_num) && sample_num %in% names(name_map)) {
      new_colnames[i] <- name_map[[sample_num]]
    }
  }
  name_counts <- ave(seq_along(new_colnames), new_colnames, FUN = seq)
  ifelse(duplicated(new_colnames) | duplicated(new_colnames, fromLast = TRUE),
         paste0(new_colnames, ".", name_counts),
         new_colnames)
}
colnames(EV[["clean"]])  <- replace_colnames(EV[["clean"]],  EV[["clinical"]])
colnames(Org[["clean"]]) <- replace_colnames(Org[["clean"]], Org[["clinical"]])
colnames(EV[["clean"]])[2]  <- "Gene"
colnames(Org[["clean"]])[2] <- "Gene"

# ----------------------------
# SECTION 6 – BATCH‑EFFECT CORRECTION (HarmonizR)
# ----------------------------
batch_effect_adjust <- function(clean.data) {
  merge_df   <- clean.data[, 4:ncol(clean.data)]
  SampleNum  <- as.numeric(str_extract(colnames(merge_df), "(?<=\\.(EV|Org)\\.)\\d+"))
  batchType  <- ifelse(SampleNum %in% c(1, 2, 3), 1, 2)
  Description <- data.frame(ID = colnames(merge_df), sample = seq_len(ncol(merge_df)), batch = batchType)
  merge_df   <- cbind(clean.data[, 1], merge_df)
  write.table(merge_df, "test.tsv", row.names = FALSE, sep = "\t")
  write.csv(Description, "Description.csv", row.names = FALSE)
  temp  <- HarmonizR::harmonizR("test.tsv", "Description.csv")
  temp$Accession <- rownames(temp)
  merge(clean.data[, 1:3], temp, by = "Accession")
}
EV[["adjust"]]  <- batch_effect_adjust(EV[["clean"]])
Org[["adjust"]] <- batch_effect_adjust(Org[["clean"]])
file.remove("./test.tsv", "./Description.csv", "./cured_data.tsv")

# ----------------------------
# SECTION 7 – EXPLORATORY PCA (PRE/POST BATCH CORRECTION)
# ----------------------------
# Function draws PCA coloured by batch assignment and saves PDF files.
batch_pca_plot <- function(my.data, my.title = "PCA plot") {
  numeric_data <- my.data[, sapply(my.data, is.numeric), drop = FALSE]
  numeric_data <- na.omit(numeric_data)
  sample_names <- colnames(numeric_data)
  splitted     <- strsplit(sample_names, "\\.")
  batch_vec    <- sapply(splitted, function(x) ifelse(length(x) < 3, "batch2", {num <- suppressWarnings(as.integer(x[3])); if (!is.na(num) && num %in% c(1, 2, 3)) "batch1" else "batch2"}))
  pca_res      <- prcomp(t(numeric_data), center = TRUE, scale. = TRUE)
  df_pca       <- as.data.frame(pca_res$x)
  var_expl     <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
  ggplot(data.frame(PC1 = df_pca$PC1, PC2 = df_pca$PC2, Group = factor(batch_vec)),
         aes(PC1, PC2, colour = Group)) +
    geom_point(size = 2) +
    labs(x = paste0("PC1: ", sprintf("%.2f", var_expl[1]), "%"),
         y = paste0("PC2: ", sprintf("%.2f", var_expl[2]), "%"),
         title = my.title) +
    theme_minimal() +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) +
    scale_colour_manual(values = c(batch1 = "#66c2a6", batch2 = "#8ea0cb"))
}
# Generate PCA plots & save
batch_pca_plot(EV[["clean"]],my.title = "PCA - EV - Before BEC")
ggsave("Figure1D_PCA_EV_BeforeBEC.pdf",width = 6, height = 6)
batch_pca_plot(EV[["adjust"]],my.title = "PCA - EV - After BEC")
ggsave("Figure1D_PCA_EV_AfterBEC.pdf",width = 6, height = 6)
batch_pca_plot(Org[["clean"]],my.title = "PCA - Org - Before BEC")
ggsave("Figure1B_PCA_Org_BeforeBEC.pdf",width = 6, height = 6)
batch_pca_plot(Org[["adjust"]],my.title = "PCA - Org - After BEC")
ggsave("Figure1B_PCA_Org_AfterBEC.pdf",width = 6, height = 6)

# ----------------------------
# SECTION 8 – MISSING‑VALUE FILTERING
# ----------------------------
# Retain proteins with ≥ `na_threshold` proportion of non‑zero values in
# any diagnostic × treatment subgroup.
filter_na_rows_improve <- function(my.nooutlier.data, na_threshold = 0.8) {
  ad.cl  <- grep("(AD|ADNPS)\\.[A-Za-z]+\\.\\d+\\.cl", colnames(my.nooutlier.data), value = TRUE)
  nci.cl <- grep("NCI\\.[A-Za-z]+\\.\\d+\\.cl", colnames(my.nooutlier.data), value = TRUE)
  ad.eo  <- grep("(AD|ADNPS)\\.[A-Za-z]+\\.\\d+\\.eo", colnames(my.nooutlier.data), value = TRUE)
  nci.eo <- grep("NCI\\.[A-Za-z]+\\.\\d+\\.eo", colnames(my.nooutlier.data), value = TRUE)
  frac_ok <- function(cols) apply(my.nooutlier.data[, cols, drop = FALSE], 1, function(x) mean(x != 0 & !is.na(x)))
  keep_rows <- frac_ok(ad.cl) > na_threshold | frac_ok(nci.cl) > na_threshold | frac_ok(ad.eo) > na_threshold | frac_ok(nci.eo) > na_threshold
  my.nooutlier.data[keep_rows, , drop = FALSE]
}
EV[["filtermiss"]]  <- filter_na_rows_improve(EV[["adjust"]])
Org[["filtermiss"]] <- filter_na_rows_improve(Org[["adjust"]])

# ----------------------------
# SECTION 9 – LOG2 TRANSFORM & MEDIAN NORMALISATION
# ----------------------------
norm_transform <- function(my.filtermiss.data) {
  meta_cols <- my.filtermiss.data[, 1:3]
  quant     <- my.filtermiss.data[, -c(1:3)]
  log2_data <- log2(quant + 1)
  norm_data <- sweep(log2_data, 2, apply(log2_data, 2, median, na.rm = TRUE), FUN = "-")
  cbind(meta_cols, norm_data)
}
EV[["norm"]]  <- norm_transform(EV[["filtermiss"]])
Org[["norm"]] <- norm_transform(Org[["filtermiss"]])

# ----------------------------
# SECTION 10 – DIFFERENTIAL EXPRESSION (CUSTOM t‑TEST WRAPPER)
# ----------------------------
perform_t_test <- function(my.data,
                           diag = FALSE,
                           treat = FALSE,
                           paired = FALSE,
                           reverse_fc = FALSE,
                           filter_diag  = NULL,
                           filter_treat = NULL) {
  if (!is.data.frame(my.data) && !is.matrix(my.data)) {
    stop("Input data must be a data frame or matrix.")
  }
  my.data <- as.data.frame(my.data)
  
  if (diag && treat) {
    stop("Please choose either diag or treat (but not both) for t-test in this function.")
  }
  
  expr_data <- my.data[, sapply(my.data, is.numeric), drop = FALSE]
  if (nrow(expr_data) == 0 || ncol(expr_data) == 0) {
    stop("No numeric expression columns found.")
  }
  
  sample_names <- colnames(expr_data) 
  splitted <- strsplit(sample_names, "\\.")
  diag_vec <- sapply(splitted, function(x) x[1])
  treat_vec <- sapply(splitted, function(x) if (length(x) >= 4) x[4] else NA)
  
  keep <- rep(TRUE, length(sample_names))
  if (!is.null(filter_diag) && length(filter_diag) > 0) {
    keep <- keep & (diag_vec %in% filter_diag)
  }
  if (!is.null(filter_treat) && length(filter_treat) > 0) {
    keep <- keep & (treat_vec %in% filter_treat)
  }
  if (!any(keep)) {
    stop("No samples left after filtering!")
  }

  expr_data <- expr_data[, keep, drop = FALSE]
  diag_vec  <- diag_vec[keep]
  treat_vec <- treat_vec[keep]
  sample_names <- sample_names[keep]
  splitted <- splitted[keep]
  
  if (diag) {
    group_labels <- ifelse(diag_vec %in% c("AD", "ADNPS"), "Dis", "NPI")
  } else if (treat) {
    group_labels <- ifelse(treat_vec == "eo", "Treatment", "NoTreatment")
  } else {
    stop("You must set either diag=TRUE or treat=TRUE for this t-test.")
  }
  if (length(unique(group_labels)) < 2) {
    stop("Only one group found after filtering; cannot perform t-test.")
  }

  results <- data.frame(
    Gene = my.data$Gene,
    Accession = my.data$Accession,
    Modifications = my.data$Modifications,
    p_value = NA_real_,
    log2FC  = NA_real_
  )

  if (paired && diag) {
    warning("Paired t-test typically makes sense for treat grouping, but diag=TRUE was set. Will attempt it anyway.")
  }
  
  if (paired) {
    individual_id <- sapply(splitted, function(x) {
      if (length(x) < 3) {
        paste(x, collapse=".")
      } else {
        paste(x[1:3], collapse=".")
      }
    })
    g1_name <- unique(group_labels)[1]
    g2_name <- unique(group_labels)[2]
    
    for (i in seq_len(nrow(expr_data))) {
      row_values <- as.numeric(expr_data[i, ])
      pairVals_g1 <- c()
      pairVals_g2 <- c()

      all_individuals <- unique(individual_id)
      for (indiv in all_individuals) {
        idx_g1 <- which(individual_id == indiv & group_labels == g1_name)
        idx_g2 <- which(individual_id == indiv & group_labels == g2_name)
        if (length(idx_g1) == 1 && length(idx_g2) == 1) {
          val1 <- row_values[idx_g1]
          val2 <- row_values[idx_g2]
          if (!is.na(val1) && !is.na(val2)) {
            pairVals_g1 <- c(pairVals_g1, val1)
            pairVals_g2 <- c(pairVals_g2, val2)
          }
        }
      }
      if (length(pairVals_g1) > 1 && length(pairVals_g2) > 1) {
        t_test_res <- tryCatch(
          t.test(pairVals_g1, pairVals_g2, paired = TRUE, var.equal = FALSE),
          error = function(e) {
            warning(paste("Paired t.test failed in row", i, ":", e$message))
            return(NULL)
          }
        )
        if (!is.null(t_test_res)) {
          p_val <- t_test_res$p.value
          log2fc_val <- mean(pairVals_g1) - mean(pairVals_g2)
          if (reverse_fc) {
            log2fc_val <- -log2fc_val
          }
          results$p_value[i] <- p_val
          results$log2FC[i]  <- log2fc_val
        }
      } else {
        results$p_value[i] <- NA
        results$log2FC[i]  <- NA
      }
    }
    
  } else {
    uni_group_labels <- unique(group_labels)
    uni_group_labels <- uni_group_labels[order(uni_group_labels)]
    g1_name <- uni_group_labels[1]
    g2_name <- uni_group_labels[2]

    for (i in seq_len(nrow(expr_data))) {
      row_values <- as.numeric(expr_data[i, ])
      group1_values <- row_values[group_labels == g1_name]
      group2_values <- row_values[group_labels == g2_name]
      if (sum(!is.na(group1_values)) > 1 && sum(!is.na(group2_values)) > 1) {
        t_test_res <- tryCatch(
          t.test(group1_values, group2_values, paired = FALSE, var.equal = FALSE),
          error = function(e) {
            warning(paste("Unpaired t.test failed in row", i, ":", e$message))
            return(NULL)
          }
        )
        if (!is.null(t_test_res)) {
          p_val <- t_test_res$p.value
          log2fc_val <- mean(group1_values, na.rm=TRUE) - mean(group2_values, na.rm=TRUE)
          if (reverse_fc) {
            log2fc_val <- -log2fc_val
          }
          results$p_value[i] <- p_val
          results$log2FC[i]  <- log2fc_val
        }
      } else {
        results$p_value[i] <- NA
        results$log2FC[i]  <- NA
      }
    }
  }
  
  results$padj <- p.adjust(results$p_value, method = "BH")
  final <- cbind(
    results[, c("Gene", "Accession", "Modifications",
                "p_value", "padj", "log2FC")],
    expr_data
  )
  
  return(final)
}

generate_DE_list <- function(norm.data){
  DEG <- list()
  DEG[["ADvsNCI"]] <- perform_t_test(norm.data, diag = T)
  DEG[["ADvsNCI_|_onlyEO"]] <- perform_t_test(norm.data, diag = T, filter_treat = "eo")
  DEG[["ADvsNCI_|_onlyCL"]] <- perform_t_test(norm.data, diag = T, filter_treat = "cl")
  DEG[["EOvsCL"]] <- perform_t_test(norm.data, treat = T, paired = T, reverse_fc = T)
  DEG[["EOvsCL_|_onlyAD"]] <- perform_t_test(norm.data, treat = T, paired = T, reverse_fc = T, filter_diag = c("AD","ADNPS"))
  DEG[["EOvsCL_|_onlyNCI"]] <- perform_t_test(norm.data, treat = T, paired = T, reverse_fc = T, filter_diag = "NCI")
  
  return(DEG)
}
EV.DEG <- generate_DE_list(EV[["norm"]])
Org.DEG <- generate_DE_list(Org[["norm"]])

EV.DEGp005 <- EV.DEG
for (i in 1:6) {
  EV.DEGp005[[i]] <- EV.DEGp005[[i]][which(EV.DEGp005[[i]]$padj < 0.05 
                                           & abs(EV.DEGp005[[i]]$log2FC)>1),]
}
Org.DEGp005 <- Org.DEG
for (i in 1:6) {
  Org.DEGp005[[i]] <- Org.DEGp005[[i]][which(Org.DEGp005[[i]]$padj < 0.05 
                                             & abs(Org.DEGp005[[i]]$log2FC)>1),]
}

# ----------------------------
# SECTION 11 – AUXILIARY ANALYSES (MARKER/TREATMENT OVERLAP)
# ----------------------------
check_treat_vs_marker <- function(DEG_list){
  temp_AD_marker <- DEG_list[["ADvsNCI_|_onlyCL"]] %>%
    filter(padj < 0.05, abs(log2FC) > log2(1.5)) %>%
    select(1:6)
  
  temp_treat_effct <- DEG_list[["EOvsCL_|_onlyAD"]] %>%
    filter(padj < 0.05, abs(log2FC) > log2(1.5)) %>%
    select(1:6)
  
  colnames(temp_AD_marker)[4:6] <- c("AD.p_value", "AD.padj", "AD.log2FC") 
  colnames(temp_treat_effct)[4:6] <- c("EO.p_value", "EO.padj", "EO.log2FC") 
  
  treat_vs_marker <- merge(temp_AD_marker, temp_treat_effct[,c(1,4:6)], by='Gene',all=F)
  treat_vs_marker$trend <- treat_vs_marker$AD.log2FC * treat_vs_marker$EO.log2FC
  
  treat_vs_marker <- treat_vs_marker[which(treat_vs_marker$trend<0),]
  return(treat_vs_marker)
}
EV.DEG[["treat_vs_marker"]] <- check_treat_vs_marker(EV.DEG)
Org.DEG[["treat_vs_marker"]] <- check_treat_vs_marker(Org.DEG)

# ----------------------------
# SECTION 12 – RESPONDER VS UNRESPONDER
# ----------------------------
unresponder <- c("AD.17","AD.18","ADNPS.19","ADNPS.20")

perform_t_test_unresponder <- function(my.data = EV[["norm"]], keep_treat ="cl") {
  
  expr_data <- my.data[, sapply(my.data, is.numeric), drop = FALSE]
  
  sample_names <- colnames(expr_data) 
  splitted <- strsplit(sample_names, "\\.")
  diag_vec <- sapply(splitted, function(x) x[1])
  treat_vec <- sapply(splitted, function(x) x[4])
  sample_id <- sapply(splitted, function(x) x[3])
  sample_name <- paste0(diag_vec,".",sample_id)
  
  keep <- rep(TRUE, length(sample_names))
  keep <- keep & (treat_vec %in% keep_treat)
  expr_data <- expr_data[, keep, drop = FALSE]
  diag_vec  <- diag_vec[keep]
  treat_vec <- treat_vec[keep]
  sample_names <- sample_names[keep]
  splitted <- splitted[keep]
  
  group_labels <- ifelse(sample_name %in% unresponder, "Unres", "Res")
  
  
  results <- data.frame(
    Gene = my.data$Gene,
    Accession = my.data$Accession,
    Modifications = my.data$Modifications,
    p_value = NA_real_,
    log2FC  = NA_real_
  )
  
  uni_group_labels <- unique(group_labels)
  uni_group_labels <- uni_group_labels[order(uni_group_labels)]
  g1_name <- uni_group_labels[1]
  g2_name <- uni_group_labels[2]
  
  for (i in seq_len(nrow(expr_data))) {
    row_values <- as.numeric(expr_data[i, ])
    group1_values <- row_values[group_labels == g1_name]
    group2_values <- row_values[group_labels == g2_name]
    if (sum(!is.na(group1_values)) > 1 && sum(!is.na(group2_values)) > 1) {
      t_test_res <- tryCatch(
        t.test(group1_values, group2_values, paired = FALSE, var.equal = FALSE),
        error = function(e) {
          warning(paste("Unpaired t.test failed in row", i, ":", e$message))
          return(NULL)
        }
      )
      if (!is.null(t_test_res)) {
        p_val <- t_test_res$p.value
        log2fc_val <- mean(group2_values, na.rm=TRUE) - mean(group1_values, na.rm=TRUE)
        results$p_value[i] <- p_val
        results$log2FC[i]  <- log2fc_val
      }
    } else {
      results$p_value[i] <- NA
      results$log2FC[i]  <- NA
    }
  }
  
  results$padj <- p.adjust(results$p_value, method = "BH")
  final <- cbind(
    results[, c("Gene", "Accession", "Modifications",
                "p_value", "padj", "log2FC")],
    expr_data
  )
  
  return(final)
}

unrespond <- list()
unrespond[["EV_DEG"]] <- perform_t_test_unresponder(EV[["norm"]])
unrespond[["EV_DEGp005"]] <- unrespond[["EV_DEG"]][which(unrespond[["EV_DEG"]]$padj<0.05 & abs(unrespond[["EV_DEG"]]$log2FC)>1),]
unrespond[["Org_DEG"]] <- perform_t_test_unresponder(Org[["norm"]])
unrespond[["Org_DEGp005"]] <- unrespond[["Org_DEG"]][which(unrespond[["Org_DEG"]]$padj<0.05 & abs(unrespond[["Org_DEG"]]$log2FC)>1),]

# ----------------------------
# SECTION 13 – VISUALISATIONS
# ----------------------------
# 13A) Batch & QC boxplots -------------------------------------------
batch_check_boxplot <- function(mydata = EV[["norm"]],
                                extract_pattern = "(?<=\\.EV\\.)\\d+"){
  long_data <- reshape2::melt(mydata, id.vars = c("Gene", "Accession", "Modifications"), 
                              variable.name = "Sample", value.name = "Expression")
  
  long_data <- long_data %>%
    mutate(
      SampleNum = as.numeric(str_extract(Sample, extract_pattern)),
      Group = case_when(
        str_detect(Sample, "^NCI") ~ "NCI",
        TRUE ~ "AD_ADNPS"
      ),
      
      Ending = ifelse(str_ends(Sample, "\\.cl"), "cl", "eo"),
      Is123 = ifelse(SampleNum %in% c(1, 2, 3), TRUE, FALSE),
      
      SortKey = paste0(Is123, "_", Ending, "_", Group, "_", sprintf("%02d", SampleNum)),
      
      ColorGroup = ifelse(Is123, "#66c2a6", "#8ea0cb")
    ) %>%
    arrange(desc(Is123), Ending, Group, SampleNum) %>%
    mutate(Sample = factor(Sample, levels = unique(Sample)))

  ggplot(long_data, aes(x = Sample, y = Expression, fill = ColorGroup)) +
    geom_boxplot(outlier.size = 0.05, na.rm = TRUE) +
    scale_fill_manual(values = c("#66c2a6" = "#66c2a6", "#8ea0cb" = "#8ea0cb")) +
    theme_minimal() +
    labs(x = "Sample", y = "Expression Level") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}
# Save boxplots
batch_check_boxplot(mydata = EV[["norm"]],
                    extract_pattern = "(?<=\\.EV\\.)\\d+")
ggsave("Figure1E_boxplot_EV_Expression_Distribution.pdf",width=14, height = 5)
batch_check_boxplot(mydata = Org[["norm"]],
                    extract_pattern = "(?<=\\.Org\\.)\\d+")
ggsave("Figure1C_boxplot_Org_Expression_Distribution.pdf",width=14, height = 5)


# 13B) Volcano plots --------------------------------------------------

volcano_plot <- function(deg_data = EV.DEG[["ADvsNCI_|_onlyCL"]], titlename, p_thre = 0.05, log2fc_thre = 1){
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
      panel.grid.major = element_blank(),  # This will remove the major grid lines
      panel.grid.minor = element_blank(),  # This will remove the minor grid lines
      legend.title = element_blank() # This will remove the legend titles if needed
    ) +
    geom_hline(yintercept = -log10(0.05),linetype=2,cex=0.5,color = "grey")+  #添加辅助线
    geom_vline(xintercept = c(-1,1),linetype=2,cex=0.5,color = "grey")
}
# Save volcano plots
volcano_plot(EV.DEG[["ADvsNCI_|_onlyCL"]], titlename = "Untreated EVs - AD vs NCI")
ggsave("Figure3D_volcano_Untreated EVs_AD vs NCI.pdf",width = 6, height = 6)
volcano_plot(Org.DEG[["ADvsNCI_|_onlyCL"]], titlename = "Untreated Organoids - AD vs NCI")
ggsave("Figure3A_volcano_Untreated Orgs_AD vs NCI.pdf",width = 6, height = 6)

volcano_plot(EV.DEG[["EOvsCL"]], titlename = "All EVs - EO vs CL",p_thre=0.001,log2fc_thre =2.5)
ggsave("Figure5A_volcano_All EVs_EO vs CL.pdf",width = 6, height = 6)
volcano_plot(Org.DEG[["EOvsCL"]], titlename = "All Organoids - EO vs CL",p_thre=0.01,log2fc_thre =1.5)
ggsave("Figure4A_volcano_All Orgs_EO vs CL.pdf",width = 6, height = 6)

volcano_plot(EV.DEG[["EOvsCL_|_onlyAD"]],titlename = "AD EVs - EO vs CL",p_thre=10e-6,log2fc_thre =2.5)
ggsave("Figure5_volcano_AD EVs_EO vs CL.pdf",width = 6, height = 6)
volcano_plot(Org.DEG[["EOvsCL_|_onlyAD"]],titlename = "AD Organoids - EO vs CL",p_thre=10e-6,log2fc_thre =2.5)
ggsave("Figure4_volcano_AD Orgs_EO vs CL.pdf",width = 6, height = 6)

volcano_plot(EV.DEG[["EOvsCL_|_onlyNCI"]], titlename = "NCI EVs - EO vs CL")
ggsave("Figure5_volcano_NCI EVs_EO vs CL.pdf",width = 6, height = 6)
volcano_plot(Org.DEG[["EOvsCL_|_onlyNCI"]], titlename = "NCI Organoids - EO vs CL")
ggsave("Figure4_volcano_NCI Orgs_EO vs CL.pdf",width = 6, height = 6)


# 13C) Violin plots ---------------------------------------------------
Violin_Protein_Expression_Facet <- function(genes,
                                            data,
                                            DEGdata,
                                            mytitle,
                                            log2 = TRUE) {
  
  # Get the data
  df_sub <- data[data$Gene %in% genes, ]
  long_data <- reshape2::melt(
    df_sub,
    id.vars      = c("Gene", "Accession", "Modifications"),
    variable.name = "Sample",
    value.name    = "Expression"
  )
  
  splitted  <- strsplit(as.character(long_data$Sample), "\\.")
  diag_vec  <- sapply(splitted, `[`, 1)
  diag_vec  <- ifelse(diag_vec == "NCI", "NCI", "AD")
  treat_vec <- sapply(splitted, `[`, 4)
  long_data$Group <- paste0(diag_vec, ".", treat_vec)
  if(log2){long_data$Expression <- log2(long_data$Expression)}
  long_data$Group <- factor(long_data$Group, levels = sort(unique(long_data$Group)))
  long_data$SubjectID <- sapply(splitted, function(x) {
    if (length(x) < 3) paste(x, collapse = ".")
    else                paste(x[1:3], collapse = ".")
  })
  
  # Get the p value
  comparisons_tbl <- tibble::tribble(
    ~pair_id, ~group1, ~group2, ~tbl_name,
    "1vs2",   "AD.cl",      "AD.eo",    "EOvsCL_|_onlyAD",
    "3vs4",   "NCI.cl",     "NCI.eo",   "EOvsCL_|_onlyNCI",
    "1vs3",   "AD.cl",      "NCI.cl",   "ADvsNCI_|_onlyCL",
    "2vs4",   "AD.eo",      "NCI.eo",   "ADvsNCI_|_onlyEO"
  )
  library(dplyr)
  library(purrr)
  
  pval_df <- purrr::map_dfr(seq_len(nrow(comparisons_tbl)), function(i) {
    comp    <- comparisons_tbl[i, ]
    deg_tbl <- DEGdata[[ comp$tbl_name ]]
    
    deg_tbl %>%
      filter(Gene %in% long_data$Gene) %>%
      transmute(
        Gene   = Gene,
        pair   = comp$pair,
        group1 = comp$group1,
        group2 = comp$group2,
        p      = coalesce(padj, 1)    # NA -> 1
      )
  })
  
  # calculate y.position
  expr_range <- long_data %>%
    group_by(Gene) %>%
    summarise(
      min_expr = min(Expression, na.rm = TRUE),
      max_expr = max(Expression, na.rm = TRUE),
      .groups  = "drop"
    )
  
  # combine p and y.position
  pval_df <- pval_df %>%
    left_join(expr_range, by = "Gene") %>%
    group_by(Gene) %>%
    mutate(
      offset_level = case_when(
        row_number() <= 2 ~ 1,
        row_number() == 3 ~ 2,
        TRUE              ~ 3
      ),
      y.position = max_expr + (max_expr - min_expr) / 6 * offset_level,
      p.signif = cut(
        p,
        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
        labels = c("***", "**", "*", "ns")
      )
    ) %>%
    ungroup() %>%
    select(-min_expr, -max_expr, -offset_level)
  
  df_lines <- long_data %>%
    filter(Group %in% c("AD.cl","AD.eo","NCI.cl","NCI.eo"))
  
  p <- ggplot(long_data,
              aes(x = Group, y = Expression, fill = Group)) +
    gghalves::geom_half_violin(position = position_nudge(x = -0.1)) +
    geom_boxplot(outlier.size = 0.1, width = 0.15) +
    geom_point(data     = df_lines,
               aes(x = Group, y = Expression),
               position = position_nudge(x = 0.15),
               shape    = 21,
               alpha    = 0.7,
               inherit.aes = FALSE) +
    
    geom_line(data     = df_lines,
              aes(x = Group, y = Expression, group = SubjectID),
              position    = position_nudge(x = 0.15),
              color       = "gray70",
              alpha       = 0.6,
              inherit.aes = FALSE) +
    ggpubr::stat_pvalue_manual( data         = pval_df,
                                label        = "p.signif",
                                xmin         = "group1",
                                xmax         = "group2",
                                y.position   = "y.position",
                                bracket.size = 0.3,
                                tip.length   = 0.01,
                                inherit.aes  = FALSE ) +
    facet_wrap(~ Gene, nrow = 2, ncol = ceiling(length(genes)/2), scales="free_y") +
    theme_bw() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c(
      "NCI.cl"    = "skyblue",
      "NCI.eo"    = "#E0F7FF",
      "AD.cl"     = "#FF9999",
      "AD.eo"     = "#FFF0F0"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    labs(title = mytitle, x=NULL, y="Protein Expression Level (log2)")
  
  return(p)
}

EV_gene <- c("NCAM1","ATP1A3","CD63","CD81","CD9","L1CAM","PDCD6IP","TSG101")
Violin_Protein_Expression_Facet(EV_gene, EV[["filtermiss"]], EV.DEG, "EVs Protein Expression Across Markers")
ggsave("Figure2C_violin_EVs_markers_before_norm.pdf",width = 12, height = 7)
Violin_Protein_Expression_Facet(EV_gene, EV[["norm"]], EV.DEG, "EVs Protein Expression Across Markers", log2 = F)
ggsave("Figure2C_violin_EVs_markers_after_norm.pdf",width = 12, height = 7)

Org_gene <- c("MAP2","TUBB3","MAPT","SYN1","GAP43","TPH1","NES","MAOA","HTR7","MAOB")
Violin_Protein_Expression_Facet(Org_gene, Org[["filtermiss"]], Org.DEG, "Organoids Protein Expression Across Markers")
ggsave("Figure2A_violin_Orgs_markers_before_norm.pdf",width = 12, height = 7)
Violin_Protein_Expression_Facet(Org_gene, Org[["norm"]], Org.DEG, "Organoids Protein Expression Across Markers", log2 = F)
ggsave("Figure2A_violin_Orgs_markers_after_norm.pdf",width = 12, height = 7)

EV_relief <- EV.DEG[["treat_vs_marker"]]$Gene
Violin_Protein_Expression_Facet(EV_relief[1:12], EV[["filtermiss"]], EV.DEG, "EVs Protein Expression Across DEPs")
#ggsave("Interest_figure_EVs_relief_before_norm1.pdf",width = 12, height = 7)
Violin_Protein_Expression_Facet(EV_relief[13:24], EV[["filtermiss"]], EV.DEG, "EVs Protein Expression Across DEPs")
#ggsave("Interest_figure_EVs_relief_before_norm2.pdf",width = 12, height = 7)

Org_relief <- Org.DEG[["treat_vs_marker"]]$Gene
Violin_Protein_Expression_Facet(Org_relief[1:13], Org[["filtermiss"]], Org.DEG, "Orgs Protein Expression Across DEPs")
#ggsave("Interest_figure_Orgs_reliefs_before_norm.pdf",width = 14, height = 7)

# 12D) PCA plots based on the marker
EO_pca_plot <- function(my.data, my.title = "PCA plot"
) {
  numeric_data <- my.data[, sapply(my.data, is.numeric), drop = FALSE]
  numeric_data <- na.omit(numeric_data)
  sample_names <- colnames(numeric_data)
  splitted <- strsplit(sample_names, "\\.")
  eo_vec <- sapply(splitted, function(x){x[4]})
  pca_res <- prcomp(t(numeric_data), center = TRUE, scale. = TRUE)
  pca_df  <- as.data.frame(pca_res$x)  
  group <- eo_vec
  group_factor <- factor(group)
  
  pc1 <- pca_df[,"PC1"]
  pc2 <- pca_df[,"PC2"]
  var_explained <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
  df_pca <- data.frame(
    PC1 = pc1,
    PC2 = pc2,
    Group = group_factor
  )
  pc1_label <- paste0("PC1: ", sprintf("%.2f", var_explained[1]), "%")
  pc2_label <- paste0("PC2: ", sprintf("%.2f", var_explained[2]), "%")
  
  ggplot(df_pca, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 2) +
    labs(
      x = pc1_label,
      y = pc2_label,
      title = my.title
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) + 
    scale_color_manual(values = c("cl" = "#FCCDE5", "eo"="#B3DE69"))
  
  
}
temp <- Org[["norm"]][which(Org[["norm"]]$Gene %in% Org_gene),]
EO_pca_plot(temp, "Org - Based on Ten Markers")
ggsave("Figure2B_PCA_Org_with 10 Markers.pdf",width = 6, height = 6)
temp <- EV[["norm"]][which(EV[["norm"]]$Gene %in% EV_gene),]
EO_pca_plot(temp, "EV - Based on Eight Markers")
ggsave("Figure2D_PCA_EV_with 8 Markers.pdf",width = 6, height = 6)



# 13D) Heatmaps -------------------------------------------------------

heatmap_plot <- function(my.data, show_protein_label = T, mytitle = "") {
  
  numeric_data <- my.data[, sapply(my.data, is.numeric), drop = FALSE]
  rownames(numeric_data) <- my.data$Gene
  
  sample_names <- colnames(numeric_data)
  splitted     <- strsplit(sample_names, "\\.")
  diag_vec <- sapply(splitted, function(x) x[1])
  treat_vec <- sapply(splitted, function(x) {
    if (length(x) >= 4) x[4] else NA
  })
  
  diag_vec2 <- ifelse(diag_vec=="NCI","NCI","AD")
  sample_id <- sapply(splitted, function(x) x[3])
  sample_label <- paste0(diag_vec2,"-",sample_id)
  
  label_df <- data.frame(SampleName = sample_label)
  label_with_cli <- label_df %>%
    left_join(sample.clinical %>% select(SampleName, Age, Sex), by = "SampleName")
  sample_age <- label_with_cli$Age
  sample_Sex <- label_with_cli$Sex
  
  batch_vec <- sapply(splitted, function(x) {
    if (length(x) < 3) {
      "batch2"
    } else {
      num_part <- suppressWarnings(as.integer(x[3]))
      if (!is.na(num_part) && num_part %in% c(1, 2, 3)) "batch1" else "batch2"
    }
  })
  
  annotation_df <- data.frame(row.names = sample_names)
  annotation_df[["Diag"]] <- factor(diag_vec)
  annotation_df[["Treat"]] <- factor(treat_vec)
  annotation_df[["Batch"]] <- factor(batch_vec)
  annotation_df[["Age"]] <- sample_age
  annotation_df[["Sex"]] <- factor(sample_Sex)
  
  annotation_colors <- list(Diag = c("AD"= "#FF9999","ADNPS" = "#8B0000","NCI"= "skyblue"),
                            Treat = c("cl" = "#FCCDE5", "eo"="#B3DE69"),
                            Batch = c("batch1" = "#66c2a6", "batch2" = "#8ea0cb"),
                            Sex = c("M"="#FB5012", "F"="#9590FF"))
  
  pheatmap(
    numeric_data,
    scale            = "row",
    show_rownames    = show_protein_label, 
    show_colnames    = TRUE, 
    cluster_rows     = TRUE,
    cluster_cols     = T,
    annotation_col   = if (ncol(annotation_df) > 0) annotation_df else NULL,
    annotation_colors = annotation_colors,
    main   = mytitle,
    na_col = "grey"
  )
}
temp_data <- Org.DEGp005[["ADvsNCI_|_onlyCL"]]
dev.off()
pdf("Figure3B_heatmap_Untreated Organoids_AD vs NCI.pdf",width = 9, height = 7)
heatmap_plot(temp_data[,c(1:3,7:ncol(temp_data))], mytitle="Untreated Organoids - AD vs NCI")
dev.off()

temp_data <- EV.DEGp005[["ADvsNCI_|_onlyCL"]]
pdf("Figure3E_heatmap_Untreated EVs_AD vs NCI.pdf",width = 9, height = 7)
heatmap_plot(temp_data[,c(1:3,7:ncol(temp_data))], mytitle="Untreated EVs - AD vs NCI")
dev.off()

temp_data <- Org.DEGp005[["EOvsCL"]]
pdf("Figure4B_heatmap_All Organoids_EO vs CL.pdf",width = 14, height = 7)
heatmap_plot(temp_data[,c(1:3,7:ncol(temp_data))], show_protein_label = F, mytitle="All Organoids - EO vs CL")
dev.off()

temp_data <- EV.DEGp005[["EOvsCL"]]
pdf("Figure5B_heatmap_All EVs_EO vs CL.pdf",width = 14, height = 7)
heatmap_plot(temp_data[,c(1:3,7:ncol(temp_data))], show_protein_label = F, mytitle="All EVs - EO vs CL")
dev.off()


# 13E) Annotated PCA for responder analysis ---------------------------

labeled_diag_treat_pca_plot <- function(my.data, label_sample = "",my.title = "PCA plot"
) {
  numeric_data <- my.data[, sapply(my.data, is.numeric), drop = FALSE]
  numeric_data <- na.omit(numeric_data)
  
  sample_names <- colnames(numeric_data)
  splitted <- strsplit(sample_names, "\\.")
  diag_vec <- sapply(splitted, function(x) x[1])
  treat_vec <- sapply(splitted, function(x) if (length(x) >= 4) x[4] else NA)
  batch_vec <- sapply(splitted, function(x) {
    if (length(x) < 3) {
      return("batch2")
    }
    num_part <- suppressWarnings(as.integer(x[3]))
    if (!is.na(num_part) && num_part %in% c(1, 2, 3)) {
      "batch1"
    } else {
      "batch2"
    }
  })
  
  pca_res <- prcomp(t(numeric_data), center = TRUE, scale. = TRUE)
  pca_df  <- as.data.frame(pca_res$x) 
  group <- paste(diag_vec, treat_vec, sep = ".")
  group_factor <- factor(group)
  
  pc1 <- pca_df[,"PC1"]
  pc2 <- pca_df[,"PC2"]
  var_explained <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
  
  df_pca <- data.frame(
    PC1 = pc1,
    PC2 = pc2,
    Group = group_factor,
    SampleName = sample_names
  )
  
  df_pca$SampleName <- ifelse(df_pca$SampleName %in% label_sample, df_pca$SampleName, NA)
  
  df_pca$SampleName <- sapply(df_pca$SampleName, function(x) {
    if (is.na(x)) return(NA) 
    
    parts <- unlist(strsplit(x, "\\."))
    
    if (length(parts) >= 3) {
      group <- parts[1]
      id <- parts[3]
      subid <- if (length(parts) >= 5) parts[5] else NULL
      
      if (!is.null(subid)) {
        return(paste0(group, id, ".", subid))
      } else {
        return(paste0(group, id))
      }
    } else {
      return(x) 
    }
  })
  
  pc1_label <- paste0("PC1: ", sprintf("%.2f", var_explained[1]), "%")
  pc2_label <- paste0("PC2: ", sprintf("%.2f", var_explained[2]), "%")
  
  ggplot(df_pca, aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = Group), shape = 21, size = 3, stroke = 0.6, color = "gray30") + 
    geom_text_repel(aes(label = SampleName), size = 3, max.overlaps = 100, color = "black") + 
    labs(
      x = pc1_label,
      y = pc2_label,
      title = my.title
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = c(
      "NCI.cl"    = "skyblue",
      "NCI.eo"    = "#E0F7FF",
      "AD.cl"     = "#FF9999",
      "AD.eo"     = "#FFF0F0",
      "ADNPS.cl"  = "#8B0000",
      "ADNPS.eo"  = "#F5CCCC"
    ))
  
}
mylabel_sample = c("ADNPS.Org.20.eo", "ADNPS.Org.19.eo", "AD.Org.18.eo", "AD.Org.17.eo", 
                   "NCI.Org.3.eo.2", "NCI.Org.3.eo.1", "NCI.Org.3.eo.3","NCI.Org.2.eo.1")
temp_data <- Org.DEGp005[["EOvsCL"]]
labeled_diag_treat_pca_plot(temp_data[,c(1:3,7:ncol(temp_data))], label_sample = mylabel_sample,
                            my.title = "All Organoids − EO vs CL")
ggsave("Figure4C_PCA_Org_EO_unresponsor.pdf",width = 8, height = 6)

mylabel_sample = c("ADNPS.EV.21.eo", "ADNPS.EV.20.eo", "ADNPS.EV.19.eo", "AD.EV.18.eo", "AD.EV.17.eo")
temp_data <- EV.DEGp005[["EOvsCL"]]
labeled_diag_treat_pca_plot(temp_data[,c(1:3,7:ncol(temp_data))], label_sample = mylabel_sample,
                            my.title = "All EVs − EO vs CL")
ggsave("Figure5C_PCA_EV_EO_unresponsor.pdf",width = 8, height = 6)


# 13F) Violin Plot for responder analysis ---------------------------
Violin_Protein_Expression_Facet_unresponder <- function(genes,
                                                        data,
                                                        DEGdata,
                                                        mytitle,
                                                        log2 = TRUE) {
  
  # Get the data
  df_sub <- data[data$Gene %in% genes, ]
  long_data <- reshape2::melt(
    df_sub,
    id.vars      = c("Gene", "Accession", "Modifications"),
    variable.name = "Sample",
    value.name    = "Expression"
  )
  
  splitted  <- strsplit(as.character(long_data$Sample), "\\.")
  diag_vec  <- sapply(splitted, `[`, 1)
  diag_vec  <- ifelse(diag_vec == "NCI", "NCI", "AD")
  treat_vec <- sapply(splitted, `[`, 4)
  long_data$Group <- paste0(diag_vec, ".", treat_vec)
  if(log2){long_data$Expression <- log2(long_data$Expression)}
  long_data$Group <- factor(long_data$Group, levels = sort(unique(long_data$Group)))
  long_data$SubjectID <- sapply(splitted, function(x) {paste0(x[1],".",x[3])})
  
  # Get the p value
  comparisons_tbl <- tibble::tribble(
    ~pair_id, ~group1, ~group2, ~tbl_name,
    "1vs2",   "AD.cl",      "AD.eo",    "EOvsCL_|_onlyAD",
    "3vs4",   "NCI.cl",     "NCI.eo",   "EOvsCL_|_onlyNCI",
    "1vs3",   "AD.cl",      "NCI.cl",   "ADvsNCI_|_onlyCL",
    "2vs4",   "AD.eo",      "NCI.eo",   "ADvsNCI_|_onlyEO"
  )
  
  pval_df <- purrr::map_dfr(seq_len(nrow(comparisons_tbl)), function(i) {
    comp    <- comparisons_tbl[i, ]
    deg_tbl <- DEGdata[[ comp$tbl_name ]]
    
    deg_tbl %>%
      filter(Gene %in% long_data$Gene) %>%
      transmute(
        Gene   = Gene,
        group1 = comp$group1,
        group2 = comp$group2,
        p      = coalesce(padj, 1)    # NA -> 1
      )
  })
  
  # calculate y.position
  expr_range <- long_data %>%
    group_by(Gene) %>%
    summarise(
      min_expr = min(Expression, na.rm = TRUE),
      max_expr = max(Expression, na.rm = TRUE),
      .groups  = "drop"
    )
  
  # combine p and y.position
  pval_df <- pval_df %>%
    left_join(expr_range, by = "Gene") %>%
    group_by(Gene) %>%
    mutate(
      offset_level = case_when(
        row_number() <= 2 ~ 1,
        row_number() == 3 ~ 2,
        TRUE              ~ 3
      ),
      y.position = max_expr + (max_expr - min_expr) / 6 * offset_level,
      p.signif = cut(
        p,
        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
        labels = c("***", "**", "*", "ns")
      )
    ) %>%
    ungroup() %>%
    select(-min_expr, -max_expr, -offset_level)
  
  df_lines <- long_data %>%
    filter(Group %in% c("AD.cl","AD.eo","NCI.cl","NCI.eo"))
  
  highlight_df <- df_lines %>%
    filter(SubjectID %in% unresponder)
  
  p <- ggplot(long_data,
              aes(x = Group, y = Expression, fill = Group)) +
    gghalves::geom_half_violin(position = position_nudge(x = -0.1)) +
    geom_boxplot(outlier.size = 0.1, width = 0.15) +
    geom_point(data     = df_lines,
               aes(x = Group, y = Expression),
               position = position_nudge(x = 0.15),
               shape    = 21,
               alpha    = 0.7,
               inherit.aes = FALSE) +
    
    geom_line(data     = df_lines,
              aes(x = Group, y = Expression, group = SubjectID),
              position    = position_nudge(x = 0.15),
              color       = "gray70",
              alpha       = 0.6,
              inherit.aes = FALSE) +
    ggpubr::stat_pvalue_manual( data         = pval_df,
                                label        = "p.signif",
                                xmin         = "group1",
                                xmax         = "group2",
                                y.position   = "y.position",
                                bracket.size = 0.3,
                                tip.length   = 0.01,
                                inherit.aes  = FALSE ) +
    geom_point(data     = highlight_df,
               aes(x = Group, y = Expression),
               position    = position_nudge(x = 0.15),
               shape       = 21,
               fill        = "darkred",
               color       = "darkred",
               size        = 2,
               inherit.aes = FALSE) +
    geom_line(
      data     = highlight_df,
      aes(x = Group, y = Expression, group = SubjectID),
      position = position_nudge(x = 0.15),
      color    = "darkred",
      size     = 0.7,
      alpha    = 0.8,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ Gene, nrow = 2, ncol = ceiling(length(genes)/2), scales="free_y") +
    theme_bw() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c(
      "NCI.cl"    = "skyblue",
      "NCI.eo"    = "#E0F7FF",
      "AD.cl"     = "#FF9999",
      "AD.eo"     = "#FFF0F0"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    labs(title = mytitle, x=NULL, y="Protein Expression Level (log2)")
  
  return(p)
}

EV_relief <- EV.DEG[["treat_vs_marker"]]$Gene
Violin_Protein_Expression_Facet_unresponder(EV_relief[1:12], EV[["filtermiss"]], EV.DEG, "EVs Protein Expression Across DEPs")
ggsave("FigureX_violin_EVs_Unresponder_relief_before_norm1.pdf",width = 12, height = 7)
Violin_Protein_Expression_Facet_unresponder(EV_relief[13:24], EV[["filtermiss"]], EV.DEG, "EVs Protein Expression Across DEPs")
ggsave("FigureX_violin_EVs_Unresponder_relief_before_norm2.pdf",width = 12, height = 7)

Org_relief <- Org.DEG[["treat_vs_marker"]]$Gene
Violin_Protein_Expression_Facet_unresponder(Org_relief[1:13], Org[["filtermiss"]], Org.DEG, "Orgs Protein Expression Across DEPs")
ggsave("FigureX_violin_Orgs_Unresponder_relief_before_norm.pdf",width = 14, height = 7)


# 13G) Volcano Plot for responder analysis ---------------------------
volcano_plot_unresponder <- function(deg_data = unrespond[["EV_DEG"]], titlename, p_thre = 0.05, log2fc_thre = 1){
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

volcano_plot_unresponder(unrespond[["EV_DEG"]], titlename = "EO Treated AD EV - Unresponder vs Responder")
ggsave("FigureX_Volcano_EV_Unresponder_DEanalysis.pdf", width=6, height = 6)
volcano_plot_unresponder(unrespond[["Org_DEG"]], titlename = "EO Treatd AD Organoid - Unresponder vs Responder")
ggsave("FigureX_Volcano_Org_Unresponder_DEanalysis.pdf", width=6, height = 6)


# 13H) Heatmap for responder analysis ---------------------------
dev.off()
pdf("FigureX_heatmap_EV_Unresponder vs Responder.pdf",width = 9, height = 7)
temp_data <- unrespond[["EV_DEGp005"]]
heatmap_plot(temp_data[,c(1:3,7:ncol(temp_data))], mytitle="Untreated EVs - Unresponder vs Responder", show_protein_label = F)
dev.off()

pdf("FigureX_heatmap_Org_Unresponder vs Responder.pdf",width = 9, height = 7)
temp_data <- unrespond[["Org_DEGp005"]]
heatmap_plot(temp_data[,c(1:3,7:ncol(temp_data))], mytitle="Untreated Organoids - Unresponder vs Responder", show_protein_label = F)
dev.off()
