###############################################################################
### Generate TPM matrix (Batch corrected) 
### Input : Combined_ARG_VFDB_Harmonized.tsv
###############################################################################
set.seed(44)
library(tidyverse)

# ---- Step 1 | Load harmonized data ----
df <- read.delim("Combined_ARG_VFDB_Harmonized.tsv", sep = "\t", header = TRUE)

# ---- Step 2 | Aggregate TPM per sample × mechanism ----
mech_tpm <- df %>%
  group_by(Sample_ID, Mechanisms) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

# ---- Step 3 | Pivot to wide matrix and keep Sample_ID ----
mech_tpm_wide <- mech_tpm %>%
  pivot_wider(
    id_cols = Sample_ID,
    names_from = Mechanisms,
    values_from = TPM,
    values_fill = 0
  )

# ---- Step 4 | Make Sample_ID the row names and keep it visible for checking ----
rownames(mech_tpm_wide) <- mech_tpm_wide$Sample_ID

# ---- Step 5 | Save clean matrix to TSV ----
write.table(mech_tpm_wide,
            file = "Mechanism_TPM_Matrix.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


###############################################################################
##  CLR Transformation of Mechanism TPM Matrix
################################################################################

# --- Load Required Libraries ---
library(compositions)
library(ggplot2)

# --- Step 1: Load TPM Data ---
# Ensure the first column contains Sample_IDs
tpm <- read.delim("Mechanism_TPM_Matrix.tsv", sep = "\t", header = TRUE, check.names = FALSE)

# If Sample_ID is the first column, set it as rownames
if ("Sample_ID" %in% colnames(tpm)) {
  rownames(tpm) <- tpm$Sample_ID
  tpm$Sample_ID <- NULL
}


# --- Step 2: CLR Transformation ---
# Add a small pseudocount to avoid log(0)
tpm_pseudo <- tpm + 1e-6

# Perform CLR transformation across columns (mechanisms)
clr_matrix <- as.data.frame(apply(tpm_pseudo, 1, function(x) clr(x)))
clr_matrix <- t(clr_matrix)


# --- Step 3: Diagnostics ---
# Distribution plot
#png("Mechanism_CLR_Distribution.png", width = 900, height = 450)
#hist(as.numeric(as.matrix(clr_matrix)),
    #breaks = 60, col = "steelblue", border = "black",
    # main = "Distribution of CLR-transformed values",
    # xlab = "CLR value", ylab = "Frequency")
#dev.off()


# --- FINAL FIGURE ---
library(ggplot2)
library(tidyverse)

# Flatten CLR matrix to a single vector
clr_vals <- as.numeric(as.matrix(clr_matrix))


ggplot(data.frame(CLR = clr_vals), aes(x = CLR)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 80,
                 fill = "gray85", 
                 color = "black", 
                 alpha = 0.9) +
  geom_density(color = "firebrick4", size = 1.2, adjust = 1.2) +
  labs(
    x = "CLR value",
    y = "Density"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.line = element_line(size = 0.8, color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )


ggsave("Mechanism_CLR_Distribution_Density.png", width = 6, height = 5, dpi = 600)



#####bOXPLOT
library(tidyverse)

# -------------------------------
# 1. Convert CLR matrix to long format
# -------------------------------
df_long <- as.data.frame(clr_matrix) %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Mechanism",
    values_to = "CLR"
  )

# -------------------------------
# 2. Define shortened mechanism labels
# -------------------------------
mechanism_short_labels <- c(
  "Adherence" = "Adherence",
  "antibiotic efflux" = "Antibiotic efflux",
  "antibiotic efflux;antibiotic target alteration" = "Efflux + TargetAlt",
  "antibiotic efflux;reduced permeability to antibiotic" = "Efflux + RedPerm",
  "antibiotic inactivation" = "Antibiotic inactivation",
  "antibiotic inactivation;antibiotic target alteration" = "Inact + TargetAlt",
  "antibiotic target alteration" = "TargetAlt",
  "antibiotic target alteration;antibiotic target replacement" = "TargetAlt + TargetRepl",
  "antibiotic target protection" = "TargetProtect",
  "antibiotic target replacement" = "TargetRepl",
  "Antimicrobial activity Competitive advantage" = "CompetitiveAdv",
  "Biofilm" = "Biofilm",
  "Effector delivery system" = "Effector delivery system",
  "Exoenzyme" = "Exoenzyme",
  "Exotoxin" = "Exotoxin",
  "Immune modulation" = "Immune modulation",
  "Invasion" = "Invasion",
  "Motility" = "Motility",
  "Nutritional Metabolic factor" = "Nutritional Metabolic factor",
  "Others" = "Others",
  "reduced permeability to antibiotic" = "RedPerm",
  "Regulation" = "Regulation",
  "Stress survival" = "Stress survival"
)

df_long$Mechanism_short <- mechanism_short_labels[df_long$Mechanism]

# -------------------------------
# 3. DEFINE FUNCTIONAL GROUPS (ARG / VF / Hybrid)
# -------------------------------
vf_mechs <- c("Adherence", "CompetitiveAdv", "Biofilm", "Effector delivery system", "Exoenzyme",
              "Exotoxin", "Immune modulation", "Invasion", "Motility", "Nutritional Metabolic factor",
              "Regulation", "Stress survival", "Others")

arg_mechs <- c("Antibiotic efflux", "Antibiotic inactivation", "TargetAlt", "TargetRepl",
               "TargetProtect", "RedPerm")

hybrid_mechs <- c("Efflux + TargetAlt", "Efflux + RedPerm",
                  "Inact + TargetAlt", "TargetAlt + TargetRepl")

df_long$Group <- case_when(
  df_long$Mechanism_short %in% vf_mechs ~ "VF",
  df_long$Mechanism_short %in% arg_mechs ~ "ARG",
  df_long$Mechanism_short %in% hybrid_mechs ~ "Hybrid",
  TRUE ~ "Others"
)

# -------------------------------
# 4. COLOR PALETTE
# -------------------------------
group_colors <- c(
  "VF" = "#3B82F6",      # blue
  "ARG" = "#DC2626",     # red
  "Hybrid" = "#7C3AED",  # purple
  "Others" = "gray50"
)

# -------------------------------
# 5. ORDER MECHANISMS BY MEDIAN CLR
# -------------------------------
df_long$Mechanism_short <- reorder(df_long$Mechanism_short, df_long$CLR, FUN = median)

# -------------------------------
# 6. PLOT (Reordered + Colored)
# -------------------------------
p <- ggplot(df_long, aes(x = Mechanism_short, y = CLR, fill = Group)) +
  geom_boxplot(color = "black", outlier.alpha = 0.35, size = 0.6) +
  scale_fill_manual(values = group_colors, name = "Mechanism Class") +
  labs(
    title = "",
    x = "",
    y = "CLR value"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    legend.position = "top"
  )

print(p)

# -------------------------------
# 5. Save the figure
# -------------------------------
ggsave("Mechanism_CLR_Boxplot.png",
       p, width = 7, height = 6, dpi = 600)


###############################################
# --- Step 4: Save Correctly (Preserve Sample IDs + Mechanism Names) ---
write.table(
  clr_matrix,
  file = "Mechanism_CLR.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA,       # ensures rownames (Sample_ID) are saved as first column
  row.names = TRUE,
  fileEncoding = "UTF-8"
)


# --- Step 5: Reload to confirm integrity ---
clr_check <- read.delim("Mechanism_CLR.tsv", sep = "\t", check.names = FALSE, row.names = 1)
cat("Recheck: ", nrow(clr_check), "samples ×", ncol(clr_check), "mechanisms.\n")
stopifnot(identical(rownames(clr_check), rownames(clr_matrix)))


###############################################################################
#PERMANOVA Variance Partitioning for Mechanism CLR Matrix
################################################################################

# --- Load Libraries ---
library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)

# --- Step 1: Load CLR Matrix ---
clr <- read.delim("Mechanism_CLR.tsv", sep = "\t", check.names = FALSE, row.names = 1)

# --- Step 2: Load Metadata ---
meta <- read.csv("Merged_metadata.csv", header = TRUE, check.names = FALSE)


# --- Match Samples ---
common_ids <- intersect(rownames(clr), meta$Sample_ID)
clr <- clr[common_ids, , drop = FALSE]
meta <- meta[meta$Sample_ID %in% common_ids, ]
rownames(meta) <- meta$Sample_ID

# --- Step 3: Prepare distance matrix (Euclidean) ---
clr_dist <- dist(clr, method = "euclidean")

# --- Step 4: Define factors for PERMANOVA ---
factors <- c("Group", "Region", "environmental_material", 
             "Instrument", "Sequencing_Center", "Project_ID", "Country", "Continent")

# Filter to only existing metadata columns
factors <- factors[factors %in% colnames(meta)]

# --- Step 5: Run PERMANOVA for each factor ---
adonis_results <- list()
for (f in factors) {
  if (length(unique(na.omit(meta[[f]]))) > 1) {  # skip if only one level
    res <- adonis2(clr_dist ~ meta[[f]], permutations = 999)
    df <- as.data.frame(res)
    df$Factor <- f
    adonis_results[[f]] <- df
  }
}

# --- Step 6: Combine all results ---
adonis_df <- bind_rows(adonis_results, .id = NULL)
adonis_summary <- adonis_df %>%
  filter(!is.na(F)) %>%
  select(Factor, R2, F, `Pr(>F)`) %>%
  rename(p_value = `Pr(>F)`) %>%
  arrange(desc(R2))

# Save results
write.csv(adonis_summary, "PERMANOVA_Variance_Summary.csv", row.names = FALSE)


# --- Step 7: Plot variance explained ---
tiff("PERMANOVA_Variance_Explained.tiff", width = 5, height = 3, units = "in", res = 600)

ggplot(adonis_summary, aes(x = reorder(Factor, R2), y = R2)) +
  geom_col(fill = "#3182bd", color = "black", width = 0.65) +
  coord_flip() +
  
  # Numerical labels at bar ends
  geom_text(aes(label = sprintf("%.2f", R2)),
            hjust = -0.1,
            size = 2,          # readable size
            fontface = "bold") +
  
  # Clean theme
  theme_classic(base_size = 16) +
  
  labs(
    title = "",
    x = "Factor",
    y = expression(R^2 ~ "(Variance Explained)")
  ) +
  
  theme(
    axis.text = element_text(color = "black", size = 8),
    axis.title = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.grid = element_blank()
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))   # extra room for labels


dev.off()


# --- Step 8: Identify potential batch effects ---
batch_candidates <- adonis_summary %>%
  filter(p_value < 0.05, R2 > 0.01, !Factor %in% c("Group", "environmental_material")) %>%
  arrange(desc(R2))

print(batch_candidates)



####################################################################################################################
# Batch Correction of Mechanism CLR Matrix using Limma (FINAL FIX)
################################################################################

# --- Load Libraries ---
library(limma)
library(vegan)
library(ggplot2)
library(dplyr)

# --- Step 1: Load Data ---
clr <- read.delim("Mechanism_CLR.tsv", sep = "\t", check.names = FALSE, row.names = 1)
meta <- read.csv("Merged_metadata.csv", header = TRUE, check.names = FALSE)

# Match sample order
common_ids <- intersect(rownames(clr), meta$Sample_ID)
clr <- clr[common_ids, ]
meta <- meta[meta$Sample_ID %in% common_ids, ]
rownames(meta) <- meta$Sample_ID

# --- Step 2: Diagnostics Before Correction ---
# PCA before correction
pca_before <- prcomp(clr, scale. = TRUE)

var_expl <- pca_before$sdev^2 / sum(pca_before$sdev^2)
pc1_var <- round(var_expl[1] * 100, 1)
pc2_var <- round(var_expl[2] * 100, 1)

pca_df <- as.data.frame(pca_before$x)
pca_df$Group <- meta$Group

library(ggplot2)

# Custom colors for environments (change as needed)
env_colors <- c(
  "Ambulance" = "#1f77b4",
  "Hosp_env" = "#2ca02c",
  "Hosp_sewage" = "#ff7f0e",
  "Public_transp" = "#9467bd"
)


p_before <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8, stroke = 0.5) +
  
  # Add ellipses (optional, remove if not desired)
  stat_ellipse(level = 0.95, linewidth = 0.8, alpha = 0.25) +
  
  scale_color_manual(values = env_colors) +
  
  labs(
    title = "",
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)")
  ) +
  
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(color = "black", size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    axis.line = element_line(color = "black", linewidth = 0.7)
  )

ggsave("PCA_Before_BatchCorrection.png",
       plot = p_before, width = 7, height = 6, dpi = 600, bg = "white")

# PERMANOVA before correction
dist_mat <- dist(clr)
adonis_before <- adonis2(dist_mat ~ Group + Instrument + Sequencing_Center + Project_ID, data = meta)
write.csv(as.data.frame(adonis_before), "PERMANOVA_before_correction.csv")

# --- Step 3: Design Matrices ---
design <- model.matrix(~ meta$Group)  # preserve biological signal
covar_v1 <- model.matrix(~ meta$Instrument + meta$Sequencing_Center)[, -1]

# --- Step 4a: Correction Model v1 (Instrument + Sequencing Center) ---
clr_corrected_v1 <- removeBatchEffect(
  as.matrix(t(clr)),        # transpose: features x samples
  covariates = covar_v1,
  design = design
)

clr_corrected_v1 <- t(clr_corrected_v1)  # transpose back: samples x features
write.csv(clr_corrected_v1, "Mechanism_CLR_batch_corrected_v1.tsv", row.names = TRUE)

# --- Step 4b: Correction Model v2 (Instrument + Sequencing Center + Project_ID) ---
clr_corrected_v2 <- removeBatchEffect(
  as.matrix(t(clr)),        # transpose: features x samples
  batch = meta$Project_ID,
  covariates = covar_v1,
  design = design
)

clr_corrected_v2 <- t(clr_corrected_v2)  # transpose back
write.csv(clr_corrected_v2, "Mechanism_CLR_batch_corrected_v2.tsv", row.names = TRUE)

# --- Step 5: Function for Diagnostics After Correction ---
run_diagnostics <- function(corrected_matrix, meta, prefix) {
  counts <- corrected_matrix
  
  
  # -------------------------------------
  # 1. PCA after correction
  # -------------------------------------
  pca <- prcomp(counts, scale. = TRUE)
  
  # Variance explained
  var_expl <- pca$sdev^2 / sum(pca$sdev^2)
  pc1_var <- round(var_expl[1] * 100, 1)
  pc2_var <- round(var_expl[2] * 100, 1)
  
  # -------------------------------------
  # 2. Build PCA dataframe
  # -------------------------------------
  # ------------------------------------------------------------
  # 6. Function for Diagnostics After Correction (Upgraded PCA)
  # ------------------------------------------------------------
  run_diagnostics <- function(corrected_matrix, meta, prefix) {
    
  
    # ===== 1. PCA =====
    pca <- prcomp(corrected_matrix, scale. = TRUE)
    
    # Variance explained
    var_expl <- pca$sdev^2 / sum(pca$sdev^2)
    pc1_var <- round(var_expl[1] * 100, 1)
    pc2_var <- round(var_expl[2] * 100, 1)
    
    # ===== 2. Build dataframe =====
    pca_df <- as.data.frame(pca$x)
    
    # Align PCA rows with metadata sample order
    pca_df$Sample_ID <- rownames(pca_df)
    pca_df$Group <- meta$Group[match(pca_df$Sample_ID, meta$Sample_ID)]
    
    # Ensure consistent factor ordering
    pca_df$Group <- factor(
      pca_df$Group,
      levels = c("Ambulance", "Hosp_env", "Hosp_sewage", "Public_transp")
    )
    
    # ===== 3. Define color palette =====
    
    # ===== 4. PCA plot =====
    p <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
      geom_point(size = 3, alpha = 0.8, stroke = 0.3) +
      stat_ellipse(level = 0.95, linewidth = 0.8, alpha = 0.18) +
      scale_color_manual(values = env_colors) +
      labs(
        title = paste("PCA After Correction:", prefix),
        x = paste0("PC1 (", pc1_var, "% variance)"),
        y = paste0("PC2 (", pc2_var, "% variance)")
      ) +
      theme_classic(base_size = 16) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(color = "black", size = 12),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 12),
        axis.line = element_line(color = "black", linewidth = 0.7)
      )
    
    # Save figure
    ggsave(
      paste0("PCA_After_", prefix, "_Nice.png"),
      plot = p,
      width = 7, height = 6, dpi = 600, bg = "white"
    )
  }
  
  # ------------------------------------------------------------
  # 7. RUN DIAGNOSTICS FOR BOTH MODELS
  # ------------------------------------------------------------
  
  # Model v1 (Instrument + Sequencing_Center corrected)
  run_diagnostics(clr_corrected_v1, meta, prefix = "BatchCorrect_v1")
  
  # Model v2 (Instrument + Sequencing_Center + Project_ID corrected)
  run_diagnostics(clr_corrected_v2, meta, prefix = "BatchCorrect_v2")
  
  
  # PERMANOVA
  dist_mat <- dist(counts)
  adonis_res <- adonis2(dist_mat ~ Group + Instrument + Sequencing_Center + Project_ID, data = meta)
  write.csv(as.data.frame(adonis_res), paste0("PERMANOVA_after_", prefix, ".csv"))
  
}

# --- Step 6: Run Diagnostics for v1 and v2 ---
clr_v1_df <- read.csv("Mechanism_CLR_batch_corrected_v1.tsv", row.names = 1)
clr_v2_df <- read.csv("Mechanism_CLR_batch_corrected_v2.tsv", row.names = 1)

run_diagnostics(clr_v1_df, meta, "v1")
run_diagnostics(clr_v2_df, meta, "v2")

# --- Step 7: Compare Variance Explained (Before vs After) ---
before <- read.csv("PERMANOVA_before_correction.csv")
after_v1 <- read.csv("PERMANOVA_after_v1.csv")
after_v2 <- read.csv("PERMANOVA_after_v2.csv")

before$Condition <- "Before"
after_v1$Condition <- "After_v1"
after_v2$Condition <- "After_v2"

combined <- bind_rows(before, after_v1, after_v2)
combined <- combined %>%
  filter(!is.na(F)) %>%
  select(Factor = X, R2, Condition)

# Plot comparison
tiff("PERMANOVA_R2_Comparison_Mechanisms.tiff", width = 4, height = 3, units = "in", res = 600)



################################################################################

library(ggplot2)

# Correct ordering
combined$Condition <- factor(combined$Condition,
                             levels = c("Before", "After_v1", "After_v2"))

cb_palette <- c(
  "Before"   = "#E69F00",
  "After_v1" = "#009E73",
  "After_v2" = "#56B4E9"
)

ggplot(combined, aes(x = R2, y = Condition, fill = Condition)) +
  geom_col(width = 0.55, color = "black") +
  
  geom_text(aes(label = sprintf("%.3f", R2)),
            hjust = -0.2, size = 1.4, fontface = "bold") +
  
  scale_fill_manual(values = cb_palette) +
  
  scale_x_continuous(limits = c(0, max(combined$R2) * 1.25)) +
  
  labs(
    x = expression(R^2 ~ "(Variance Explained)"),
    y = "",
    title = ""
  ) +
  
  theme_classic(base_size = 5) +
  theme(
    plot.title = element_text( hjust = 0.5, size = 5),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("PERMANOVA_R2_Comparison_THIN.png",
       width = 2, height = 1, dpi = 600, bg = "white")

###########################################################

