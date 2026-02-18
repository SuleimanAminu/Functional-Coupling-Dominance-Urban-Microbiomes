###############################################################################
###  Mechanism-level Differential Abundance 
### Method: limma on batch-corrected CLR
###############################################################################

suppressPackageStartupMessages({
  library(limma)
  library(tidyverse)
})

set.seed(42)

###############################################################################
# 1. LOAD CLR MATRIX (CSV FORMAT, ROW NAMES = SAMPLE IDs)
###############################################################################

clr_raw <- read.csv(
  "Mechanism_CLR_batch_corrected_v1.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1
)

# Keep only numeric mechanism columns (drops any stray index column)
clr <- clr_raw %>% select(where(is.numeric))

cat(
  "CLR matrix:",
  nrow(clr), "samples ×", ncol(clr), "mechanisms\n"
)

stopifnot(ncol(clr) == 23)

# Convert to numeric matrix
expr <- as.matrix(clr)
mode(expr) <- "numeric"

###############################################################################
# 2. LOAD & ALIGN METADATA
###############################################################################

meta <- read.csv("Merged_metadata.csv", check.names = FALSE)
meta$Sample_ID <- trimws(meta$Sample_ID)
rownames(meta) <- meta$Sample_ID

common_ids <- intersect(rownames(expr), meta$Sample_ID)
expr <- expr[common_ids, , drop = FALSE]
meta <- meta[common_ids, ]
rownames(meta) <- meta$Sample_ID

stopifnot(nrow(expr) == nrow(meta))

###############################################################################
# 3. DEFINE GROUP FACTOR
###############################################################################

meta$Group <- factor(
  meta$Group,
  levels = c("Ambulance", "Hosp_env", "Hosp_sewage", "Public_transp")
)

###############################################################################
# 4. DESIGN MATRIX (NO INTERCEPT)
###############################################################################

design <- model.matrix(~ 0 + Group, data = meta)
colnames(design) <- levels(meta$Group)

###############################################################################
# 5. FIT LIMMA MODEL (FEATURES × SAMPLES)
###############################################################################

fit <- lmFit(t(expr), design)

cat(
  "limma fit:",
  nrow(t(expr)), "features ×", ncol(t(expr)), "samples\n"
)

###############################################################################
# 6. DEFINE CONTRASTS
###############################################################################

contrast_matrix <- makeContrasts(
  Ambulance_vs_HospEnv     = Ambulance - Hosp_env,
  Ambulance_vs_Sewage     = Ambulance - Hosp_sewage,
  Ambulance_vs_Public     = Ambulance - Public_transp,
  HospEnv_vs_Sewage       = Hosp_env - Hosp_sewage,
  HospEnv_vs_Public       = Hosp_env - Public_transp,
  Sewage_vs_Public        = Hosp_sewage - Public_transp,
  levels = design
)

###############################################################################
# 7. APPLY CONTRASTS + EMPIRICAL BAYES
###############################################################################

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

###############################################################################
# 8. EXTRACT DIFFERENTIAL ABUNDANCE RESULTS
###############################################################################

da_results <- lapply(colnames(contrast_matrix), function(cn) {
  topTable(
    fit2,
    coef = cn,
    number = Inf,
    adjust.method = "BH"
  ) %>%
    rownames_to_column("Mechanism") %>%
    mutate(
      Contrast  = cn,
      Direction = case_when(
        logFC > 0 & adj.P.Val < 0.05 ~ "Up",
        logFC < 0 & adj.P.Val < 0.05 ~ "Down",
        TRUE                        ~ "NS"
      )
    )
}) %>% bind_rows()

write.csv(
  da_results,
  "Mechanism_Differential_Abundance_limma.csv",
  row.names = FALSE
)

###############################################################################
# 9. DIRECTIONAL SUMMARY (FOR DBI / FID)
###############################################################################

da_direction_summary <- da_results %>%
  filter(Direction != "NS") %>%
  separate(Contrast, into = c("Env1", "Env2"), sep = "_vs_") %>%
  mutate(Environment = ifelse(Direction == "Up", Env1, Env2)) %>%
  group_by(Environment, Direction) %>%
  summarise(n_mechanisms = n(), .groups = "drop")

write.csv(
  da_direction_summary,
  "Mechanism_Directional_Enrichment_by_Environment.csv",
  row.names = FALSE
)

##########################################################################################
####### Environment_specific Functions
###################################################################################



library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(grid)

# -------------------------------
# INPUT MATRIX
# -------------------------------
mat <- as.matrix(env_mat)

# -------------------------------
# FORCE CORRECT COLUMN NAMES + ORDER
# -------------------------------
# Rename columns explicitly (DO NOT rely on existing names)
colnames(mat) <- c("Ambulance", "Hosp_env", "Public_transp", "Hosp_sewage")

# Enforce order
mat <- mat[, c("Ambulance", "Hosp_env", "Public_transp", "Hosp_sewage")]

# -------------------------------
# DEFINE FUNCTIONAL CLASSES
# -------------------------------
arg_mechs <- c("RedPerm","TargetAlt","Antibiotic efflux",
               "TargetRepl","TargetProtect","Antibiotic inactivation")

hybrid_mechs <- c("Inact + TargetAlt","TargetAlt + TargetRepl",
                  "Efflux + RedPerm","Efflux + TargetAlt")

vf_mechs <- c("Adherence","Biofilm","CompetitiveAdv","Exotoxin",
              "Immune modulation","Invasion","Motility","Others",
              "Regulation","Stress survival",
              "Effector delivery system","Nutritional Metabolic factor",
              "Exoenzyme")

mech_class <- case_when(
  rownames(mat) %in% arg_mechs    ~ "ARG",
  rownames(mat) %in% hybrid_mechs ~ "Hybrid",
  rownames(mat) %in% vf_mechs     ~ "VF",
  TRUE                            ~ "Other"
)

row_split <- factor(mech_class, levels = c("ARG", "Hybrid", "VF"))

# Order rows by class only (no clustering chaos)
row_order <- order(row_split)
mat <- mat[row_order, ]
row_split <- row_split[row_order]

# -------------------------------
# HEATMAP COLOR (REFERENCE-LIKE RED PALETTE)
# -------------------------------
col_fun <- colorRamp2(
  c(0, max(mat) * 0.4, max(mat)),
  c("#f7fbff", "#c6dbef", "#2171b5")
)

# -------------------------------
# ROW CLASS STRIP
# -------------------------------
row_ha <- rowAnnotation(
  Class = row_split,
  col = list(
    Class = c(
      "ARG"    = "#2166ac",
      "Hybrid" = "#7C3AED",
      "VF"     = "#CB181D"
    )
  ),
  annotation_width = unit(4, "mm"),
  gp = gpar(col = NA)
)

# -------------------------------
# HEATMAP (CLEAN, NO INTERNAL GRID)
# -------------------------------
ht <- Heatmap(
  mat,
  name = "Mean |logFC|",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_split,
  
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 45,
  
  rect_gp = gpar(col = NA),      # REMOVE internal borders
  border = TRUE,                 # KEEP outer border
  border_gp = gpar(col = "black", lwd = 0.6),
  
  left_annotation = row_ha,
  
  heatmap_legend_param = list(
    title = "Mean |logFC|",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    grid_width = unit(4, "mm")
  )
)

# -------------------------------
# SAVE FIGURE (NATURE QUALITY)
# -------------------------------
tiff(
  "Figure_Functional_Mechanism_Heatmap_FINAL_Nature.tiff",
  width = 7,
  height = 5,
  units = "in",
  res = 600
)
draw(ht)
dev.off()
##############################################################

#### Pairwise

#############################################

library(tidyverse)
library(ggplot2)

# ============================================================
# INPUT
# ============================================================
da <- read.csv(
  "Mechanism_Differential_Abundance_limma.csv",
  check.names = FALSE
)

# ============================================================
# CANONICAL ENVIRONMENT NAME MAP (CRITICAL FIX)
# ============================================================
env_name_map <- c(
  "Ambulance" = "Ambulance",
  "HospEnv" = "Hosp_env",
  "Hosp_env" = "Hosp_env",
  "Public" = "Public_transp",
  "Public_transp" = "Public_transp",
  "Sewage" = "Hosp_sewage",
  "Hosp_sewage" = "Hosp_sewage"
)

# ============================================================
# MECHANISM RENAMING (CONSISTENT WITH ALL PREVIOUS FIGURES)
# ============================================================
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

# ============================================================
# FUNCTIONAL CLASSES
# ============================================================
arg_mechs <- c(
  "Antibiotic efflux","Antibiotic inactivation",
  "TargetAlt","TargetRepl","TargetProtect","RedPerm"
)

hybrid_mechs <- c(
  "Efflux + TargetAlt","Efflux + RedPerm",
  "Inact + TargetAlt","TargetAlt + TargetRepl"
)

vf_mechs <- c(
  "Adherence","Biofilm","Effector delivery system",
  "Exoenzyme","Exotoxin","Immune modulation",
  "Invasion","Motility","Nutritional Metabolic factor",
  "CompetitiveAdv","Stress survival","Regulation","Others"
)

# ============================================================
# PREPARE DATA (WITH ENV NAME NORMALIZATION)
# ============================================================
df <- da %>%
  filter(adj.P.Val < 0.05) %>%
  separate(Contrast, into = c("Env1_raw", "Env2_raw"), sep = "_vs_") %>%
  mutate(
    Env1 = env_name_map[Env1_raw],
    Env2 = env_name_map[Env2_raw],
    Contrast_clean = paste0(Env1, "_vs_", Env2),
    Mechanism_short = mechanism_short_labels[Mechanism],
    Class = case_when(
      Mechanism_short %in% arg_mechs    ~ "ARG",
      Mechanism_short %in% hybrid_mechs ~ "Hybrid",
      Mechanism_short %in% vf_mechs     ~ "VF",
      TRUE                              ~ NA_character_
    )
  ) %>%
  drop_na(Mechanism_short, Class, Env1, Env2)

# ============================================================
# COLOR SCHEME (NATURE-STYLE)
# ============================================================
class_colors <- c(
  "ARG"    = "#2166AC",
  "Hybrid" = "#7C3AED",
  "VF"     = "#B2182B"
)

# ============================================================
# OUTPUT DIRECTORY
# ============================================================
dir.create("Forest_Plots_By_Contrast", showWarnings = FALSE)

# ============================================================
# AUTOMATED FOREST PLOTS (ALL CLEAN CONTRASTS)
# ============================================================
for (cn in unique(df$Contrast_clean)) {
  
  sub <- df %>%
    filter(Contrast_clean == cn) %>%
    arrange(logFC)
  
  if (nrow(sub) == 0) next
  
  sub$Mechanism_short <- factor(
    sub$Mechanism_short,
    levels = sub$Mechanism_short
  )
  
  envs <- strsplit(cn, "_vs_")[[1]]
  
  p <- ggplot(sub, aes(x = logFC, y = Mechanism_short)) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.6
    ) +
    
    geom_point(
      aes(color = Class),
      size = 3.4
    ) +
    
    scale_color_manual(
      values = class_colors,
      name = "Mechanism class"
    ) +
    
    labs(
      title = paste(envs[1], "vs", envs[2]),
      x = expression(log[2]~fold~change),
      y = NULL
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 11),
      axis.title.x = element_text(size = 12),
      axis.line = element_line(linewidth = 0.7),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 8),
      legend.text  = element_text(size = 7),
      panel.grid = element_blank()
    )
  
  ggsave(
    filename = file.path(
      "Forest_Plots_By_Contrast",
      paste0("Forest_", cn, "_Mechanisms_ClassColored.tiff")
    ),
    plot = p,
    width = 5.4,
    height = max(4, 0.28 * nrow(sub)),
    dpi = 600,
    bg = "white"
  )
}

##################################################################
########### Representative gene extraction per enriched mechanism
###########################################################
library(tidyverse)

# ============================================================
# 1. LOAD MECHANISM-LEVEL DA (ANCHOR)
# ============================================================
da <- read.csv(
  "Mechanism_Differential_Abundance_limma.csv",
  check.names = FALSE
)

# Canonical environment name map (SINGLE SOURCE OF TRUTH)
env_name_map <- c(
  "Ambulance" = "Ambulance",
  "HospEnv" = "Hosp_env",
  "Hosp_env" = "Hosp_env",
  "Public" = "Public_transp",
  "Public_transp" = "Public_transp",
  "Sewage" = "Hosp_sewage",
  "Hosp_sewage" = "Hosp_sewage"
)

# Extract enriched mechanisms per environment (canonical names)
sig_mechs <- da %>%
  filter(adj.P.Val < 0.05) %>%
  separate(Contrast, into = c("Env1_raw", "Env2_raw"), sep = "_vs_") %>%
  mutate(
    Env1 = env_name_map[Env1_raw],
    Env2 = env_name_map[Env2_raw],
    Enriched_env = ifelse(logFC > 0, Env1, Env2)
  ) %>%
  select(Mechanism, Enriched_env) %>%
  distinct()

# ============================================================
# 2. LOAD HARMONIZED GENE TABLE
# ============================================================
genes <- read.delim(
  "Combined_ARG_VFDB_Harmonized.tsv",
  check.names = FALSE
)

stopifnot(all(c(
  "Gene_ID", "Mechanisms", "TPM", "Sample_ID", "Group", "Species_name"
) %in% colnames(genes)))

# ============================================================
# 3. FILTER GENES UNDER ENRICHED MECHANISMS
#    (ALLOW MANY-TO-MANY — THIS IS BIOLOGICALLY EXPECTED)
# ============================================================
genes_sig <- genes %>%
  inner_join(
    sig_mechs,
    by = c("Mechanisms" = "Mechanism"),
    relationship = "many-to-many"
  ) %>%
  filter(Group == Enriched_env)

# ============================================================
# 4. GENE-LEVEL SUMMARY (DESCRIPTIVE ONLY)
# ============================================================
gene_summary <- genes_sig %>%
  group_by(
    Enriched_env,
    Mechanisms,
    Gene_ID
  ) %>%
  summarise(
    Total_TPM  = sum(TPM, na.rm = TRUE),
    Prevalence = mean(TPM > 0),
    .groups = "drop"
  )

# ============================================================
# 5. SELECT REPRESENTATIVE GENES (TOP 5)
# ============================================================
rep_genes <- gene_summary %>%
  arrange(
    Enriched_env,
    Mechanisms,
    desc(Total_TPM),
    desc(Prevalence)
  ) %>%
  group_by(Enriched_env, Mechanisms) %>%
  slice_head(n = 5) %>%
  ungroup()

# ============================================================
# 6. SAVE
# ============================================================
write.csv(
  rep_genes,
  "Representative_Genes_by_Mechanism_and_Environment.csv",
  row.names = FALSE
)

###############################

###############################################################################
###############################################################################
### GRUMB–FPP | Genes & Species supporting environment-structured mechanisms
### FINAL — naming fixed, ALL environments retained
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

set.seed(42)

###############################################################################
# 1. LOAD MECHANISM-LEVEL DA (ANCHOR ONLY)
###############################################################################

da <- read.csv(
  "Mechanism_Differential_Abundance_limma.csv",
  check.names = FALSE
)

sig_mechanisms <- da %>%
  filter(adj.P.Val < 0.05) %>%
  select(Mechanism) %>%
  distinct()

###############################################################################
# 2. LOAD HARMONIZED GENE TABLE
###############################################################################

genes <- read.delim(
  "Combined_ARG_VFDB_Harmonized.tsv",
  check.names = FALSE
)

# ---- CRITICAL FIX: harmonise environment names ----
genes <- genes %>%
  mutate(
    Group = case_when(
      Group == "Public" ~ "Public_transp",
      TRUE              ~ Group
    )
  )

# Sanity check (must show 4 groups)
print(unique(genes$Group))

stopifnot(all(c(
  "Gene_ID",
  "Mechanisms",
  "TPM",
  "Sample_ID",
  "Group",
  "Species_name"
) %in% colnames(genes)))

###############################################################################
# 3. KEEP GENES BELONGING TO SIGNIFICANT MECHANISMS
###############################################################################

genes_filt <- genes %>%
  inner_join(
    sig_mechanisms,
    by = c("Mechanisms" = "Mechanism")
  )

###############################################################################
# 4. SUMMARISE GENE SUPPORT PER ENV × MECHANISM × SPECIES
###############################################################################

gene_summary <- genes_filt %>%
  group_by(
    Group,
    Mechanisms,
    Gene_ID,
    Species_name
  ) %>%
  summarise(
    Mean_TPM   = mean(TPM, na.rm = TRUE),
    Prevalence = mean(TPM > 0),
    .groups = "drop"
  )

###############################################################################
# 5. SELECT REPRESENTATIVE GENES (TOP 2 PER ENV × MECHANISM)
###############################################################################

rep_genes <- gene_summary %>%
  arrange(
    Group,
    Mechanisms,
    desc(Mean_TPM),
    desc(Prevalence)
  ) %>%
  group_by(Group, Mechanisms) %>%
  slice_head(n = 3) %>%
  ungroup()

###############################################################################
# 6. SAVE — MASTER TABLE (GENES + SPECIES + ENV)
###############################################################################

write.csv(
  rep_genes,
  "Representative_Genes_and_Species_by_Environment_and_Mechanism.csv",
  row.names = FALSE
)

cat("Representative genes and species extracted for ALL environments.\n")
###############################################################################


