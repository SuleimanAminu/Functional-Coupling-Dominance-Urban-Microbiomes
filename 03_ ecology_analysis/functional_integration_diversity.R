###############################################################################
###  Functional Integration Diversity (FID)
### ARG–VF coupling analysis (sample-level)
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggplot2)
  library(FSA)
})

set.seed(45)

# -------------------------------------------------------------------------
# 1. LOAD METADATA
# -------------------------------------------------------------------------
meta <- read.csv("Merged_metadata.csv", check.names = FALSE)
meta$Sample_ID <- trimws(meta$Sample_ID)
rownames(meta) <- meta$Sample_ID

# -------------------------------------------------------------------------
# 2. LOAD ARG + VF MECHANISM TABLES (GENE-LEVEL)
# -------------------------------------------------------------------------
ARG <- read.delim("ARG_mechanism.tsv", check.names = FALSE)
VF  <- read.delim("VF_mechanism.tsv",  check.names = FALSE)

ARG$Sample_ID <- trimws(ARG$Sample_ID)
VF$Sample_ID  <- trimws(VF$Sample_ID)

# -------------------------------------------------------------------------
# 3. COLLAPSE GENES → MECHANISMS PER SAMPLE
# -------------------------------------------------------------------------
ARG_mech <- ARG %>%
  group_by(Sample_ID, ARG_Mechanism) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = ARG_Mechanism,
              values_from = TPM,
              values_fill = 0)

VF_mech <- VF %>%
  group_by(Sample_ID, VF_Mechanism) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = VF_Mechanism,
              values_from = TPM,
              values_fill = 0)

# -------------------------------------------------------------------------
# 4. ALIGN SAMPLES
# -------------------------------------------------------------------------
common_ids <- Reduce(intersect, list(
  ARG_mech$Sample_ID,
  VF_mech$Sample_ID,
  meta$Sample_ID
))

ARG_mech <- ARG_mech %>% filter(Sample_ID %in% common_ids)
VF_mech  <- VF_mech  %>% filter(Sample_ID %in% common_ids)
meta     <- meta     %>% filter(Sample_ID %in% common_ids)

ARG_mat <- ARG_mech %>% column_to_rownames("Sample_ID")
VF_mat  <- VF_mech  %>% column_to_rownames("Sample_ID")

# -------------------------------------------------------------------------
# 5. FUNCTIONAL ENTROPY CALCULATIONS
# -------------------------------------------------------------------------
H_ARG <- diversity(ARG_mat, index = "shannon")
H_VF  <- diversity(VF_mat,  index = "shannon")
H_TOT <- diversity(cbind(ARG_mat, VF_mat), index = "shannon")

J_ARG <- H_ARG / log(specnumber(ARG_mat))
J_VF  <- H_VF  / log(specnumber(VF_mat))

FID <- data.frame(
  Sample_ID = names(H_ARG),
  H_ARG,
  H_VF,
  H_TOT,
  J_ARG,
  J_VF,
  FCR = log2(H_ARG / H_VF),
  Integration_Index = (H_TOT - abs(H_ARG - H_VF)) / H_TOT
)

FID <- FID %>% left_join(meta, by = "Sample_ID")

write.csv(FID, "Functional_Integration_Diversity.csv", row.names = FALSE)

# -------------------------------------------------------------------------
# 6. NATURE-STYLE VISUALIZATION — INTEGRATION SPACE
# -------------------------------------------------------------------------
group_colors <- c(
  "Ambulance"     = "#1f78b4",
  "Hosp_env"      = "#33a02c",
  "Hosp_sewage"   = "#ff7f00",
  "Public_transp" = "#6a3d9a"
)

p1 <- ggplot(
  FID,
  aes(
    x = H_TOT,
    y = Integration_Index,
    color = Group
  )
) +
  
  # Points: slightly lighter so ellipses don’t overpower
  geom_point(
    size = 1.3,
    alpha = 0.8
  ) +
  
  # Ellipses: dashed, thinner → contextual, not dominant
  stat_ellipse(
    type = "norm",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  
  scale_color_manual(values = group_colors) +
  
  labs(
    x = "Total functional entropy (ARG + VF)",
    y = "ARG–VF integration index"
  ) +
  
  # Clean Nature-style theme
  theme_classic(base_size = 10) +
  
  theme(
    # Single outer border only
    panel.border = element_rect(
      color = "black",
      fill  = NA,
      linewidth = 1.0
    ),
    
    # REMOVE axis lines to avoid double-border effect
    axis.line = element_blank(),
    
    # Typography
    axis.title = element_text(
      face = "bold",
      size = 7
    ),
    axis.text = element_text(
      size = 6,
      color = "black"
    ),
    
    # Legend: clean and unobtrusive
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    
    # No grids
    panel.grid = element_blank(),
    
    # Tight margins (journal look)
    plot.margin = margin(6, 6, 6, 6)
  )

ggsave(
  "FID_Integration_Space.tiff",
  p1,
  width  = 4.0,
  height = 3.0,
  units  = "in",
  dpi    = 600,
  bg     = "white"
)



# -------------------------------------------------------------------------
# 7. FUNCTIONAL COUPLING RATIO (ARG vs VF)
# -------------------------------------------------------------------------
p2 <- ggplot(FID, aes(x = Group, y = FCR, fill = Group)) +
  
  geom_violin(
    trim = FALSE,
    alpha = 0.75,
    color = "black",
    linewidth = 0.2
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.1,
    color = "red"
  ) +
  
  stat_summary(
    fun = median,
    geom = "point",
    size = 1.0,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.6
  ) +
  
  scale_fill_manual(values = group_colors) +
  
  ## CONTROL Y TICKS
  scale_y_continuous(breaks = c(-2, 0, 2, 4)) +
  
  labs(
    y = bquote(bold("FCR (" * log[2] * " ARG / VF)")),
    x = NULL
  ) +
  
  theme_classic(base_size = 10) +
  
  theme(
    # Outer border only
    panel.border = element_rect(
      color = "black",
      fill  = NA,
      linewidth = 0.5
    ),
    
    axis.line = element_blank(),
    
    # 🔹 AXIS TITLES
    axis.title.y = element_text(
      face = "bold",
      size = 3
    ),
    
    # 🔹 AXIS TICK LABELS (THIS IS WHAT YOU WANT)
    axis.text.x = element_text(
      face = "bold",
      size = 3,
      angle = 90,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 3,
      color = "black"
    ),
    
    # No legend
    legend.position = "none",
    
    panel.grid = element_blank(),
    plot.margin = margin(6, 6, 6, 6)
  )



ggsave(
  "FID_Functional_Coupling.tiff",
  p2,
  width  = 1.0,
  height = 1.1,
  units  = "in",
  dpi    = 600,
  bg     = "white"
)


# -------------------------------------------------------------------------
# 8. STATISTICS
# -------------------------------------------------------------------------
kw_int <- kruskal.test(Integration_Index ~ Group, data = FID)
kw_fcr <- kruskal.test(FCR ~ Group, data = FID)

dunn_int <- dunnTest(Integration_Index ~ Group, data = FID, method = "bh")$res
dunn_fcr <- dunnTest(FCR ~ Group, data = FID, method = "bh")$res

write.csv(dunn_int, "Dunn_Integration_Index.csv", row.names = FALSE)
write.csv(dunn_fcr, "Dunn_Functional_Coupling.csv", row.names = FALSE)

###############################################################################


































































