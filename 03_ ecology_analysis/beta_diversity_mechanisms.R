###############################################################################
### Beta Diversity (Aitchison + PERMANOVA + CAP)
###############################################################################
set.seed(42)

library(vegan)
library(ggplot2)
library(dplyr)

###############################################################################
### 1. LOAD DATA
###############################################################################

meta <- read.csv("Merged_metadata.csv", check.names = FALSE)
rownames(meta) <- meta$Sample_ID

clr_v1 <- read.csv("Mechanism_CLR_batch_corrected_v1.tsv",
                   row.names = 1, check.names = FALSE)

common <- intersect(rownames(clr_v1), rownames(meta))
clr_v1 <- clr_v1[common, ]
meta   <- meta[common, ]

meta$Group <- factor(meta$Group)

cat("Loaded:", nrow(clr_v1), "samples ×", ncol(clr_v1), "mechanisms\n")

###############################################################################
### 2. AITCHISON DISTANCE
###############################################################################

dist_v1 <- dist(clr_v1, method = "euclidean")
cat("Aitchison distance computed.\n")

###############################################################################
### 3. GLOBAL PERMANOVA
###############################################################################

adonis_v1 <- adonis2(dist_v1 ~ Group, data = meta, permutations = 999)
write.csv(as.data.frame(adonis_v1), "PERMANOVA_Aitchison_v1.csv", row.names = FALSE)
print(adonis_v1)

###############################################################################
### 4. PAIRWISE PERMANOVA 
###############################################################################

pairwise_permanova <- function(dist, group, p.adjust.method = "BH") {
  
  group <- factor(group)
  groups <- levels(group)
  results <- data.frame()
  
  for(i in 1:(length(groups)-1)){
    for(j in (i+1):length(groups)){
      
      g1 <- groups[i]
      g2 <- groups[j]
      
      idx <- group %in% c(g1, g2)
      dist_sub  <- as.dist(as.matrix(dist)[idx, idx])
      group_sub <- droplevels(group[idx])
      
      ad <- adonis2(dist_sub ~ group_sub, permutations = 999)
      
      results <- rbind(results, data.frame(
        Group1 = g1,
        Group2 = g2,
        R2      = ad$R2[1],
        p_value = ad$`Pr(>F)`[1]
      ))
    }
  }
  
  results$p_adj <- p.adjust(results$p_value, method = p.adjust.method)
  return(results)
}

pairwise_res <- pairwise_permanova(dist_v1, meta$Group)
write.csv(pairwise_res, "Pairwise_PERMANOVA_Aitchison_v1.csv", row.names = FALSE)
print(pairwise_res)

###############################################################################
### 5. BETADISPER (distance-to-centroid)
###############################################################################

bd <- betadisper(dist_v1, meta$Group)
anova_bd <- anova(bd)
write.csv(as.data.frame(anova_bd), "Betadisper_Aitchison_v1.csv")

centroid_df <- data.frame(
  Sample_ID = names(bd$distances),
  Distance  = bd$distances,
  Group     = meta$Group
)

write.csv(centroid_df, "Distance_to_Centroid_Aitchison_v1.csv", row.names = FALSE)

###############################################################################
### 6. CAP — Constrained PCoA 
###############################################################################

cap <- capscale(dist_v1 ~ Group, data = meta)
anova_cap <- anova(cap, permutations = 999)
write.csv(as.data.frame(anova_cap), "CAP_ANOVA_v1.csv")

# Extract site scores
cap_scores <- scores(cap, display = "sites", choices = 1:2)
cap_df <- data.frame(
  CAP1  = cap_scores[,1],
  CAP2  = cap_scores[,2],
  Group = meta$Group
)

# Variance explained
eig_vals <- cap$CCA$eig
cap1_var <- round(100 * eig_vals[1] / sum(eig_vals), 2)
cap2_var <- round(100 * eig_vals[2] / sum(eig_vals), 2)
###############################################################################
### 7. INTERNAL GROUP CENTROIDS (for inside labels)
###############################################################################

centroids <- cap_df %>%
  group_by(Group) %>%
  summarise(
    cx = median(CAP1),
    cy = median(CAP2)
  )
################################################################################
### 8. PLOTTING 
###############################################################################

group_colors <- c(
  "Ambulance"     = "#1f77b4",
  "Hosp_env"      = "#2ca02c",
  "Hosp_sewage"   = "#ff7f0e",
  "Public_transp" = "#9467bd"
)

cap_plot <- ggplot(cap_df, aes(CAP1, CAP2, color = Group)) +
  
  geom_point(size = 1.6, alpha = 0.85) +
  stat_ellipse(type = "norm", linewidth = 0.5, linetype = 2) +
  scale_color_manual(values = group_colors) +
  
  # Add internal text labels
  geom_text(data = centroids,
            aes(x = cx, y = cy, label = Group),
            color = "black",
            fontface = "bold",
            size = 2) +
  
  theme_minimal(base_size = 11) +
  labs(
    title = "",
    x = paste0("CAP1 (", cap1_var, "%)"),
    y = paste0("CAP2 (", cap2_var, "%)")
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    panel.grid   = element_blank(),          
    plot.title   = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black"), 
    axis.title   = element_text(face = "bold")
  )

tiff("Aitchison_CAP_v1_INTERNAL_LABELS_600dpi.tiff",
     width = 5, height = 4, units = "in", res = 600)
print(cap_plot)
dev.off()

###############################################################################
## Distance-to-Centroid Boxplots (Aitchison Beta Dispersion)
###############################################################################

library(ggplot2)
library(ggpubr)
library(FSA)
library(dplyr)

###############################################################################
### 1. EXTRACT DISTANCES FROM BETADISPER
###############################################################################

bd_v1 <- betadisper(dist_v1, meta$Group)

centroid_df <- data.frame(
  Sample_ID = names(bd_v1$distances),
  Distance  = bd_v1$distances,
  Group     = meta$Group
)

write.csv(centroid_df, "Distance_to_Centroid_Aitchison_v1.csv", row.names = FALSE)

###############################################################################
### 2. STATISTICAL TESTING — KW + Dunn
###############################################################################

# Global test
kw_dist <- kruskal.test(Distance ~ Group, data = centroid_df)
kw_p <- kw_dist$p.value

# Pairwise posthoc
dunn_dist_res <- dunnTest(Distance ~ Group, data = centroid_df, method = "bh")$res
write.csv(dunn_dist_res, "Dunn_Distance_to_Centroid.csv", row.names = FALSE)

kw_label <- paste0(
  "italic('Kruskal–Wallis  p = ", 
  formatC(kw_p, format='e', digits=2),
  "')"
)

###############################################################################
### Distance-to-Centroid Boxplot 
###############################################################################

library(ggplot2)
library(ggpubr)
library(dplyr)

env_cols <- c(
  "Ambulance"     = "#1f77b4",
  "Hosp_env"      = "#2ca02c",
  "Hosp_sewage"   = "#ff7f0e",
  "Public_transp" = "#9467bd"
)

# Kruskal–Wallis p-label
kw <- kruskal.test(Distance ~ Group, data = centroid_df)
kw_label <- paste0(
  "bolditalic('Kruskal–Wallis:  p = ",
  formatC(kw$p.value, format='e', digits=2),
  "')"
)

p_centroid <- ggplot(centroid_df, aes(x = Group, y = Distance, fill = Group)) +
  
  # Jitter points (behind the boxplot)
  #geom_jitter(width = 0.05, size = 0.05, alpha = 0.35, color = "black") +
  
  # Clean boxplot
  geom_boxplot(
    width = 0.6,
    alpha = 0.85,
    color = "black",
    linewidth = 0.4,
    outlier.size = 1.1,
    outlier.alpha = 0.5
  ) +
  
  scale_fill_manual(values = env_cols) +
  
  # Kruskal annotation
  annotate(
    "text",
    x = 0.7,
    y = max(centroid_df$Distance) * 1.05,
    label = kw_label,
    parse = TRUE,
    color = "red",
    size = 0.73,
    hjust = 0
  ) +
  
  labs(
    x = "",
    y = "Distance to Group Centroid"
  ) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    
    # Uniform thickness → remove axis.line
    axis.line = element_blank(),
    
    # Keep border only
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    
    # X-axis
    axis.text.x = element_text(
      face = "bold",
      size = 3,
      angle = 90,
      hjust = 1,
      color = "black"
    ),
    
    # Y-axis
    axis.text.y = element_text(
      size = 3,
      color = "black",
      face = "bold"
    ),
    
    axis.title.y = element_text(face = "bold", size = 3),
    
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.3)
  )

# Save
ggsave(
  "Distance_to_Centroid_Boxplot_with_Jitter_600dpi.tiff",
  plot = p_centroid,
  width = 1.0, height = 1.9, dpi = 600, bg = "white"
)

p_centroid

