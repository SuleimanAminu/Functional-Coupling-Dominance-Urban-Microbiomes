###############################################################################
### Alpha Diversity 
### 
###############################################################################
set.seed(42)
library(vegan)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(FSA)
library(tidyr)

###############################################################################
### LOAD DATA
###############################################################################

tpm <- read.delim("Mechanism_TPM_Matrix.tsv", sep = "\t", row.names = 1, check.names = FALSE)
meta <- read.csv("Merged_metadata.csv", header = TRUE, check.names = FALSE)
rownames(meta) <- meta$Sample_ID

# Align
common <- intersect(rownames(tpm), rownames(meta))
tpm <- tpm[common, ]
meta <- meta[common, ]

cat("Using", nrow(tpm), "samples ×", ncol(tpm), "mechanisms\n")

###############################################################################
### COMPUTE ALPHA DIVERSITY METRICS
###############################################################################

alpha_df <- data.frame(
  Sample_ID = rownames(tpm),
  Richness = specnumber(tpm),
  Shannon  = diversity(tpm, index = "shannon"),
  Simpson  = diversity(tpm, index = "simpson"),
  Evenness = diversity(tpm, index = "shannon") / log(specnumber(tpm))
)

alpha_df <- merge(alpha_df, meta, by = "Sample_ID")

write.csv(alpha_df, "Mechanism_Alpha_Diversity_Table.csv", row.names = FALSE)

###############################################################################
### KRUSKAL–WALLIS TESTS FOR ALL METRICS
###############################################################################

metrics <- c("Richness", "Shannon", "Simpson", "Evenness")

kw_summary <- lapply(metrics, function(m){
  kw <- kruskal.test(alpha_df[[m]] ~ alpha_df$Group)
  data.frame(
    Metric = m,
    Chi_Square = as.numeric(kw$statistic),
    df = as.numeric(kw$parameter),
    P_value = kw$p.value
  )
})

kw_summary <- do.call(rbind, kw_summary)
write.csv(kw_summary, "Supplementary_KruskalWallis_AllMetrics.csv", row.names = FALSE)

###############################################################################
### DUNN POST-HOC TESTS FOR ALL METRICS (BH-adjusted)
###############################################################################

dunn_summary <- lapply(metrics, function(m){
  d <- dunnTest(alpha_df[[m]] ~ alpha_df$Group, method = "bh")$res
  data.frame(
    Metric = m,
    Comparison = d$Comparison,
    Z = d$Z,
    P_uncorrected = d$P.unadj,
    P_adjusted = d$P.adj
  )
})

dunn_summary <- do.call(rbind, dunn_summary)
write.csv(dunn_summary, "Supplementary_DunnTest_AllMetrics.csv", row.names = FALSE)


###############################################################################
### BOXPLOTS 
###############################################################################
### OUTPUT DIRECTORY
###############################################################################
outdir <- "Nature_Boxplots"

dir.create(outdir, showWarnings = FALSE)

###############################################################################
### ORDER + COLORS
###############################################################################

env_order <- c(
  "Hosp_sewage",
  "Public_transp",
  "Hosp_env",
  "Ambulance"
)

env_colors <- c(
  "Ambulance"     = "#1f77b4",
  "Hosp_env"      = "#2ca02c",
  "Hosp_sewage"   = "#ff7f0e",
  "Public_transp" = "#9467bd"
)

###############################################################################
### FUNCTION
###############################################################################

make_boxplot <- function(metric){
  
  cat("\nMaking plot for:", metric, "\n")
  
  ###########################################################################
  ### DATA
  ###########################################################################
  
  df <- alpha_df %>%
    select(Group, value = all_of(metric))
  
  df$Group <- factor(
    df$Group,
    levels = env_order
  )
  
  levels(df$Group) <- c(
    "Hosp_sewage",
    "Public_transp",
    "Hosp_env",
    "Ambulance"
  )
  
  ###########################################################################
  ### KRUSKAL-WALLIS TEST
  ###########################################################################
  
  kw <- kruskal.test(
    value ~ Group,
    data = df
  )
  
  p_label <- paste0(
    "Kruskal-Wallis  p = ",
    formatC(
      kw$p.value,
      format = "e",
      digits = 2
    )
  )
  
  ###########################################################################
  ### PLOT
  ###########################################################################
  
  p <- ggplot(
    df,
    aes(
      x = Group,
      y = value,
      fill = Group
    )
  ) +
    
    #########################################################################
  ### BOXPLOTS ONLY
  #########################################################################
  
  geom_boxplot(
    width = 0.58,
    linewidth = 0.55,
    alpha = 0.90,
    color = "black",
    
    outlier.shape = 16,
    outlier.size = 0.8,
    outlier.alpha = 0.35
  ) +
    
    #########################################################################
  ### COLORS
  #########################################################################
  
  scale_fill_manual(
    values = env_colors
  ) +
    
    #########################################################################
  ### KRUSKAL LABEL
  #########################################################################
  
  annotate(
    "text",
    x = 4,
    y = max(df$value, na.rm = TRUE) * 1.06,
    label = p_label,
    hjust = 1,
    size = 3.2,
    fontface = "italic",
    family = "sans"
  ) +
    
    #########################################################################
  ### LABELS
  #########################################################################
  
  labs(
    x = NULL,
    y = metric
  ) +
    
    #########################################################################
  ### PREVENT LABEL CLIPPING
  #########################################################################
  
  coord_cartesian(
    clip = "off"
  ) +
    
    #########################################################################
  ### THEME
  #########################################################################
  
  theme_classic(
    base_size = 10,
    base_family = "sans"
  ) +
    
    theme(
      
      #######################################################################
      ### REMOVE LEGEND
      #######################################################################
      
      legend.position = "none",
      
      #######################################################################
      ### AXIS TITLES
      #######################################################################
      
      axis.title.y = element_text(
        face = "bold",
        size = 12,
        color = "black",
        margin = margin(r = 8)
      ),
      
      axis.title.x = element_blank(),
      
      #######################################################################
      ### X AXIS LABELS
      #######################################################################
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        face = "plain",
        size = 9,
        color = "black",
        margin = margin(t = 4)
      ),
      
      #######################################################################
      ### Y AXIS LABELS
      #######################################################################
      
      axis.text.y = element_text(
        size = 9,
        color = "black"
      ),
      
      #######################################################################
      ### AXIS LINES
      #######################################################################
      
      axis.line = element_line(
        linewidth = 0.5,
        color = "black"
      ),
      
      #######################################################################
      ### AXIS TICKS
      #######################################################################
      
      axis.ticks = element_line(
        linewidth = 0.5,
        color = "black"
      ),
      
      axis.ticks.length = unit(
        0.15,
        "cm"
      ),
      
      #######################################################################
      ### PLOT MARGINS
      #######################################################################
      
      plot.margin = margin(
        t = 10,
        r = 12,
        b = 10,
        l = 10
      )
    )
  
  ###########################################################################
  ### SHOW PLOT
  ###########################################################################
  
  print(p)
  
  ###########################################################################
  ### SAVE PNG
  ###########################################################################
  
  png_file <- file.path(
    outdir,
    paste0(metric, "_boxplot.png")
  )
  
  ggsave(
    filename = png_file,
    plot = p,
    width = 3.6,
    height = 3.1,
    dpi = 600,
    bg = "white"
  )
  
  ###########################################################################
  ### SAVE TIFF
  ###########################################################################
  
  tiff_file <- file.path(
    outdir,
    paste0(metric, "_boxplot.tiff")
  )
  
  ggsave(
    filename = tiff_file,
    plot = p,
    width = 3.6,
    height = 3.1,
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
  
  return(p)
}

###############################################################################
### RUN
###############################################################################

metrics <- c(
  "Richness",
  "Evenness",
  "Shannon",
  "Simpson"
)

plots <- lapply(
  metrics,
  make_boxplot
)

###############################################################################
### CHECK FILES
###############################################################################

list.files(outdir)

###############################################################################
###  NOTE: We replaced the violin plot with Boxplot (As shown in the above code)
###############################################################################
### VIOLIN PLOTS
###############################################################################

env_colors <- c(
  "Ambulance" = "#1f77b4",
  "Hosp_env" = "#2ca02c",
  "Hosp_sewage" = "#ff7f0e",
  "Public_transp" = "#9467bd"
)

make_violin_plot <- function(metric){
  
  df <- alpha_df %>% select(Group, value = all_of(metric))
  
  kw <- kruskal.test(value ~ Group, data = df)
  kw_p <- kw$p.value
  
  # Italic red Kruskal label
  # 
  kw_label <- paste0(
    "bolditalic('Kruskal–Wallis,  p = ",
    formatC(kw_p, format='e', digits=2),
    "')"
  )
  
  p <- ggplot(df, aes(y = Group, x = value, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.85, color = "black", linewidth = 0.34) +
    geom_boxplot(width = 0.12, fill = "white", color = "black",
                 outlier.shape = NA, linewidth = 0.24) +
    scale_fill_manual(values = env_colors) +
    
    annotate(
      "text",
      x = min(df$value),
      y = 0.2,
      label = kw_label,
      parse = TRUE,
      color = "red",
      size = 1.45,
      hjust = -0.08,
      vjust = -0.2
    ) +
    
    labs(x = metric, y = "") +
    theme_classic(base_size = 16)  +
    theme(
      legend.position = "none",
      axis.title.x = element_text(face = "bold", size = 6),
      axis.title.y = element_blank(),
      axis.text.y  = element_text(face = "bold", size = 5),
      axis.text.x  = element_text(size = 6),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.35)
    )
  
  # Save plot
  ggsave(
    paste0("Violin_", metric, "_KWonly_600dpi.tiff"),
    plot = p, width = 2, height = 1.5, dpi = 600, bg = "white"
  )
  
  return(p)
}

plots <- lapply(metrics, make_violin_plot)

