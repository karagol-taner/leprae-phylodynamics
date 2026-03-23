library(ggtree)
library(treeio)
library(ggplot2)
library(ape)
library(dplyr)

# ── Paths ─────────────────────────────────────────────────────
tree_file <- "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/lep_MCC_dated.tree"
meta_file <- "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/leprosy_metadata.csv"
out_pdf   <- "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/M_leprae_Publication_Detailed.pdf"

# ── Load & prune ──────────────────────────────────────────────
beast_tree <- read.beast(tree_file)
phy        <- beast_tree@phylo

get_tips <- function(node, phy) {
  if (node <= Ntip(phy)) return(node)
  unlist(lapply(phy$edge[phy$edge[,1]==node, 2], get_tips, phy=phy))
}
root_node   <- Ntip(phy) + 1
children    <- phy$edge[phy$edge[,1] == root_node, 2]
tips1       <- get_tips(children[1], phy)
tips2       <- get_tips(children[2], phy)
out_idx     <- if (length(tips1) <= length(tips2)) tips1 else tips2
tree_pruned <- treeio::drop.tip(beast_tree, phy$tip.label[out_idx])
print(paste("Tips after pruning:", Ntip(tree_pruned@phylo)))

# ── Load metadata ─────────────────────────────────────────────
metadata <- read.csv(meta_file, stringsAsFactors = FALSE)
id_col   <- colnames(metadata)[1]
metadata$CleanSRR <- gsub(".final.bam", "", basename(as.character(metadata[[id_col]])))

# ── Build base tree plot (no labels yet) ──────────────────────
p_base <- ggtree(tree_pruned, mrsd = 2020, linewidth = 0.5)

# ── Extract real tip coordinates from the plot data ───────────
tip_data <- p_base$data %>%
  filter(isTip == TRUE) %>%
  select(label, x, y)

# Clean the label to get SRR ID for matching
tip_data$CleanSRR <- gsub(".final.bam", "", basename(tip_data$label))

# Join metadata
tip_data <- left_join(tip_data, metadata, by = "CleanSRR")

# Build display label
tip_data$DisplayLabel <- ifelse(
  is.na(tip_data$Location) | tip_data$Location == "",
  tip_data$CleanSRR,
  paste(tip_data$CleanSRR, tip_data$Location, sep = " | ")
)

# Host color
if (!"Host" %in% colnames(tip_data)) tip_data$Host <- "Global"
tip_data$Host[is.na(tip_data$Host)] <- "Global"

print(paste("Tips matched to metadata:", sum(!is.na(tip_data$Location))))
print(head(tip_data[, c("CleanSRR", "DisplayLabel", "Host", "x", "y")], 5))

# ── Define label x position (fixed column to the right of 2020) ──
LABEL_X <- 2300  # fixed x where all labels will be aligned

# ── Build final plot ──────────────────────────────────────────
p <- p_base +

  # HPD bars
  geom_range(range = 'height_0.95_HPD', color = 'steelblue', alpha = 0.35, size = 2) +

  # Anchor lines
  geom_vline(xintercept = 1050, linetype = "dotted", color = "#8b0000", linewidth = 0.7) +
  geom_vline(xintercept = 2020, linetype = "dashed", color = "darkgray",  linewidth = 0.5) +

  # Tip points
  geom_point(data = tip_data, aes(x = x, y = y, color = Host), size = 2) +

  # Dotted leader lines from tip to label column
  geom_segment(
    data     = tip_data,
    aes(x = x + 5, xend = LABEL_X - 5, y = y, yend = y, color = Host),
    linetype = "dotted",
    linewidth = 0.3
  ) +

  # Labels at fixed x column
  geom_text(
    data  = tip_data,
    aes(x = LABEL_X, y = y, label = DisplayLabel, color = Host),
    size  = 2.5,
    hjust = 0   # left-align text at LABEL_X
  ) +

  scale_color_manual(
    name   = "Host Species",
    values = c("Red Squirrel" = "#FF0000", "Human" = "#0000FF", "Global" = "#808080"),
    na.value = "#808080"
  ) +

  theme_tree2() +

  labs(
    title    = "Evolutionary Timeline of M. leprae",
    subtitle = "Bayesian Molecular Dating · Red dotted line = 1050 AD Winchester epidemic peak",
    x        = "Calendar Year (AD)"
  ) +

  scale_x_continuous(breaks = seq(800, 2000, by = 200)) +
  coord_cartesian(xlim = c(800, 3200), clip = "off") +

  theme(
    legend.position = "bottom",
    legend.title    = element_text(face = "bold", size = 12),
    legend.text     = element_text(size = 11),
    plot.title      = element_text(face = "bold", size = 16),
    plot.subtitle   = element_text(size = 11, color = "#8b0000"),
    axis.title.x    = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.text.x     = element_text(color = "black", size = 10),
    plot.margin     = margin(t = 10, r = 20, b = 10, l = 10, unit = "mm")
  )

# ── Save ──────────────────────────────────────────────────────
ggsave(out_pdf, p, width = 24, height = 30, limitsize = FALSE)
print(paste("✅ Saved to:", out_pdf))
