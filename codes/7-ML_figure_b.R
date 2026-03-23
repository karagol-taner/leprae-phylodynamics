library(ggtree)
library(treeio)
library(ggplot2)
library(ape)

print("--- Constructing Upgraded Panel A: Macro Evolutionary Divergence ---")

tree_file <- "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/lep_phylogeny_final.treefile"
out_pdf_A <- "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/ggtree_panel_A_macro_improved.pdf"

# 1. Load the full, unpruned tree
tree <- read.tree(tree_file)
tree$tip.label <- gsub(".final.bam", "", tree$tip.label)
tree$tip.label <- sapply(strsplit(tree$tip.label, "/"), tail, 1)

N <- length(tree$tip.label)
dist_to_root <- dist.nodes(tree)[N + 1, 1:N]

# 2. Dynamically classify every single tip (Host + Species)
classify_tip <- function(node_name, dist) {
  if (dist > 0.05) {
    return("M. lepromatosis (UK Squirrel)")
  } else if (grepl("SRR367", node_name)) {
    return("M. leprae (UK Squirrel)")
  } else if (grepl("SRR847", node_name)) {
    return("M. leprae (Medieval UK Human)")
  } else {
    return("M. leprae (Modern Global Human)")
  }
}

# Apply the classification to all branches
classifications <- mapply(classify_tip, tree$tip.label, dist_to_root)

macro_meta <- data.frame(
  node = 1:N,
  label = tree$tip.label,
  Group = classifications
)

# Order the legend logically
macro_meta$Group <- factor(macro_meta$Group, levels=c(
  "M. lepromatosis (UK Squirrel)",
  "M. leprae (UK Squirrel)",
  "M. leprae (Medieval UK Human)",
  "M. leprae (Modern Global Human)"
))

# Define our high-contrast color palette
my_colors <- c(
  "M. lepromatosis (UK Squirrel)" = "#800080",  # Purple
  "M. leprae (UK Squirrel)"       = "#FF0000",  # Red
  "M. leprae (Medieval UK Human)" = "#0000FF",  # Blue
  "M. leprae (Modern Global Human)"= "#B0B0B0"  # Gray
)

# 3. Draw the Upgraded Macro Tree
pA <- ggtree(tree, linewidth=0.6) %<+% macro_meta +
  # Add colored points (alpha=0.8 makes them slightly see-through to reveal overlaps)
  geom_tippoint(aes(color=Group), size=3, alpha=0.8) +
  scale_color_manual(values=my_colors) +
  theme_tree2() +
  labs(title="Macroscopic Evolutionary Divergence",
       subtitle="Highlighting the dual-infection of UK squirrels with distinct leprosy bacilli",
       x="Genetic Distance (Substitutions per site)",
       color="Pathogen & Host Identity") +
  # Move the legend inside the plot area to make the figure wider and cleaner
  theme(legend.position=c(0.25, 0.85),
        legend.background=element_rect(fill="white", color="black", linewidth=0.5),
        legend.title=element_text(face="bold", size=12),
        legend.text=element_text(size=11),
        plot.title=element_text(face="bold", size=18),
        plot.subtitle=element_text(size=14, color="gray30"))

# 4. Save the beautiful new Panel A
ggsave(out_pdf_A, pA, width=12, height=9)

print(paste("✅ Upgraded Panel A saved to:", out_pdf_A))
plot(pA)
