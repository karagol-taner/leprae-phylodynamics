# 1. Load Required Packages (Assuming they are already installed from the last run)
library(ggtree)
library(treeio)
library(ggplot2)
library(ape)

print("--- Constructing cleanly aligned annotated phylogeny ---")

# Define paths
tree_file <- "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/lep_phylogeny_final.treefile"
meta_file <- "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/leprosy_metadata.csv"
out_pdf <- "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/ggtree_annotated_publication_aligned.pdf"

# 2. Load the Tree and Clean Tip Labels
tree <- read.tree(tree_file)
tree$tip.label <- gsub(".final.bam", "", tree$tip.label)
tree$tip.label <- sapply(strsplit(tree$tip.label, "/"), tail, 1)

# 3. The Digital Scalpel: Remove M. lepromatosis
N <- length(tree$tip.label)
dist_to_root <- dist.nodes(tree)[N + 1, 1:N]
outgroups <- tree$tip.label[dist_to_root > 0.05]
tree_pruned <- drop.tip(tree, outgroups)

# 4. Load the Metadata Dictionary
metadata <- read.csv(meta_file)

# 5. Draw and Annotate the Tree with ALIGNED TEXT
p <- ggtree(tree_pruned, linewidth=0.5) %<+% metadata +

  # Add colored points at the end of the branches based on the Host
  geom_tippoint(aes(color=Host), size=2) +

  # THE FIX: Align the text cleanly to the right, shrink font, add dotted lines
  geom_tiplab(aes(label=Location, color=Host),
              size=2.5,          # Smaller font to prevent vertical collision
              align=TRUE,        # Aligns all text into a neat column on the far right
              linesize=0.3,      # Adds a thin dotted line connecting the branch to the text
              offset=0.0002) +   # Adds a tiny gap between the point and the dotted line

  # Apply our specific color palette
  scale_color_manual(values=c("Red Squirrel"="#FF0000", "Human"="#0000FF", "Global"="#808080")) +

  # Add the x-axis scale and clean up the theme
  theme_tree2() +
  labs(title="Evolutionary Phylogeny of M. leprae",
       subtitle="Evidence of Zoonotic Transmission Dynamics",
       x="Genetic Distance (Substitutions per site)",
       color="Host Species") +

  # Expand x-axis further to fit the newly aligned right-hand text column
  xlim(0, 0.015)

# 6. Save the High-Resolution Publication PDF (Increased height to 30 inches!)
ggsave(out_pdf, p, width=12, height=30, limitsize=FALSE)

print(paste("✅ Perfectly aligned vector figure saved to:", out_pdf))

# Preview the plot in Colab (Note: the inline preview might look slightly compressed, but the saved PDF will be perfect)
plot(p)
