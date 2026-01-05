library(ggtree)
library(ape)
library(stringr)
library(dplyr)
library(tidytree)
library(readr)

othersp_color_group <- function(tip_labels) {
  sapply(tip_labels, function(x) {
    parts <- str_split(x, "_")[[1]]
    return(parts[1])  
  })
}

plot_colored_tree <- function(treepath, root_genes = NULL) {
  tree <- read.tree(treepath)
  tip_labels <- tree$tip.label
  color_groups <- othersp_color_group(tip_labels)
  
  if (!is.null(root_genes) && length(root_genes) > 0) {
    matching_tips <- tip_labels[sapply(tip_labels, function(tip) {
      any(sapply(root_genes, function(gene) grepl(gene, tip, fixed = TRUE)))
    })]
    
    if (length(matching_tips) > 0) {
      if (length(matching_tips) > 1) {
        message("  Found ", length(matching_tips), " matching tips, using first one: ", matching_tips[1])
        matching_tips <- matching_tips[1]
      } else {
        message("  Rooting on: ", matching_tips)
      }
      tree <- root(tree, outgroup = matching_tips, resolve.root = TRUE)
    } else {
      message("  No matching tips found for rooting")
    }
  }
  
  tip_data <- data.frame(label = tip_labels, group = color_groups, stringsAsFactors = FALSE)
  
  # Calculate maximum label length to adjust offset
  max_label_length <- max(nchar(tip_labels))
  # Adjust offset based on label length (less offset = more space)
  label_offset <- 0.001  # Very small offset
  
  # Plot with legend on LEFT and more space for labels
  p <- ggtree(tree, layout = "rectangular") %<+% tip_data +
    geom_tiplab(aes(color = group), 
                offset = label_offset, 
                size = 6,
                hjust = 0,  # Left-align labels
                align = FALSE) +  # Don't align tips
    scale_color_discrete(name = "Species") +
    theme(
      legend.position = "left",  # LEGEND ON LEFT
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      plot.margin = margin(10, 150, 10, 10)  # Extra margin on RIGHT for labels
    ) +
    xlim(NA, max(node.depth.edgelength(tree)) * 1.5)  # Extend x-axis for labels
  
  return(p)
}

mapping_file <- "/Volumes/Scratch/vef25hok/pcap_cluster_up_down/round2/othersp/Orthrun/OrthoFinder/Results_Dec21/Orthologues/Orthologues_pcap_top12tf/pcap_top12tf__v__FungiDB-68_PcapsiciLT1534_AnnotatedProteins.tsv"
ortho_mapping <- read_tsv(mapping_file, col_types = cols())

root_gene_list <- ortho_mapping %>%
  group_by(Orthogroup) %>%
  summarise(root_genes = list(unique(pcap_top12tf)[1]), .groups = 'drop') %>%
  deframe()

treedir_path <- "/Volumes/Scratch/vef25hok/pcap_cluster_up_down/round2/othersp/Orthrun/OrthoFinder/Results_Dec21/"
othertree_files <- list.files(
  path = treedir_path, 
  pattern = "_tree_dedup\\.txt$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(othertree_files) == 0) {
  stop("No tree files found.")
}

output_dir <- "~/Documents/scwgcna/dev28/"
output_dir <- path.expand(output_dir)

suppressWarnings({
  for (file in othertree_files) {
    og_id <- str_extract(basename(file), "OG[0-9]+")
    message("\nPlotting: ", basename(file), " (", og_id, ")")
    
    root_genes <- root_gene_list[[og_id]]
    
    if (!is.null(root_genes)) {
      message("  Root gene: ", root_genes)
    }
    
    tryCatch({
      tree <- read.tree(file)
      n_tips <- length(tree$tip.label)
      
     
      pdf_height <- max(10, min(40, n_tips * 0.4))
      pdf_width <- 24  
      otherp <- plot_colored_tree(file, root_genes = root_genes)
      
      pdf_filename <- paste0(og_id, "_tree.pdf")
      pdf_path <- file.path(output_dir, pdf_filename)
      
      ggsave(
        filename = pdf_path,
        plot = otherp,
        width = pdf_width,
        height = pdf_height,
        units = "in",
        limitsize = FALSE  # Allow large sizes
      )
      
      print(otherp)
      
    }, error = function(e) {
      message("check ", basename(file), ": ", e$message)
    })
  }
})
