library(tidyverse)
process_ortholog_file <- function(csv_path, 
                                  geneprot, 
                                  genielinks_12_motifs, 
                                  reg_counts_motif_rename,
                                  suffix_to_remove = '_RNA-p1' #or -t26.. for other species
                                 ) {
 
  df <- read.csv(csv_path, sep = '\t', header = TRUE)

  pcap_col <- names(df)[grep('pcap', names(df), ignore.case = TRUE)][1]
  protein_col <- names(df)[grep('FungiDB|Annotated', names(df))][1]
  
  message(paste("Using pcap column:", pcap_col))
  message(paste("Using protein column:", protein_col))
  
  
  df <- df %>%
    separate_rows(!!sym(protein_col), sep = ", ")  
  
  df[[protein_col]] <- sub(suffix_to_remove, '', df[[protein_col]])
  
 
  df <- df %>%
    left_join(geneprot, by = setNames('V1', pcap_col))
  
 
  df <- df %>%
    mutate(V2 = sub('-', '_', V2))
  df_motifs <- df %>%
    left_join(genielinks_12_motifs, by = c('V2' = 'regulatoryGene'))
  
  
  match1 <- df_motifs %>%
    inner_join(reg_counts_motif_rename, 
               by = c('targetGene' = 'pcap_id', 'motif' = 'motif'))
  
 
  return(match1)
}
para_match1 <- process_ortholog_file(
  csv_path = '/Volumes/Scratch/vef25hok/pcap_cluster_up_down/round2/othersp/Orthrun/OrthoFinder/Results_Dec21/Orthologues/Orthologues_FungiDB-68_PinfestansT30-4_AnnotatedProteins/FungiDB-68_PinfestansT30-4_AnnotatedProteins__v__pcap_top12tf.tsv',
  geneprot = geneprot,
  genielinks_12_motifs = genielinks_12_motifs,
  reg_counts_motif_rename = reg_counts_motif_rename
)
ramo_match1=process_ortholog_file(
  csv_path = '/Volumes/Scratch/vef25hok/pcap_cluster_up_down/round2/othersp/Orthrun/OrthoFinder/Results_Dec21/Orthologues/Orthologues_FungiDB-68_Pramorum14567_AnnotatedProteins/FungiDB-68_Pramorum14567_AnnotatedProteins__v__pcap_top12tf.tsv',
  geneprot = geneprot,
  genielinks_12_motifs = genielinks_12_motifs,
  reg_counts_motif_rename = reg_counts_motif_rename
)
