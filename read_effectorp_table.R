hyaeff=read.csv('/Volumes/Scratch/vef25hok/pcap_cluster_up_down/round2/othersp/Hyaloperonospora_arabidopsidis/FungiDB-68_HarabidopsidisEmoy2_AnnotatedProteins.efp.txt',sep='\t')
read_effector_file <- function(file_path, 
                               identifier_col = "X..Identifier",
                               extract_gene = TRUE,
                               remove_identifier = TRUE) {
  
  library(tidyverse)
  
  # Read file
  df <- read.csv(file_path, sep = '\t')
  
  # Extract gene if requested
  if (extract_gene) {
    df <- df %>%
      mutate(gene = str_extract(.data[[identifier_col]], "gene=([^|]+)") %>%
               str_remove("gene=") %>%
               str_trim())
  }
  
  # Remove identifier column if requested
  if (remove_identifier) {
    df <- df %>%
      select(-all_of(identifier_col))
  }
  
  return(df)
}


hyaeff <- read_effector_file('/Volumes/Scratch/vef25hok/pcap_cluster_up_down/round2/othersp/Hyaloperonospora_arabidopsidis/FungiDB-68_HarabidopsidisEmoy2_AnnotatedProteins.efp.txt')
sojaeff=read_effector_file('/Volumes/Scratch/vef25hok/pcap_cluster_up_down/round2/othersp/Phytophthora_sojae/FungiDB-68_PsojaeP6497_AnnotatedProteins.efp.txt')
paraeff=read_effector_file('/Volumes/Scratch/vef25hok/pcap_cluster_up_down/round2/othersp/Phytophthora_parasitica/FungiDB-68_PparasiticaINRA-310_AnnotatedProteins.efp.txt')
ramoeff=read.csv('/Volumes/Scratch/vef25hok/pcap_cluster_up_down/round2/othersp/Orthrun/FungiDB-68_Pramorum14567_AnnotatedProteins.efp.txt',sep='\t',header=T)%>%
  mutate(gene = str_extract(X..Identifier, "KAH[0-9]+\\.[0-9]+"))



