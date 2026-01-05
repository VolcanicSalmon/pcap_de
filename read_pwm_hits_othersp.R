library(purrr)
library(dplyr)
library(stringr)
library(readr)
#before running this, read othersp_dna_binding_data.csv with tfogdat=read.csv()
genielinks_12=genielinks_us[genielinks_us$regulatoryGene%in%tfogdat$DataFrame,]
genielinks_12=genielinks_us[genielinks_us$weight>=0.05,]
genielinks_12_motifs=genielinks_12%>%left_join(reg_counts_motif,by=c('targetGene'='seq_id'))
genielinks_12_motifs
dir_path='/Volumes/Scratch/vef25hok/pcap_cluster_up_down/round2/othersp/Phytophthora_sojae/pcap_motifscans/'
files=list.files(dir_path,pattern='ft\\.txt$',full.names=T)
read_hitfile <- function(path) {
  read_tsv(
    path,
    comment = "#",
    col_names = c("seq_id", "start", "end", "strand", "score"),
    show_col_types = FALSE
  ) %>%
    mutate(
      file = basename(path),
      motif = str_remove(basename(path), "\\.ft\\.txt$")
    )
}
soja_hits <- map_df(files, read_hitfile)
soja_hits_dedu=soja_hits %>%
  select(-start, -end, -strand, -score, -file) %>%  
  distinct()  
soja_hits_homo=soja_hits_dedu[soja_hits_dedu$motif %in% genielinks_12_motifs$motif,]
soja_hits_homo_uniq=soja_hits_homo%>%mutate(seq_id = str_replace(seq_id, "-t26.*", ""))%>%distinct()
ramo_hits_homo=ramo_hits_dedu[ramo_hits_dedu$motif %in% genielinks_12_motifs$motif,]
ramo_hits_homo_uniq=ramo_hits_homo%>%mutate(seq_id = str_replace(seq_id, "-CDS.*", ""))%>%distinct()
infes_hits_homo=infes_hits_dedu[infes_hits_dedu$motif %in% genielinks_12_motifs$motif,]
infes_hits_homo_uniq=infes_hits_homo%>%mutate(seq_id = str_replace(seq_id, "-t26.*", ""))%>%distinct()
