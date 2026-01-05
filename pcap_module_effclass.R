#before running this, read p.capsici effectorP outputs and name the variable pcapeff
process_module_effectors <- function(file_path, effector_df) {
 
  links <- links %>%filter(weight>=0.05)%>%
    mutate(
      regulatoryGene = str_replace_all(regulatoryGene, "-", "_"),
      targetGene = str_replace_all(targetGene, "-", "_")
    )
  

  links_with_effector <- links %>%
    left_join(effector_df, by = c("targetGene" = "gene")) %>%
    mutate(
      effector_type = case_when(
        is.na(Prediction) ~ "Not predicted",
        str_detect(Prediction, "Non-effector") ~ "Non-effector",
        str_detect(Prediction, "Cytoplasmic effector") ~ "Cytoplasmic effector",
        str_detect(Prediction, "Apoplastic") ~ "Apoplastic effector",
        str_detect(Prediction, "/apoplastic effector") ~ "Cytoplasmic/apoplastic effector",
        TRUE ~ "Other"
      )
    )
  
  
  effector_counts <- links_with_effector %>%
    count(effector_type) %>%
    arrange(desc(n))
  
  
  module <- str_extract(basename(file_path), "m[0-9]+")
  

  effector_counts <- effector_counts %>%
    mutate(module = module) %>%
    select(module, effector_type, n)
  
  return(list(
    counts = effector_counts,
    full_data = links_with_effector
  ))
}


base_path <- "~/Downloads/genie3_module_08121201"
modules <- paste0("m", 1:9)

# Store results
all_counts <- list()
all_data <- list()

for (module in modules) {
  file_path <- file.path(base_path, paste0("genie3_links_Haustorium-", module, ".csv"))
  
  if (file.exists(file_path)) {
   
    
    result <- process_module_effectors(file_path, pcapeff)  
    
    all_counts[[module]] <- result$counts
    all_data[[module]] <- result$full_data
    
  } else {
    cat("File not found:", file_path)
  }
}


combined_counts <- bind_rows(all_counts)


print(combined_counts)


counts_wide_thres005 <- combined_counts %>%
  pivot_wider(
    names_from = effector_type,
    values_from = n,
    values_fill = 0
  )

print(counts_wide_thres005)


#Visualisation
ggplot(counts_long_thres005, aes(x = module, y = count, fill = effector_type)) +
  geom_col() +
  labs(title = "Binding Target Effector Classification",
       x = "Module",
       y = "Count",
       fill = "Effector Classes") +
  theme_minimal() +
  theme(legend.position = "right")
