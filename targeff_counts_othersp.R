hya_targeff_counts <- hya_match1 %>%
  left_join(hyaeff, 
            by = c("FungiDB.68_HarabidopsidisEmoy2_AnnotatedProteins" = "gene")) %>%
  mutate(
    effector_type = case_when(
      is.na(Prediction) ~ "Not predicted",
      str_detect(Prediction, "Non-effector") ~ "Non-effector",
      str_detect(Prediction, "Cytoplasmic effector") ~ "Cytoplasmic effector",
      str_detect(Prediction, "Apoplastic") ~ "Apoplastic effector",
      str_detect(Prediction, "/apoplastic effector") ~ "Cytoplasmic/apoplastic effector",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Orthogroup, effector_type) %>%
  summarise(
    count = n(),
    proteins = paste(unique(FungiDB.68_HarabidopsidisEmoy2_AnnotatedProteins), collapse = ", "),
    .groups = 'drop'
  ) %>%
  arrange(Orthogroup, desc(count))
