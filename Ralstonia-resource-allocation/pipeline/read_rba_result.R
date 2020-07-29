# READ RBA RESULT TABLES
# ------------------------------------------------------------------------------
# generalized function to load and combine data from multiple result tables
library(stringi)
library(tidyverse)
read_rba_result <- function(files) {
  lapply(files, function(filename) {
    # read table
    df <- read_tsv(
      filename, 
      col_names = FALSE,
      col_types = cols()
      ) %>%
      rename(key = X1, value = X2)
    
    # make name for simulation column
    new_name = gsub("M_", "", filename)
    if (grepl("iteration", new_name)) {
      simulation = stri_extract_first(new_name, 
        regex = "iteration_[0-9]*")
      into = c("type", "iteration")
    } else {
      simulation = stri_extract_first(new_name,
        regex = "(succ|fru|for)_[0-9]*(\\.[0-9]*)?_nh4_[0-9]*(\\.[0-9]*)?_[0-9]*")
      into = c("carbon_source", "carbon_conc", "nitrogen_source", "nitrogen_conc", "sim_run")
    }
    
    # read data
    df %>% mutate(simulation = simulation) %>%
    separate(simulation, into = into, 
      sep = "_", remove = FALSE) %>%
    mutate_at(vars(matches("conc|iteration")), as.numeric)
  }) %>% dplyr::bind_rows()
}