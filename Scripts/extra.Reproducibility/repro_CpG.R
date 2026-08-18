library(data.table)
library(dplyr)

mod <- "5mC"
dir <- file.path("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_proportions/", mod)

files <- list.files(
  path = "path/to/your/files",
  pattern = "\\.bed$",
  full.names = TRUE
)

files


# At some point, check how many sites were dropped because of coverage