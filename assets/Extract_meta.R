###
#
# Let's generate the medata file
#
##

# 0) Modules
# ------- ------- ------ ------ ----- #

library(tidyverse)

# 1) Path and Input
# ------- ------- ------ ------ ----- #

main.wd.p <- "/Users/stubbe/Desktop/Ressources/CUTandRUN/Run1/fastq/"
setwd(main.wd.p)

#Reads
read.f <- list.files(pattern = ".fastq.gz")

#Path of reads
path_for_read <- "/mnt/d/Academics/NGS/CR_reads/"
path_for_read <- main.wd.p
  
# 2) Construct the metadata file
# ------- ------- ------ ------ ----- #

#Get the librarires to work with
Libraries <- 
  read.f %>%
  str_split("_") %>%
  map_chr(~ str_c(.x[1:3], collapse = "_")) %>%
  unique()

#Construct the metadata
my_metadata <- map_dfr(Libraries, function(lib) {
  reads <- str_subset(read.f, lib)
  
  meta <- str_split(lib, "_")[[1]]
  
  tibble(
    id   = lib,
    group    = meta[2],
    replicate = meta[3],
    fastq_path_1   = str_c(path_for_read, reads[1], sep = ""),
    fastq_path_2   = str_c(path_for_read, reads[2], sep = "")
  )
})

my_metadata

# 3) Export the file
# ------- ------- ------ ------ ----- #

my_metadata %>% write_csv("/Users/stubbe/Desktop/CR_1_metadata_mac.csv")



