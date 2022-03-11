library(dplyr)
library(stringr)
library(tidyr)
library(tibble)

#te annotated set of all identified insertions
complete_all_final <- read.csv(snakemake@input[["complete_all_final"]], header = TRUE)

#other three tables not yet subtype/element annotated
germinal_insertions <- read.csv(snakemake@input[["germinal_insertions"]], header = TRUE)
all_insertions_annotated <- read.csv(snakemake@input[["all_insertions_annotated"]], header = TRUE)
germinal_insertions_annotated <- read.csv(snakemake@input[["germinal_insertions_annotated"]], header = TRUE)


#read.csv("", header = TRUE)

#using all these columns of the un-/annotated germinal and annotated all tables might be overkill..
#but this way I can spot if sth. upstream destroyed the coordinates read counts

dplyr::left_join(abcd_final[200:400,1:6], abcd_final[,1:10], keep = FALSE, na_matches = "never", by=c("Chr"="Chr", 
                                                       "InsertionStart"="InsertionStart", 
                                                       "InsertionEnd"="InsertionEnd", 
                                                       "Sample"="Sample", 
                                                       "StartReads"="StartReads", 
                                                       "EndReads"="EndReads")) %>% 
  #rows not also found in our big type annotated table will hava NA under TotalReads ;D
  filter(!is.na(TotalReads)) %>% 
  nrow()
