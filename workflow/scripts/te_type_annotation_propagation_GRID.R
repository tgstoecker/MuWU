library(dplyr)
library(stringr)
library(tidyr)
library(tibble)

#te annotated set of all identified insertions
complete_all_final <- read.csv(snakemake@input[["complete_all_final"]], header = TRUE) %>% as_tibble() %>% mutate(Chr = as.character(Chr))

#other three tables not yet subtype/element annotated
germinal_insertions <- read.csv(snakemake@input[["germinal_identified_insertions"]], header = TRUE) %>% as_tibble() %>% mutate(Chr = as.character(Chr))
all_insertions_annotated <- read.csv(snakemake@input[["all_identified_insertions_annotated"]], header = TRUE) %>% as_tibble() %>% mutate(Chr = as.character(Chr))
germinal_insertions_annotated <- read.csv(snakemake@input[["germinal_identified_insertions_annotated"]], header = TRUE) %>% as_tibble() %>% mutate(Chr = as.character(Chr))


#using all these columns of the un-/annotated germinal and annotated all tables might be overkill..
#but this way I can spot if sth. upstream destroyed the coordinates read counts


### 1 - Germinal insertions

germinal_ait_annotated_final <- dplyr::left_join(germinal_insertions, 
                        complete_all_final, 
                        keep = FALSE, 
                        na_matches = "never", 
                        by=c("Chr"="Chr", 
                             "InsertionStart"="InsertionStart", 
                             "InsertionEnd"="InsertionEnd", 
                             "Sample"="Sample", 
                             "StartReads"="StartReads", 
                             "EndReads"="EndReads")) %>% 
  #rows not also found in our big type annotated table will hava NA under TotalReads ;D
  filter(!is.na(TotalReads)) 


write.csv(germinal_ait_annotated_final, "test.output")
