library(dplyr)
library(fuzzyjoin)
library(readxl)
library(IRanges)


snake_all_ins <- snakemake@input[["all"]]
snake_annotation <- snakemake@input[["annotation"]]


## Read in tables with all identified insertions; 
all_ins <- read.csv(snake_all_ins, header=TRUE)

## Read Annotation File
annotation <- read.delim(snake_annotation, header=FALSE, comment.char="#")

# add column names
colnames(annotation) <- c("GeneID", "Chr", "Start", "End")

# reduce genes to ONE entry which is the longest
# multiple entries per gene likely, since after gffread, transcripts of different lengths are included (share name)
# we find the min and max positions (start & end coordinate with greatest distance) per gene and keep this record  
annotation <- annotation %>%
  group_by(GeneID, Chr) %>%
  summarise(Start = min(Start), End = max(End))

# compute gene_length
annotation <- annotation %>% 
  mutate(Gene_length = End - Start + 1)

# join MuGerminal table with annotation
all_ins_annotated <- fuzzyjoin::genome_inner_join(all_ins, 
                                annotation, 
                                by=c("Chr", "Start", "End") 
                                )



# rename, relocate and drop the columns to create a useful table
all_ins_annotated <- all_ins_annotated %>%
  dplyr::rename(
    Chr = "Chr.y",
    Start = "Start.y",
    End = "End.y",
    InsertionStart = "Start.x",
    InsertionEnd = "End.x"
  ) %>%
  relocate(
    Sample, InsertionStart, InsertionEnd, StartReads, EndReads, .before = Gene_length
  ) %>%
  select(
    -Chr.x
  )



# order the table as to have intersecting samples next to one another
all_ins_annotated <- all_ins_annotated %>%
  arrange(., Chr, Start, InsertionStart)


write.csv(all_ins_annotated, 
          "results/insertions_table_final/all_identified_insertions_annotated.csv", 
          row.names=F)
