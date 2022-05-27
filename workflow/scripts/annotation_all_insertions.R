library(dplyr)
library(fuzzyjoin)
library(readxl)
library(IRanges)


snake_all_ins <- snakemake@input[["all"]]
snake_annotation <- snakemake@input[["annotation"]]
snake_extension <- snakemake@params[["extension"]]


## Read in tables with all identified insertions; 
all_ins <- read.csv(snake_all_ins, header=TRUE)
all_ins <- all_ins %>%
  as_tibble() %>%
  mutate(Chr = as.character(Chr))  

## Read Annotation File
annotation <- read.delim(snake_annotation, header=FALSE, comment.char="#")

# add column names
colnames(annotation) <- c("GeneID", "Chr", "GeneStart", "GeneEnd")

# reduce genes to ONE entry which is the longest
# multiple entries per gene likely, since after gffread, transcripts of different lengths are included (share name)
# we find the min and max positions (start & end coordinate with greatest distance) per gene and keep this record  
annotation <- annotation %>%
  group_by(GeneID, Chr) %>%
  summarise(GeneStart = min(GeneStart), GeneEnd = max(GeneEnd))

# compute gene_length
annotation <- annotation %>% 
  mutate(Gene_length = GeneEnd - GeneStart + 1)


# add bp extension as defined in the config.yaml
# this is useful for inclusion of upstream/downstream regions or if UTRs are poorly characterized
annotation <- annotation %>%
  mutate(GeneStart = GeneStart - snake_extension, GeneEnd = GeneEnd + snake_extension)


#head(as_tibble(all_ins))
#head(annotation)

# join MuGerminal table with annotation
all_ins_annotated <- fuzzyjoin::genome_inner_join(all_ins, 
                                annotation, 
                                by=c("Chr", "InsertionStart"="GeneStart", "InsertionEnd"="GeneEnd") 
                                )

print("all good!")

# rename, relocate and drop the columns to create a useful table
all_ins_annotated <- all_ins_annotated %>%
  dplyr::rename(
    Chr = "Chr.y",
#    GeneStart = "GeneStart.y",
#    GeneEnd = "GeneEnd.y",
#    InsertionStart = "GeneStart.x",
#    InsertionEnd = "GeneEnd.x"
  ) %>%
  relocate(
    Sample, InsertionStart, InsertionEnd, StartReads, EndReads, .before = Gene_length
  ) %>%
  select(
    -Chr.x
  )



# order the table as to have intersecting samples next to one another
# + remove extension so as to have correct gene coordinates for start & end
all_ins_annotated <- all_ins_annotated %>%
  arrange(., Chr, GeneStart, InsertionStart) %>%
  mutate(GeneStart = GeneStart + snake_extension, GeneEnd = GeneEnd - snake_extension)


write.csv(all_ins_annotated, 
          snakemake@output[[1]], 
          row.names=F,
          quote=F)
