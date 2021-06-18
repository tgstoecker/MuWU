library(dplyr)
library(fuzzyjoin)
library(readxl)
library(IRanges)


snake_germinal_ins <- snakemake@input[["germinal"]]
snake_grid <- snakemake@input[["grid_table"]]
snake_annotation <- snakemake@input[["annotation"]]
snake_extension <- snakemake@params[["extension"]]

## Read in tables with all identified insertions; 
germinal_ins <- read.csv(snake_germinal_ins, header=TRUE)

## Read in stock matrix
MY_STOCK <- list.files("config/stock_matrix/", pattern="\\.xlsx$")
stock <- as.data.frame(read_excel(paste0("config/stock_matrix/", MY_STOCK)))
names(stock)[1] <- "X__1"
rownames(stock) <- stock$X__1
stock <- stock[, -1]
stock <- as.matrix(stock)


#grid table - to associate arbitrary names with row/col
grid_table <- read.csv(snake_grid, sep = "\t", header=TRUE)
grid_table <- grid_table %>%
                select(dim, base_name)

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

# add bp extension as defined in the config.yaml
# this is useful for inclusion of upstream/downstream regions or if UTRs are poorly characterized
annotation <- annotation %>%
  mutate(Start = Start - snake_extension, End = End + snake_extension)

# join MuGerminal table with annotation
germinal_ins_annotated <- fuzzyjoin::genome_inner_join(germinal_ins,
                                annotation,
                                by=c("Chr", "Start", "End")
                                )



# rename, relocate and drop the columns to create a useful table
germinal_ins_annotated <- germinal_ins_annotated %>%
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


#merge annotated germinal ins. table with grid table
germinal_ins_grid_info <- merge(germinal_ins_annotated, grid_table, by.x = "Sample", by.y = "base_name")

# order the table as to have intersecting samples next to one another
germinal_ins_grid_info <- germinal_ins_grid_info %>%
  arrange(., Chr, Start, InsertionStart)


## assigning Stock
#create seperate dataframes - just with Cols and just with Rows
just_row <- germinal_ins_grid_info %>%
  filter(dim == "row") %>%
  select(-dim)
  
just_col <- germinal_ins_grid_info %>%
  filter(dim == "col") %>%
  select(-dim)

full <- c()
for (i in 1:nrow(just_row))  {
  full[i] <- stock[cbind(just_row$Sample[i], just_col$Sample[i])]
}

just_col$stock <- full
just_row$stock <- full

germinal_ins_annotated_stocks <- rbind(just_row, just_col)

# order the table as to have intersecting samples next to one another; again..
germinal_ins_annotated_stocks <- germinal_ins_annotated_stocks %>%
  arrange(., Chr, Start, InsertionStart) %>%
  relocate(Sample, .before = InsertionStart)

write.csv(germinal_ins_annotated_stocks,
          snakemake@output[[1]],
          row.names=F)
