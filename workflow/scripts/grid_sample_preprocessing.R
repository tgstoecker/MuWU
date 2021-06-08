library(readxl)
library(stringr)
library(dplyr)

# read in stock matrix and transform to usable matrix dataset
MY_STOCK <- list.files("config/stock_matrix/", pattern="\\.xlsx$")
stock <- as.data.frame(read_excel(paste0("config/stock_matrix/", MY_STOCK)))
names(stock)[1] <- "X__1"
rownames(stock) <- stock$X__1
stock <- stock[, -1]
stock <- as.matrix(stock)

# compute grid dimension length of row and column
row_dim_length <- length(rownames(stock))
col_dim_length <- length(colnames(stock))

# function returning all files sharing pattern
list_all_files <- function(x){
  list.files(path="rawreads/", pattern=x)  
}

# apply "list_all_files" function to all row & column samples
all_row_samples <- sapply(X = base::rownames(stock), FUN = list_all_files)
all_col_samples <- sapply(X = base::colnames(stock), FUN = list_all_files)


# decider condition in which mode to continue
# first we check if row/column sample files are themselves consistently matrixes -> thus PE
# in else if, we check if we truly have one sample per row or column position -> thereby SE
# for both SE or PE we write a read_type.yaml containing info of SE/PE - that we can parse later
if (is.matrix(all_row_samples) & is.matrix(all_col_samples)) {
    print("Paired-end data detected - continuing in PE mode")
    PE <- TRUE
    read_type_file <- file("config/read_type.yaml")
    writeLines(c("read_type: PE"), read_type_file)
    close(read_type_file)
} else if (length(all_row_samples) == row_dim_length & length(all_col_samples) == col_dim_length) {
    print("Single-end data detected - continuing in SE")
    PE <- FALSE
    read_type_file <- file("config/read_type.yaml")
    writeLines(c("read_type: SE"), read_type_file)
    close(read_type_file)
} else {
    print("Input data seems incorrect - inconsistent # of fastq files for either SE or PE input set")
}


## create row and column tibbles
## also add info of row/column ("dim") + and sample position on the grid ("num")

# whether PE or SE (so PE == FALSE) - create appropiate sample and mode sheet
if (PE == TRUE) {
  row_tibble <- tibble(
    dim="row",
    num=seq.int(ncol(all_row_samples)),
    base_name = colnames(all_row_samples),
    fq_1_file = all_row_samples[1,],
    fq_2_file = all_row_samples[2,],
    row.names = NULL  
  )

  col_tibble <- tibble(
    dim="col",
    num=seq.int(ncol(all_col_samples)),
    base_name = colnames(all_col_samples),
    fq_1_file = all_col_samples[1,],
    fq_2_file = all_col_samples[2,],
    row.names = NULL  
  )

  row_col_tibble <- rbind(row_tibble, col_tibble)
  row_col_tibble

  row_col_tibble_final <- row_col_tibble %>%
    mutate(fq_1_end = str_remove(fq_1_file, pattern = c(base::rownames(stock), base::colnames(stock)))) %>%
    mutate(fq_2_end = str_remove(fq_2_file, pattern = c(base::rownames(stock), base::colnames(stock))))
    
  # write generated sample sheet to 'config/grid_sample_sheet.tsv'
  write.table(row_col_tibble_final, 
              file='config/grid_sample_sheet.tsv', 
              quote=FALSE, 
              sep='\t', 
              col.names = TRUE,
              row.names = FALSE)
    
} else if (PE == FALSE) {
    row_tibble <- tibble(
      dim = "row",
      num = seq.int(length(all_row_samples)),
      base_name = names(all_row_samples),
      fq_1_file = all_row_samples
    )

  col_tibble <- tibble(
    dim="col",
    num=seq.int(length(all_col_samples)),
    base_name = names(all_col_samples),
    fq_1_file = all_col_samples
  )

  row_col_tibble <- rbind(row_tibble, col_tibble)

  row_col_tibble_final <- row_col_tibble %>%
    mutate(fq_1_end = str_remove(fq_1_file, pattern = c(base::rownames(stock), base::colnames(stock))))  
    
  # write generated sample sheet to 'config/grid_sample_sheet.tsv'
  write.table(row_col_tibble_final, 
              file='config/grid_sample_sheet.tsv', 
              quote=FALSE, 
              sep='\t', 
              col.names = TRUE,
              row.names = FALSE)
    
} else {
    print("Sth. somewhere went horribly wrong. You should not be able to be here.")
    print("Edge case encountered - contact maintainer of MuWU.")
}
