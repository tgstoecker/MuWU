library(dplyr)
library(fuzzyjoin)
library(readxl)
library(IRanges)
library(stringr)
library(tidyr)
library(tibble)
library(data.table)
library(foreach)
library(doParallel)

#get snakemake variables
message("Getting snakemake variables:")

samples <- snakemake@params[["samples"]]
#samples <- c('Row_01', 'Col_01',
#             'Row_02', 'Col_02',
#             'Row_03', 'Col_03',
#             'Row_04', 'Col_04',
#             'Row_05', 'Col_05',
#             'Row_06', 'Col_06',
#             'Row_07', 'Col_07')
all_types <- snakemake@params[["all_types"]]
te_typing_cluster_cores <- snakemake@params[["te_typing_cluster_cores"]]
insertion_table_file <- snakemake@input[["all_identified_insertions"]]
#insertion_table_name <- snakemake@params[["insertion_table_name"]]

#germinal_annotated <- snakemake@input[["germinal_annotated"]]
#germinal <- snakemake@input[["germinal"]]
#all_annotated <- snakemake@input[["all_annotated"]]
#all_NOT_annotated <- snakemake@input[["all"]]


#input_list <- list(germinal_annotated, 
#                   germinal_NOT_annotated, 
#                   all_annotated, 
#                   all_NOT_annotated)

#names(input_list) <- c("germinal_annotated", 
#                       "germinal_NOT_annotated", 
#                       "all_annotated", 
#                       "all_NOT_annotated")





#for (i in 1:length(input_list)) {
#  print(names(input_list[i]))
#  print(head(as.character(input_list[i])))
#}

base_all_types <- all_types %>% 
  str_remove(., c("_L")) %>%
  str_remove(., c("_R")) %>%
  unique()

message("All supplied TE types:")
base_all_types


#OMG - this was really annoying xD
#depending on where one uses data.table one has to consider the local amount of cores..
#https://stackoverflow.com/questions/23107237/r-3-1-0-crashes-with-segfault-when-loading-data-table-package-1-9-2
#apparently the default number of threads is set by data.table to the max number available
#on a system like our work compute server with 128 cores and some packages you are bound to run into memory and address problems

#one can either set global params in bash etc. e.g.: 
#setenv OMP_NUM_THREADS 1

#the other option is to use the openmp utils - https://jangorecki.gitlab.io/data.cube/library/data.table/html/openmp-utils.html,
#to quickly set the number of threads to sth. feasible
#data.table::getDTthreads()
data.table::setDTthreads(1)
#data.table::getDTthreads()


### Create cluster function for on-demand parallelization ###
#FORK cluster since I expect a Linux machine
#autostop=TRUE since I don't want to handle this manually
#with my.cluster & stopCluster(my.cluster) I could check the status
message("Create cluster function for on-demand parallelization:")

setup_cluster <- function(){

  #define cluster
  parallel::detectCores()
  n.cores <- te_typing_cluster_cores
  n.cores

  #create the cluster - FORK because this way libraries, variables etc. are copied to the clusters!
  my.cluster <- parallel::makeForkCluster(
    n.cores, 
    type = "FORK",
    autostop=TRUE
  )

  #check cluster definition (optional)
  print(my.cluster)

  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)

  #check if it is registered (optional)
  print(
    foreach::getDoParRegistered()
  )
  #how many workers are available? (optional)
  print(
    foreach::getDoParWorkers()
  )

}


### Create a function to completely get rid of a used FORK cluster & free up memory
message("Create a function to completely get rid of a used FORK cluster & free up memory")

#a function to get and close connection (to cluster sockets in our case)
burn_socks <- function(x){
  close.connection(getConnection(x))
}


#function to truly get rid of old Cluster/sockets
rm_cluster <- function(){
  stopImplicitCluster()

  connections <- showConnections(all = FALSE) 
  socket_connections <- as.data.frame(connections) %>%
    filter(class == "sockconn") %>%
    rownames()

  message("Removing all unwanted FORK connections - purging closed cluster sockets")
  message("This will kill zombie proesses and free up RAM")

  lapply(X = socket_connections, 
         FUN = burn_socks)
}



### Get sample level files with te type information ###
message("Get sample level files with te type information:")

#setup & start FORK cluster
setup_cluster()

#read-in final sample files
type_files_list <- foreach(s = samples, .final = function(s) setNames(s, paste0(samples, "_type_file"))) %dopar% {
  tmp_type_file <- fread(paste0("results/te_typing/pre_sorting/", s, "/", s, "_te_types_merged.tsv"), header = TRUE, sep="\t")
  tmp_type_file$Strand <- as.character(tmp_type_file$Strand)
  assign(paste0(s, "_type_file"), tmp_type_file)
    
  return(get(paste0(s, "_type_file")))
}

##really cool - extracting dfs/tibbles from list object by simply sending them to the global environment
#https://stackoverflow.com/questions/59169631/split-a-list-into-separate-data-frame-in-r
list2env(type_files_list, envir = .GlobalEnv)

#get rid of cluster
rm_cluster()


### Get SAM files and perform some formatting ###
message("Get SAM files and perform some formatting")

#setup & start FORK cluster
setup_cluster()

system.time(

#need to load all sam files
sam_object_list <- foreach(s = samples, .final = function(s) setNames(s, paste0("sam_", samples))) %dopar% {
    SAM_path <- paste0("results/dedup_sam/", s, ".dedup.sam")
    SAM_object <- paste0("sam_", s)
#    assign(SAM_object, read.csv(SAM_path, sep = "\t", header = FALSE, row.names = NULL))
    assign(SAM_object, fread(SAM_path, select=c(1,2,3,4,10), sep = "\t", header=FALSE, fill=TRUE))
    #take only the first four columns and the sequence
#    assign(SAM_object, get(SAM_object)[,c(1:4,10)])
    assign(SAM_object,
           get(SAM_object) %>%
             #some renaming of the header
             dplyr::rename(Name=V1, Flag=V2, Chr=V3, Start=V4) %>%
             #set type of chromosome column to character
             mutate(Chr = as.character(Chr)) %>%
             #compute lengths of the alignments
             mutate(Length = str_length(V10)) %>%
             #compute and add End coordinate of reads
             mutate(End = Start+Length-1) %>%
             #drop Lengh & V10 (Sequence) columns
             select(-V10, -Length)
    )
    return(get(SAM_object))
}

)
                           
##really cool - extracting dfs/tibbles from list object by simply sending them to the global environment
#https://stackoverflow.com/questions/59169631/split-a-list-into-separate-data-frame-in-r
list2env(sam_object_list, envir = .GlobalEnv)

#get rid of cluster
rm_cluster()



### Function to read SAM bitwise encoding - and return whether read is forward or reverse mate ###

#quick function to analyse bitwise encoding of SAM flag
#modified, but core idea idea from here - https://stackoverflow.com/a/12088263
number2binary = function(number, noBits) {
       binary_vector = rev(as.numeric(intToBits(number)))
       if(missing(noBits)) {
          return(rev(binary_vector))
       } else {
          return(rev(binary_vector[-(1:(length(binary_vector) - noBits))]))
       }
}



### get insertion files ###
## looping over all 4 files (germinal/all un-/annotated)  ##

message("Getting insertion files:")


#for (i in 1:length(input_list)) {

  #get name of current input file (input files were named at the start of the script)
#  name_input <- names(input_list[i])

insertions <- read.csv(as.character(insertion_table_file), header=TRUE)
insertions <- insertions %>%
  filter(Sample %in% samples) %>%
  #set type of chromosome column to character
  mutate(Chr = as.character(Chr))

#create insertions table with additional columns per (Mu) element/type
insertions_typed <- insertions

insertions_typed[,all_types]=0
#add additional column for uncategorized reads
insertions_typed <- insertions_typed %>%
  add_column(uncategorized = 0)


### Annotating with type information  ###
message("Annotating with type information:")

setup_cluster()

pre_ait_annotated <- foreach(row = 1:nrow(insertions_typed)) %dopar% {
  
    #take row to work on
    tmp_ins <- insertions_typed[row,]

    #match Sample with sam object
    tmp_sam <- get(paste0("sam_", tmp_ins$Sample))

    #match sample of insertion file with te typing file
    tmp_type_file <- get(paste0(tmp_ins$Sample, "_type_file"))


    # fuzzyjoin for overlap + end or start needs to match
    tmp_merge_ins_sam <- fuzzyjoin::genome_inner_join(tmp_ins,
                             tmp_sam, 
                              by=c("Chr", "InsertionStart"="Start", "InsertionEnd"="End")
                            ) %>%
                       select(Name, Flag, Chr.y, Start, End) %>%
                       dplyr::rename(Chr=Chr.y) %>%
                       filter(Start == tmp_ins$InsertionStart | 
                              #Start == tmp_ins$InsertionEnd | 
                              #End == tmp_ins$InsertionStart | 
                              End == tmp_ins$InsertionEnd
                             ) %>%
                       #translate Flag into forward (+) or reverse (-) strand mapping/read
                       rowwise() %>%
                       mutate(
                         Strand = case_when(
                           number2binary(Flag, 8)[7] == 1 & number2binary(Flag, 8)[8] == 0 ~ "1",
                           number2binary(Flag, 8)[7] == 0 & number2binary(Flag, 8)[8] == 1 ~ "2",
                           TRUE ~ "LOST"
                         )
                       ) %>%
                       select(-Flag) %>%
                       relocate(Strand, .after=Name)

    #associate sam reads with type information
    tmp_merge_ins_sam_typed_summary <- dplyr::left_join(tmp_merge_ins_sam, 
                                                      tmp_type_file, 
                                                      by = c("Name"="Name", "Strand"="Strand"), 
                                                      na_matches = "na") %>%
    group_by(Type) %>%
    summarize(n = n()) %>%
    #change NA to uncategorized
    tidyr::replace_na(list(Type = "uncategorized"))

    #add number of type association to current row
    for (t in tmp_merge_ins_sam_typed_summary$Type) {
      tmp_ins <- tmp_ins %>%
      mutate({{t}} := tmp_merge_ins_sam_typed_summary %>%
                 filter(Type == {{t}}) %>%
                 pull(n)
      )
    }
    
    #switch out new for old row in table   
    
    return(tmp_ins)
  }

#get rid of cluster
rm_cluster()

#combine to one dataframe
ait_annotated <- bind_rows(pre_ait_annotated)

head(ait_annotated)


### Annotating with likely/best candidate/s, uncategorized reads and percentages etc. ###
#compute TotalReads per insertion and total (left+right) reads per type per insertion
message("Annotating with best candidates, uncategorized label and percentages:")

for (t in base_all_types) {
    
    tmp_left <- paste0(t, "_L")
    tmp_right <- paste0(t, "_R")
    
    ait_annotated <- ait_annotated %>%
      #add total amount of reads supporting insertion
      mutate(TotalReads = StartReads + EndReads) %>%
      dplyr::relocate(TotalReads, .after = EndReads) %>%
      #add total amount of support reads per (Mu) type
      mutate("{t}_total" := .data[[tmp_left]] + .data[[tmp_right]]) 
  }

#move uncategorized to the end of the df
ait_annotated <- ait_annotated %>%
  relocate(uncategorized, .after = last_col())


#isolate columns for all and type (without uncategorized) and write into vector
all_sss_count_cols <- ait_annotated %>%
  select(ends_with(c("_total", "uncategorized"))) %>%
  names()

type_sss_count_cols <- ait_annotated %>%
  select(ends_with("_total")) %>%
  names()


#start cluster
setup_cluster()

pre_ait_annotated_final <- foreach(row = 1:nrow(ait_annotated)) %dopar% {

  #compute the maximum value found among types/types+uncategorized per insertion
  row_ait_annotated_final <- ait_annotated[row,] %>%
    rowwise() %>%
    #max value for all total and uncategorized columns
    mutate(all_max_value = max(across(ends_with(c("_total", "uncategorized"))))) %>%
    #also compute max value for all types excl. uncategorized
    mutate(type_max_value = max(across(ends_with("_total")))) %>%
    #number of max values for all total and uncategorized columns
    #also, catch edge case when no type has any associated values/reads
    mutate(all_n_max_value = case_when(
           all_max_value > 0 ~ length(which(across(ends_with(c("_total", "uncategorized", "all_max_value")))==all_max_value)),
           TRUE ~ as.integer(0)
      )
    ) %>%
    #number of max values for all total columns excl. uncategorized
    #again, catch edge case when no type has any associated values/reads
    mutate(type_n_max_value = case_when(
           type_max_value > 0 ~ length(which(across(ends_with(c("_total", "type_max_value")))==type_max_value)),
           TRUE ~ as.integer(0)
         )
    )


  #loop over df to identify uncear cases in which the max support for type classification is tied
  #create column with info on whether or not one type is more often than others
  row_ait_annotated_final <- row_ait_annotated_final %>%
    mutate(all_TIES = 
      case_when(
        all_n_max_value > 2 ~ "TIED",
        all_n_max_value == 2 ~ "WINNER",
        TRUE ~ "NO READS"
      ) 
    ) %>%
    mutate(type_TIES = 
      case_when(
        type_n_max_value > 2 ~ "TIED",
        type_n_max_value == 2 ~ "WINNER",
        TRUE ~ "NO READS"
      ) 
    )


  #extract the name of the max value type/uncategorized or infer that it is unclear ("TIED")
  row_ait_annotated_final <- row_ait_annotated_final %>%
    rowwise() %>%
    #probably a better way to do this - case_when() strict; "WINNER" code is a "hidden" list thus we have to add [1] on top
    mutate(all_max_name = 
      case_when(
        all_TIES == "WINNER" ~ names(.[,c(all_sss_count_cols)])[which(across(ends_with(c("_total", "uncategorized")))==all_max_value)][1],
        all_TIES == "TIED" ~ "TIED",
        TRUE ~ "NO READS"
      )
    ) %>%
    mutate(type_max_name = 
      case_when(
        type_TIES == "WINNER" ~ names(.[,c(type_sss_count_cols)])[which(across(ends_with("_total"))==type_max_value)][1],
        type_TIES == "TIED" ~ "TIED",
        TRUE ~ "NO READS"
      )
    ) %>%
    mutate(all_max_name = str_remove(all_max_name, "_total")) %>%
    mutate(type_max_name = str_remove(type_max_name, "_total"))


  #compute stats:
  #percent uncategorized of all reads (TotalReads)
  #percent of the best type/s of all reads (TotalReads) 
  row_ait_annotated_final <- row_ait_annotated_final %>%
    mutate(perc_uncategorized = uncategorized/TotalReads) %>%
    mutate(perc_best_type_of_types = type_max_value/TotalReads)


  #add list elements as df entries per insertion - what are likely candidates (important in unclear cases):
  #excl. uncategorized
  #incl. uncategorized
  row_ait_annotated_final <- row_ait_annotated_final %>%
    mutate(
      all_candidates = case_when(
        type_max_value > 0 ~ list(names(.[,c(all_sss_count_cols)])[which(across(ends_with(c("_total", "uncategorized"))) >= type_max_value)]),
        type_max_value == 0 & all_max_value > 0 ~ list(c("uncategorized")),
        TRUE ~ list(NA)
      )
    ) %>%
    mutate(
      type_candidates = case_when(
        type_max_value > 0 ~ list(names(.[,c(type_sss_count_cols)])[which(across(ends_with("_total")) >= type_max_value)]),
        TRUE ~ list(NA)
      )
    ) %>%
    mutate(all_candidates = list(str_remove(all_candidates, "_total"))) %>%
    mutate(type_candidates = list(str_remove(type_candidates, "_total")))


  return(row_ait_annotated_final)

#close cluster loop
}

#combine to final dataframe
ait_annotated_final <- bind_rows(pre_ait_annotated_final)


#get rid of cluster
rm_cluster()



### write complete table to output ###
data.table::fwrite(ait_annotated_final, 
                   "results/insertions_table_final_te_typed/complete_all_identified_insertions.csv", 
                   row.names=F)


short_ait_annotated_final <- ait_annotated_final %>%
  select(Chr, InsertionStart, InsertionEnd, Sample, StartReads, EndReads,
         perc_uncategorized, perc_best_type_of_types, all_candidates, type_candidates)

data.table::fwrite(short_ait_annotated_final,
                   "results/insertions_table_final_te_typed/short_all_identified_insertions.csv",
                   row.names=F)



  ### create assessment table for collaborators ###
  #- simple extension of insertion table with of te typing information
  #- reduced information (more to the point)

  #check the type of input table
  #three possibilites:
  #1 germinal-annotated
  #2 all-annotated
  #3 germinal/all UNannotated
  #data.table fwrite can deal with coluns containing lists (elements become seperated by "|")

#  if ( all(c("GeneID", "stock") %in% names(ait_annotated)) ) {
#    out <- ait_annotated %>%
#      select(GeneID, Chr, GeneStart, GeneEnd, Sample, InsertionStart, InsertionEnd, StartReads, EndReads, Gene_length, stock,
#             perc_uncategorized, perc_best_type_of_types, all_candidates, type_candidates)
#    data.table::fwrite(out,
#                       paste0("results/insertions_table_final_te_typed/short_", insertion_table_name, ".csv"), 
#                       row.names=F)
#  } else if ( "GeneID" %in% names(ait_annotated) && !("stock" %in% names(ait_annotated)) ) {
#    out <- ait_annotated %>%
#      select(GeneID, Chr, GeneStart, GeneEnd, Sample, InsertionStart, InsertionEnd, StartReads, EndReads, Gene_length,
#             perc_uncategorized, perc_best_type_of_types, all_candidates, type_candidates)
#    data.table::fwrite(out,
#                       paste0("results/insertions_table_final_te_typed/short_", insertion_table_name, ".csv"),
#                       row.names=F)
#  } else {
#    out <- ait_annotated %>%
#      select(Chr, Sample, InsertionStart, InsertionEnd, StartReads, EndReads,
#             perc_uncategorized, perc_best_type_of_types, all_candidates, type_candidates)
#      data.table::fwrite(out,
#                         paste0("results/insertions_table_final_te_typed/short_", insertion_table_name, ".csv"),
#                         row.names=F)
#  }


  ### create stats file for library run
  #- insertions categorized
  #- insertions uncategorized
  #- unclear vs clear

  type_candidates_stats <- ait_annotated_final %>%
    group_by(type_candidates) %>%
    summarize(count = n()) %>%
    arrange(desc(count))

  all_candidates_stats <- ait_annotated_final %>%
    group_by(all_candidates) %>%
    summarize(count = n()) %>%
    arrange(desc(count))


  data.table::fwrite(type_candidates_stats,
                     "results/insertions_table_final_te_typed/type_candidate_stats_all_identified_insertions.csv", 
                     row.names=F)

  data.table::fwrite(all_candidates_stats, 
                     "results/insertions_table_final_te_typed/all_candidates_stats_all_identified_insertions.csv", 
                     row.names=F)



  ### create files for insertions with uncategorized reads ###

  #create subset df with only uncategorized insertions / or with a percentage cutoff
  unc_ait_annotated_final <- ait_annotated_final %>%
    filter(type_max_name == "NO READS")

  #perhaps make this a config.yaml parameter
  #unc_ait_annotated <- ait_annotated %>%
  #  filter(perc_uncategorized > 0.5)


  # create empty dataframe/tibble
  all_uncategorized_ins <- tibble(
                           Name = character(),
                           Strand = character(),
                           Chr = character(),
                           InsertionStart = integer(),
                           InsertionEnd = integer(),
                           Sample = character()
                         )

  #loop through all rows
  setup_cluster()

  pre_unc_ins <- foreach(cur_row = 1:nrow(unc_ait_annotated_final)) %dopar% {

    #match Sample with sam object
    tmp_sam <- get(paste0("sam_", unc_ait_annotated_final[cur_row,]$Sample))

    #match sample of insertion file with te typing file
    tmp_type_file <- get(paste0(unc_ait_annotated_final[cur_row,]$Sample, "_type_file"))

    # fuzzyjoin for overlap + end or start needs to match
    tmp_merge_ins_sam <- fuzzyjoin::genome_inner_join(unc_ait_annotated_final[cur_row,],
                                                      tmp_sam, 
                                                      by=c("Chr", "InsertionStart"="Start", "InsertionEnd"="End")
                          ) %>%
    select(Name, Flag, Chr.y, Start, End, Sample, InsertionStart, InsertionEnd) %>%
    dplyr::rename(Chr=Chr.y) %>%
    filter(Start == unc_ait_annotated_final[cur_row,]$InsertionStart | 
           #Start == unc_ait_annotated[cur_row,]$InsertionEnd | 
           #End == unc_ait_annotated[cur_row,]$InsertionStart | 
           End == unc_ait_annotated_final[cur_row,]$InsertionEnd
          ) %>%
    #translate Flag into forward (+) or reverse (-) strand mapping/read
    rowwise() %>%
    mutate(
      Strand = case_when(
        number2binary(Flag, 8)[7] == 1 & number2binary(Flag, 8)[8] == 0 ~ "1",
        number2binary(Flag, 8)[7] == 0 & number2binary(Flag, 8)[8] == 1 ~ "2",
        TRUE ~ "LOST"
      )
    ) %>%
    select(-Flag) %>%
    relocate(Strand, .after=Name)
    


#    all_uncategorized_ins <- bind_rows(all_uncategorized_ins, tmp_merge_ins_sam)
    return(tmp_merge_ins_sam)
  }

#merge to final dataframe
unc_ins <- bind_rows(pre_unc_ins)

#get rid of cluster
rm_cluster()



#split into forward and reverse reads
strand_1_uncategorized_ins <- unc_ins %>%
  filter(Strand == "1")

strand_2_uncategorized_ins <- unc_ins %>%
  filter(Strand == "2")


#create files only with headers
headers_all_uncategorized_ins <- unc_ins %>%
  select(Name)

headers_strand_1_uncategorized_ins <- strand_1_uncategorized_ins %>%
  select(Name)

headers_strand_2_uncategorized_ins <- strand_2_uncategorized_ins %>%
  select(Name)

#write table for all_uncategorized_ins
write.csv(unc_ins, 
          "results/insertions_table_final_te_typed/all_uncategorized_all_identified_insertions.csv", 
          quote=FALSE, 
          row.names=FALSE)

#header files for strand 1 and 2 seperate
write.csv(headers_strand_1_uncategorized_ins, 
          "results/insertions_table_final_te_typed/headers_strand_1_uncategorized_all_identified_insertions.csv", 
          quote=FALSE,
          row.names=FALSE)
write.csv(headers_strand_2_uncategorized_ins, 
          "results/insertions_table_final_te_typed/headers_strand_2_uncategorized_all_identified_insertions.csv", 
          quote=FALSE,
          row.names=FALSE)


#for loop - input_list
#}

message("All done!! ;D")
