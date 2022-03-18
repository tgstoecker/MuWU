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
all_types <- snakemake@params[["all_types"]]



#te annotated set of all identified insertions
complete_all_final <- read.csv(snakemake@input[["complete_all_final"]], header = TRUE) %>% as_tibble() %>% mutate(Chr = as.character(Chr))

#other three tables not yet subtype/element annotated
germinal_insertions <- read.csv(snakemake@input[["germinal_identified_insertions"]], header = TRUE) %>% as_tibble() %>% mutate(Chr = as.character(Chr))
all_insertions_annotated <- read.csv(snakemake@input[["all_identified_insertions_annotated"]], header = TRUE) %>% as_tibble() %>% mutate(Chr = as.character(Chr))
germinal_insertions_annotated <- read.csv(snakemake@input[["germinal_identified_insertions_annotated"]], header = TRUE) %>% as_tibble() %>% mutate(Chr = as.character(Chr))

#using all these columns of the un-/annotated germinal and annotated all tables might be overkill..
#but this way I can spot if sth. upstream destroyed the coordinates read counts

#set DT threads to 1
data.table::setDTthreads(1)


#create my cluster functions
te_typing_cluster_cores <- snakemake@params[["te_typing_cluster_cores"]]

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

  ## some time before next command - making sure all sockets are completely wiped before new cluster start
  Sys.sleep(60)
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






### 1 - Germinal insertions NON-annotated

germinal_ait_final <- dplyr::left_join(germinal_insertions,
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


write.csv(germinal_ait_final,
          "results/insertions_table_final_te_typed/complete_germinal_identified_insertions.csv",
          row.names=F,
          quote=F)


#short version of the table
short_germinal_ait_final <- germinal_ait_final %>%
  select(Chr, InsertionStart, InsertionEnd, Sample, StartReads, EndReads,
         perc_uncategorized, perc_best_type_of_types, all_candidates, type_candidates)

data.table::fwrite(short_germinal_ait_final,
                   "results/insertions_table_final_te_typed/short_germinal_identified_insertions.csv",
                   row.names=F)


### create stats file for library run
#- insertions categorized
#- insertions uncategorized
#- unclear vs clear

type_candidates_stats <- germinal_ait_final %>%
  group_by(type_candidates) %>%
  summarize(count = n()) %>%
  arrange(desc(count))

all_candidates_stats <- germinal_ait_final %>%
  group_by(all_candidates) %>%
  summarize(count = n()) %>%
  arrange(desc(count))


data.table::fwrite(type_candidates_stats,
                   "results/insertions_table_final_te_typed/type_candidate_stats_germinal_identified_insertions.csv", 
                   row.names=F)

data.table::fwrite(all_candidates_stats, 
                   "results/insertions_table_final_te_typed/all_candidates_stats_germinal_identified_insertions.csv", 
                   row.names=F)



### create files for insertions with uncategorized reads ###

#create subset df with only uncategorized insertions / or with a percentage cutoff
germinal_unc_ait_annotated_final <- germinal_ait_final %>%
  filter(type_max_name == "NO READS")


# create empty dataframe/tibble
#germinal_uncategorized_ins <- tibble(
#                           Name = character(),
#                           Strand = character(),
#                           Chr = character(),
#                           InsertionStart = integer(),
#                           InsertionEnd = integer(),
#                           Sample = character()
#                         )

#loop through all rows
setup_cluster()

pre_unc_ins <- foreach(cur_row = 1:nrow(germinal_unc_ait_annotated_final)) %dopar% {

    #match Sample with sam object
    tmp_sam <- get(paste0("sam_", germinal_unc_ait_annotated_final[cur_row,]$Sample))

    #match sample of insertion file with te typing file
    tmp_type_file <- get(paste0(germinal_unc_ait_annotated_final[cur_row,]$Sample, "_type_file"))

    # fuzzyjoin for overlap + end or start needs to match
    tmp_merge_ins_sam <- fuzzyjoin::genome_inner_join(germinal_unc_ait_annotated_final[cur_row,],
                                                      tmp_sam, 
                                                      by=c("Chr", "InsertionStart"="Start", "InsertionEnd"="End")
                          ) %>%
    select(Name, Flag, Chr.y, Start, End, Sample, InsertionStart, InsertionEnd) %>%
    dplyr::rename(Chr=Chr.y) %>%
    filter(Start == germinal_unc_ait_annotated_final[cur_row,]$InsertionStart | 
           #Start == unc_ait_annotated[cur_row,]$InsertionEnd | 
           #End == unc_ait_annotated[cur_row,]$InsertionStart | 
           End == germinal_unc_ait_annotated_final[cur_row,]$InsertionEnd
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

    return(tmp_merge_ins_sam)
  }

#merge to final dataframe
germinal_unc_ins <- bind_rows(pre_unc_ins)


#get rid of cluster
rm_cluster()



#split into forward and reverse reads
strand_1_uncategorized_ins <- germinal_unc_ins %>%
  filter(Strand == "1")

strand_2_uncategorized_ins <- germinal_unc_ins %>%
  filter(Strand == "2")


#create files only with headers
headers_all_uncategorized_ins <- germinal_unc_ins %>%
  select(Name)

headers_strand_1_uncategorized_ins <- strand_1_uncategorized_ins %>%
  select(Name)

headers_strand_2_uncategorized_ins <- strand_2_uncategorized_ins %>%
  select(Name)

#write table for all_uncategorized_ins
write.csv(germinal_unc_ins,
          "results/insertions_table_final_te_typed/all_uncategorized_germinal_identified_insertions.csv", 
          quote=FALSE, 
          row.names=FALSE)

#header files for strand 1 and 2 seperate
write.csv(headers_strand_1_uncategorized_ins, 
          "results/insertions_table_final_te_typed/headers_strand_1_uncategorized_germinal_identified_insertions.csv", 
          quote=FALSE,
          row.names=FALSE)
write.csv(headers_strand_2_uncategorized_ins, 
          "results/insertions_table_final_te_typed/headers_strand_2_uncategorized_germinal_identified_insertions.csv", 
          quote=FALSE,
          row.names=FALSE)







########################################
#%>%
#  #order so that row-column intersections are paired
#  arrange(., Chr, GeneStart, InsertionStart)
