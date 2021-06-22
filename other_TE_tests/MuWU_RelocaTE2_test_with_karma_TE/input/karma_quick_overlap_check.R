library(fuzzyjoin)
library(IRanges)


rTE2_karma <- read.csv("input/karma_TE_relocaTE2_results.gff", header=FALSE, sep="\t")
rTE2_karma_short <- rTE2_karma[, c(1,4,5)]
names(rTE2_karma_short) <- c("Chr" , "Start", "End")


MuWU_karma <- read.csv("results/insertions_table_final/all_identified_insertions.csv", header=TRUE, sep=",")


message("RelocaTE2 finds 23 karma insertions in this dataset:")
message("")
rTE2_karma

message("")
message("Overlapping with MuWU results:")
message("")
joined_sets <- fuzzyjoin::genome_join(rTE2_karma_short, MuWU_karma, by=c("Chr", "Start", "End"), mode="inner")
joined_sets

message("")
message("Consider, that MuWU would normally be used after a targeted PCR approach.")
message("We nevertheless identify ", length(unique(joined_sets$Start.x)), " of the 23 karma insertions as identified by RelocaTE2")
