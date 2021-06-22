library(dplyr)
library(fuzzyjoin)
library(IRanges)


itis_mping <- read.csv("input/test.mping.filtered.bed", header=FALSE, sep="\t")
itis_mping_short <- itis_mping[, c(1,2,3)]
names(itis_mping_short) <- c("Chr" , "Start", "End")


MuWU_karma <- read.csv("results/insertions_table_final/all_identified_insertions.csv", header=TRUE, sep=",")


message("ITIS finds 6 mping insertions in this dataset.")
message("")


message("")
message("Overlapping with MuWU results:")
message("")
itis_mping

message("")
joined_sets <- fuzzyjoin::genome_join(itis_mping_short, MuWU_karma, by=c("Chr", "Start", "End"), mode="inner")
joined_sets

message("")
message("Consider, that MuWU would normally be used after a targeted PCR approach.")
message("We nevertheless identify ", length(unique(joined_sets$Start.x)), " of the 6 high confidence mping insertions as identified by ITIS")

message("")
message("Also, compared to the test with RelocaTE2 the data test data here is especially enriched for one TE - mping.")
message("Looking at the top supported candidates in our MuWU analysis - StartEnd + EndReads")
message("")

MuWU_karma_top <- MuWU_karma %>%
  mutate(Support = StartReads + EndReads) %>%
  arrange(desc(Support)) %>%
  head(n = 8)
MuWU_karma_top

message("The overlap - 8 topHits MuWU & 6 high confidence candidates from ITIS:")
joined_sets_top <- fuzzyjoin::genome_join(itis_mping_short, MuWU_karma_top, by=c("Chr", "Start", "End"), mode="inner")
joined_sets_top

message("")
message("MuWU overlaps quite well with the results of ITIS.")
message("In the top 8 MuWU candidates we find all 6 ITIS top candidates.")
message("(Note that ITIS and MuWU insertions positions slightly differently - MuWU presents TSDs.)")
