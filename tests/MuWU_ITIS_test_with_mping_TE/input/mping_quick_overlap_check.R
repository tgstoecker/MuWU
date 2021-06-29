library(dplyr)
library(fuzzyjoin)
library(IRanges)


itis_mping <- read.csv("input/test.mping.filtered.bed", header=FALSE, sep="\t")
itis_mping_short <- itis_mping[, c(1,2,3)]
names(itis_mping_short) <- c("Chr" , "Start", "End")


MuWU_karma <- read.csv("results/insertions_table_final/all_identified_insertions.csv", header=TRUE, sep=",")


message("ITIS finds 6 high-confidence mping insertions in this dataset:")
message("")
itis_mping

message("")
message("MuWU intersect with ITIS results:")

message("")
joined_sets <- fuzzyjoin::genome_join(itis_mping_short, MuWU_karma, by=c("Chr", "Start", "End"), mode="inner")

joined_sets <- joined_sets %>%
                 dplyr::select(-Chr.y) %>%
                 dplyr::rename(Chr = Chr.x,
                               Start.ITIS = Start.x, 
                               End.ITIS =  End.x,
                               Start.MuWU = Start.y,
                               End.MuWU =  End.y,
                               MuWU.StartReads = StartReads,
                               MuWU.EndReads = EndReads)

joined_sets

message("")
message("Top 8 MuWU insertions contain the 6 high confidence candidates from ITIS:")

message("")
MuWU_karma_top <- MuWU_karma %>%
  mutate(Support = StartReads + EndReads) %>%
  arrange(desc(Support)) %>%
  head(n = 8)
MuWU_karma_top

#message("The overlap - 8 topHits MuWU & 6 high confidence candidates from ITIS:")
#joined_sets_top <- fuzzyjoin::genome_join(itis_mping_short, MuWU_karma_top, by=c("Chr", "Start", "End"), mode="inner")
#joined_sets_top

