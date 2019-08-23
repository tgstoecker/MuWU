## This adds gene annotation (and transcript annotation) to the previously created MuSeq .csv table

library(GenomicRanges)
library(ChIPpeakAnno)

options(warn=-1)

## Read table with identified Mu single line Insertions
MuSingle <- read.csv("MuSeq_table/SLI-MuSeq_FGS.csv", header=T)
MuSingle$Insertionsites <- paste(MuSingle$Chr, MuSingle$Start, MuSingle$End)

## Read table with identified Mu single line Insertions
MuAll <- read.csv("MuSeq_table/MuSeq_FGS.csv", header=T)
MuAll$Insertionsites <- paste(MuAll$Chr, MuAll$Start, MuAll$End)

## Read GTF File
setwd("FGS/")  
MY_gtf <- list.files(getwd(), pattern="\\.gtf$") 
FGS <- read.delim(MY_gtf, header=FALSE, comment.char="#")
FGS$gene <- paste(FGS$V1, FGS$V4, FGS$V5)

## Calculate Ranges for Mu single line and genes in the FGS
Ranges <- GRanges(seqnames=MuSingle$Chr, ranges=IRanges(start=MuSingle$Start, end=MuSingle$End))
GTF <- GRanges(seqnames = FGS$V1, ranges=IRanges(start=FGS$V4, end=FGS$V5))

## Calculate Ranges for all Mu insertions in the FGS
Ranges_all <- GRanges(seqnames=MuAll$Chr, ranges=IRanges(start=MuAll$Start, end=MuAll$End))


## Indentification of Insertions that are in exons of genes in the FGS ##
annotatedPeak <- annotatePeakInBatch(Ranges, AnnotationData=GTF, output="inside")
annotatedPeak_all <- annotatePeakInBatch(Ranges_all, AnnotationData=GTF, output="inside")

Insertions_FGS <- as.data.frame(annotatedPeak)
Insertions_FGS <- na.omit(Insertions_FGS)
Insertions_FGS_all <- as.data.frame(annotatedPeak_all)
Insertions_FGS_all <- na.omit(Insertions_FGS_all)


Insertions_FGS$Insertionsites <- paste(Insertions_FGS$seqnames, Insertions_FGS$start, Insertions_FGS$end)
Insertions_FGS$Insertionsites <- gsub("chr", "", Insertions_FGS$Insertionsites)
Insertions_FGS$gene <- paste(Insertions_FGS$seqnames, Insertions_FGS$start_position, Insertions_FGS$end_position)
Insertions_FGS$gene <- gsub("chr", "", Insertions_FGS$gene)

Insertions_FGS_all$Insertionsites <- paste(Insertions_FGS_all$seqnames, Insertions_FGS_all$start, Insertions_FGS_all$end)
Insertions_FGS_all$Insertionsites <- gsub("chr", "", Insertions_FGS_all$Insertionsites)
Insertions_FGS_all$gene <- paste(Insertions_FGS_all$seqnames, Insertions_FGS_all$start_position, Insertions_FGS_all$end_position)
Insertions_FGS_all$gene <- gsub("chr", "", Insertions_FGS_all$gene)


Insertions_FGSinsidemerged <- merge(Insertions_FGS,FGS, by="gene")
Insertions_FGSinside_geneIDs <- Insertions_FGSinsidemerged[,c(1,16,25)]

Insertions_FGS_allinsidemerged <- merge(Insertions_FGS_all,FGS, by="gene")
Insertions_FGS_allinside_geneIDs <- Insertions_FGS_allinsidemerged[,c(1,16,25)]


Mu_single_GeneIds <- merge(MuSingle, Insertions_FGSinside_geneIDs, by="Insertionsites")
Mu_single_GeneIds <- unique(Mu_single_GeneIds[,c(2,3,4,5,6,7,9)])
colnames(Mu_single_GeneIds) <- c("Chromosome", "Start", "End", "Sample", "StartReads", "EndReads", "GeneID")

Mu_all_GeneIds <- merge(MuAll, Insertions_FGS_allinside_geneIDs, by="Insertionsites")
Mu_all_GeneIds <- unique(Mu_all_GeneIds[,c(2,3,4,5,6,7,9)])
colnames(Mu_all_GeneIds) <- c("Chromosome", "Start", "End", "Sample", "StartReads", "EndReads", "GeneID")

setwd("../MuSeq_table")
write.csv(Mu_single_GeneIds, "SLI-MuSeq_FGS_annotated.csv", row.names=F)
write.csv(Mu_all_GeneIds, "MuSeq_FGS_annotated.csv", row.names=F)
