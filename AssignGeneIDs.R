## This adds gene annotation to the previously created MuSeq .csv table

library(GenomicRanges)
library(ChIPpeakAnno)

## Read table with identified Mu single line Insertions
MuSingle <- read.csv("MuSeq_table/SLI_MuSeq_FGS.csv", header=T)
MuSingle$Insertionsites <- paste(MuSingle$Chr, MuSingle$Start, MuSingle$End)

## Read GTF File
setwd("FGS/")  
MY_gtf <- list.files(getwd(), pattern="\\.gtf$") 
FGS <- read.delim(MY_gtf, header=FALSE, comment.char="#")
FGS$gene <- paste(FGS$V1, FGS$V4, FGS$V5)

## Calculate Ranges for Mu single line and genes in the FGS
Ranges <- GRanges(seqnames=MuSingle$Chr, ranges=IRanges(start=MuSingle$Start, end=MuSingle$End))
GFF <- GRanges(seqnames = FGS$V1, ranges=IRanges(start=FGS$V4, end=FGS$V5))

## Indentification of Insertions that are in exons of genes in the FGS ##
annotatedPeak <- annotatePeakInBatch(Ranges, AnnotationData=GFF, output="inside")

Insertions_FGS <- as.data.frame(annotatedPeak)
Insertions_FGS <- na.omit(Insertions_FGS)

Insertions_FGS$Insertionsites <- paste(Insertions_FGS$seqnames, Insertions_FGS$start, Insertions_FGS$end)
Insertions_FGS$Insertionsites <- gsub("chr", "", Insertions_FGS$Insertionsites)
Insertions_FGS$gene <- paste(Insertions_FGS$seqnames, Insertions_FGS$start_position, Insertions_FGS$end_position)
Insertions_FGS$gene <- gsub("chr", "", Insertions_FGS$gene)

Insertions_FGSinsidemerged <- merge(Insertions_FGS,FGS, by="gene")
Insertions_FGSinside_geneIDs <- Insertions_FGSinsidemerged[,c(1,16,25)]

Mu_single_GeneIds <- merge(MuSingle, Insertions_FGSinside_geneIDs, by="Insertionsites")
Mu_single_GeneIds <- unique(Mu_single_GeneIds[,c(2,3,4,5,6,7,9)])
colnames(Mu_single_GeneIds) <- c("Chromosome", "Start", "End", "Sample", "StartReads", "EndReads", "GeneID")

setwd("../MuSeq_table")
write.csv(Mu_single_GeneIds, "SLI_MuSeq_FGS_annotated.csv", row.names=F)
