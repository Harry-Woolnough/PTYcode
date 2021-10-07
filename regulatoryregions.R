library(biomaRt)
library(GenomicRanges)
library(IRanges)
library(Repitools)
library(dplyr)
library(tidyverse)
library(rtracklayer) 


# using biomart to access the regulatory regions dataset from ensembl
ensembl <- useEnsembl(biomart = "regulation", dataset = "mmusculus_regulatory_feature")
attributes <- c( "chromosome_name", "chromosome_start", "chromosome_end", "feature_type_name", "regulatory_stable_id", "so_accession")
regulatoryregions <- getBM(attributes=attributes , mart=ensembl)

# loading the annotated dataset, splitting the chromosome and such
dataset<-read.csv("annotations/DMPsInteractionModel_J20_sig_genotypeANNOTATIONS.csv")
# taking the chr:position column, and splitting it into Chr and position
start_positions<-dataset$X
temp <- strsplit(start_positions, ":")
mat  <- matrix(unlist(temp), ncol=2, byrow=TRUE)
df   <- as.data.frame(mat)
colnames(df)<- c("chromosome", "start_position")
# taking the Chr1, and turning it into chr and 1 for annotation
positions<-df$chromosome
temp <- strsplit(positions, "r")
mat2  <- matrix(unlist(temp), ncol=2, byrow=TRUE)
df2   <- as.data.frame(mat2)
# adding the start and end position colunmns
df2["start_position"]<-df["start_position"]
df2["end_position"]<-df["start_position"]
colnames(df2)<-c("ch", "chromosome", "start_position", "end_position")

for (i in 1:nrow(df2)){
  if (df2[i,2] == "M"){
    df2[i,2] <- "MT"
  }
  
}
# using granges to find which of the chromosome positions overlap
gr <- GRanges(seqnames = regulatoryregions$chromosome_name, ranges = IRanges(start = regulatoryregions$chromosome_start, end = regulatoryregions$chromosome_end), mcol = regulatoryregions$feature_type_name)
gr2 <- GRanges(seqnames = df2$chromosome, ranges = IRanges(start = as.numeric(df2$start_position), end = as.numeric(df2$end_position)))
olap <- findOverlaps(gr2, gr, type = "within", ignore.strand=TRUE)
gr2.matched <- gr2[queryHits(olap)]
mcols(gr2.matched) <- cbind.data.frame(
  mcols(gr2.matched),
  mcols(gr[subjectHits(olap)]))

# turn the GRanges into a dataframe, so that it can be merged with the original data
gr3<- annoGR2DF(gr2.matched)
colnames(gr3)<- c("chromosome", "start_position", "end_position","Width" , "regulatory_region")
gr4<-merge(gr3, df2, by = intersect(names(gr3), names(df2)), all.x = TRUE, all.y = TRUE)

# get rid of unecessary columns
dataset["start_position"]<-df["start_position"]
dataset<-merge(x = dataset, y = gr4, by.x = "start_position", by.y = "start_position", all.y = TRUE)
dataset<-dataset[-c(1,2,25,26,27,29)]
dataset[is.na(dataset)] = ""

# used where there are multiple genes/regulatory regions annotated to one position, to put them all in one cell as a list
dataset2 <- dataset %>%
  group_by(X) %>%
  # Replace NA with ""
  mutate_all(funs(replace(., is.na(.), ""))) %>%
  # Combine all strings
  summarize_all(funs(toString(unique(.)))) %>%
  # Replace the strings ended with ", "
  mutate_all(funs(str_replace(., ", $", ""))) %>%
  ungroup()

write.csv(dataset2,"final/DMPsInteractionModel_J20_sig_genotypeFINAL.csv")