library(biomaRt)
library(GenomicRanges)
library(IRanges)
library(Repitools)
library(dplyr)
library(tidyverse)
# loading the data, turning it into a dataframe
data<- read.csv("starters/DMPsInteractionModel_J20_sig_genotype.csv")
# taking the chr:position column, and splitting it into Chr and position
start_positions<-data$X
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

# a small loop that cheks for M chromosome and replaces it with MT
for (i in 1:nrow(df2)){
  if (df2[i,2] == "M"){
    df2[i,2] <- "MT"
  }
  
}

# using biomart the get he mouse genome and printing out the chromosome positions of each gene
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
attributes <- c( "chromosome_name",
                 "start_position", "end_position", "external_gene_name")
genes <- getBM(attributes=attributes , mart=mart)

# using genomic ranges to turn the data frames into GRanges, so that I can use the findoverlap feature
# In turn allowing me to find which positions are in which range, allowing me to allocate them to a gene
gr <- GRanges(seqnames = genes$chromosome_name, ranges = IRanges(start = genes$start_position, end = genes$end_position), mcol = genes$external_gene_name)
gr2 <- GRanges(seqnames = df2$chromosome, ranges = IRanges(start = as.numeric(df2$start_position), end = as.numeric(df2$end_position)))
olap <- findOverlaps(gr2, gr, type = "within", ignore.strand=TRUE)
gr2.matched <- gr2[queryHits(olap)]
mcols(gr2.matched) <- cbind.data.frame(
  mcols(gr2.matched),
  mcols(gr[subjectHits(olap)]))


# find the distance to the nearest gene
DistanceToGene<-as.data.frame(distanceToNearest(gr2,gr))


# using the repitools package I converted the GRanges table into a dataframe, so that i can use the Merge function to annotate it back to the original dataframe
# in the future this could be replaced with a loop that checks the distance to the gene
gr3<- annoGR2DF(gr2.matched)
colnames(gr3)<- c("chromosome", "start_position", "end_position","Width" , "GENE_ID")
gr4<-merge(gr3, df2, by = intersect(names(gr3), names(df2)), all.x = TRUE, all.y = TRUE)
gr4<-gr4[-c(4,6,7)]

# annotate the nearest gene
nearest<-nearest(gr2,gr)
for (i in 1:length(nearest)){
  j<-nearest[i] 
  data$NearestGene[i] = genes[j,4]
}

data["Distance to Gene"]<-DistanceToGene$distance


data["start_position"]<-df["start_position"]
data<-merge(x = data, y = gr4, by.x = "start_position", by.y = "start_position", all.y = TRUE)
data<-data[-c(1,23,24)]
data[is.na(data)] = ""


dat2 <- data %>%
  group_by(X) %>%
  # Replace NA with ""
  mutate_all(funs(replace(., is.na(.), ""))) %>%
  # Combine all strings
  summarize_all(funs(toString(unique(.)))) %>%
  # Replace the strings ended with ", "
  mutate_all(funs(str_replace(., ", $", ""))) %>%
  ungroup()


# exporting the dataframe with Gene IDs as a csv file
write.csv(dat2,"annotations/DMPsPathology_J20_sig_pathologyANNOTATIONS.csv")
