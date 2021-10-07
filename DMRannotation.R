library(biomaRt)
library(GenomicRanges)
library(IRanges)
library(Repitools)
library(tidyverse)

#loading the data, turning it into a dataframe
data<- read.csv("DMRsPathologyTg4510.csv")
#taking the chr:position column, and splitting it into Chr and position
seqnames<-data$seqnames
temp <- strsplit(seqnames, "r")
mat  <- matrix(unlist(temp), ncol=2, byrow=TRUE)
df   <- as.data.frame(mat)
colnames(df)<- c("ch", "chromosome")
#taking the Chr1, and turning it into chr and 1 for annotation
#adding the start and end position colunmns
df["start_position"]<-data["start"]
df["end_position"]<-data["end"]
df["strand"]<-data["strand"]
colnames(df)<-c("ch", "chromosome", "start_position", "end_position", "strand")

#a small loop that cheks for M chromosome and replaces it with MT
for (i in 1:nrow(df)){
  if (df[i,2] == "M"){
    df[i,2] <- "MT"
  }
  
}

#using biomart the get he mouse genome and printing out the chromosome positions of each gene
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
attributes <- c( "chromosome_name",
                 "start_position", "end_position", "external_gene_name", "name_1006")
genes <- getBM(attributes=attributes , mart=mart)

#using genomic ranges to turn the data frames into GRanges, so that I can use the findoverlap feature
#In turn allowing me to find which positions are in which range, allowing me to allocate them to a gene
gr <- GRanges(seqnames = genes$chromosome_name, ranges = IRanges(start = genes$start_position, end = genes$end_position), mcol = genes$external_gene_name)
gr2 <- GRanges(seqnames = df$chromosome, ranges = IRanges(start = as.numeric(df$start_position), end = as.numeric(df$end_position)))
olap <- findOverlaps(gr2, gr, type = "within")
gr2.matched <- gr2[queryHits(olap)]
mcols(gr2.matched) <- cbind.data.frame(
  mcols(gr2.matched),
  mcols(gr[subjectHits(olap)]))


#find the distance to the nearest gene
DistanceToGene<-as.data.frame(distanceToNearest(gr2,gr))


#using the repitools package I converted the GRanges table into a dataframe, so that i can use the Merge function to annotate it back to the original dataframe
#in the future this could be replaced with a loop that checks the distance to the gene
gr3<- annoGR2DF(gr2.matched)
colnames(gr3)<- c("chromosome", "start_position", "end_position","Width" , "GENE ID")
gr4<-merge(gr3, df, by = intersect(names(gr3), names(df)), all.y = TRUE)
gr4<-gr4[-c(4,6,7)]

gr5 <- gr4 %>%
  group_by(start_position) %>%
  # Replace NA with ""
  mutate_all(funs(replace(., is.na(.), ""))) %>%
  # Combine all strings
  summarize_all(funs(toString(unique(.)))) %>%
  # Replace the strings ended with ", "
  mutate_all(funs(str_replace(., ", $", ""))) %>%
  ungroup()

order.pop<-order(gr5$chromosome)
gr6<-gr5[order.pop,]

data["Gene"]<-gr6["GENE ID"]
data["Distance to Gene"]<-DistanceToGene$distance
data[is.na(data)] = ""


#annotate the nearest gene
nearest<-nearest(gr2,gr)
for (i in 1:length(nearest)){
  j<-nearest[i] 
  data$NearestGene[i] = genes[j,4]
}

#exporting the dataframe with Gene IDs as a csv file
write.csv(data,"DMRsPathologyTg4510J20ANNOTATIONS.csv")
