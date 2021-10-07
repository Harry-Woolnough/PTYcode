#Number of sites vs number of annotated genes
#Number of annotated genes in gene body
#Proportion of DMPs within 1kb of a gene
#Proportion of DMPs within 100 kb of a gene

library(ggplot2)
library(scales)
library("ggVennDiagram")
library(ggpubr)
options(scipen = n)

#a function which reads the file and saves the nearest gene column, and figures out how many unique annotated genes there are. 
#Also works out the proportion of genes distances that is under <100 and <1000kb 
SitesVsGenes<- function(file){
data<-read.csv(file)
neargenes<-unique(data$Nearest_Gene)
ngenes<- unique(data$Gene)
pdata<-data$Distance_to_gene
pdata1000<-length(pdata[pdata<1000 & pdata !=0])
p1000<-pdata1000/nrow(data)
pdata100<-length(pdata[pdata<100 & pdata !=0])
p100<-pdata100/nrow(data)
pdatamore<-length(pdata[pdata>100000])
pmore<-pdatamore/nrow(data)
pdata100000<-length(pdata[pdata<100000 & pdata !=0])
p100000<-pdata100000/nrow(data)

#Prints the data into the console
print(paste("Number of annotated genes:",length(neargenes), "| Number of sites:",nrow(data)  ))
print(paste("Number of annotated genes in gene body:", length(ngenes)))
print(paste("Proprtion of genes within 1000kb: ", p1000))
print(paste("Proprtion of genes within 100kb: ", p100))

#creating a bar chart to compare number of annotated genes and number of sites
par(mfrow=c(2,2))

bardata<-data.frame(
  value= c(length(neargenes), nrow(data)),
  group= c("Annotated Genes", "Sites")
)


bc<-ggplot(bardata,aes(x = group, y = value, fill = group)) + geom_bar(stat="identity") + ggtitle("Number of annotated genes vs Number of sites") + 
  geom_text(aes(label = value))

#using the distance proportions a data frame is made for the piechart ,which is made with ggplot
piedata<-data.frame(
  value = c(p100,(p1000-p100), (p100000-p1000), pmore ),
  group = c("n<100bp", "100bp<n<1kb", "1kb<n<100Kb","n>100000")
)

pc<-ggplot(piedata, aes(x="", y=value, fill = group))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + geom_text(aes(label = percent(value)), size=3, position=position_stack(vjust=0.5))+
theme_void() + ggtitle("Proportion of DMPs within 100kb and 1000kbs of a gene")

#using the ggarrange package to combine more than on plot into one image
ggarrange(bc,pc, labels = c("A", "B"), ncol=2, nrow= 1)
  }
