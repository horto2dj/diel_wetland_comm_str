library(permute)
library(vegan)
library(MASS)
library(Hmisc)
library(ggbiplot)
library(nlme)
library(ggplot2)
library(ggrepel)
library(psych)
library(xlsx)
library(gridExtra)
library(phyloseq)
library(WGCNA)
library(survival)
library(multcomp)
require(DESeq2)
require(ggmap)
require(tidyverse)
library(RColorBrewer)

memory.size(max = TRUE)
memory.limit(size = 15000)
memory.limit()

setwd("~/Grad School/Central Michigan/Dissertation/Diel/2015 DNA study/stats")


##################
### WGCNA MM_B ###
##################

NUT <- read.table("NUT_MM_B.txt", header = TRUE, row.names = 1, sep ="\t")
OTU <- read.table("OTU_MM_B.txt", header = TRUE, row.names = 1, sep="\t")
TAX <- read.table("diel_tax_singdoubrem_cols.txt", header = TRUE, row.names = 1, sep = "\t")

rownames(OTU) <- paste0("OTU", 1:nrow(OTU))
rownames(OTU)

TAX <- as.matrix(TAX, rownames.force = NA)
rownames(TAX) <- paste0("OTU", 1:nrow(TAX))
rownames(TAX)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)

physeq = phyloseq(OTU,TAX)

META = sample_data(NUT)
rownames(META) <- sample_names(physeq)

META = sample_data(META)

# get all the data in a phyloseq instance, or whatever
ALL = phyloseq(OTU,TAX,META)
ALL

### Prepping data to remove rare or erronious OTUs #####

# We need to de-noise the data by plotting the number of reads on a curve and look for the inflection point

at.least.n.in.m <- function(x, n, m){
  all(x[x>0]>=n)&length(x[x>0])>=m
}
counts<- rep(0,10)
for (i in 1:length(counts)){
  rows.to.keep<- apply(otu_table(ALL, taxa_are_rows = TRUE), 1,at.least.n.in.m, n=i, m=2)
  counts[i]<-sum(rows.to.keep)
}

plot(1:10, counts, xlab= 'Min sequences in 2 samples', ylab= 'Number of taxa remaining')

### Inflection point was 2 ###

# Filter taxa that arent seen more than twice in greater than 25% of the data.
CUT2<-filter_taxa(ALL, function(x) sum(x > 2) > (0.25*length(x)), TRUE)
CUT2
write.csv(tax_table(CUT2), "MM_B_seqnames.csv")
write.csv(otu_table(CUT2), "MM_B_seqOTU.csv")

### Deseq2 Normalization
### convert the counts to integer mode and group data to normalize together(e.g. region)
dds_WGCNA <- phyloseq_to_deseq2(CUT2, ~ Diel)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans_WGCNA = apply(counts(dds_WGCNA), 1, gm_mean)
dds_WGCNA = estimateSizeFactors(dds_WGCNA, geoMeans=geoMeans_WGCNA)
dds_WGCNA = estimateDispersions(dds_WGCNA)

### variance stabilizing transformation
#Make it a OTU again for phyloseq (no tax or sample data)
vst_WGCNA <- varianceStabilizingTransformation(dds_WGCNA, blind=FALSE)
vstMat_WGCNA <- assay(vst_WGCNA)
vstMat_WGCNA[vstMat_WGCNA<0]<-0
vst.otu.WGCNA <- otu_table(vstMat_WGCNA, taxa_are_rows=TRUE)
write.csv(otu_table(vst.otu.WGCNA), "vst.otu.WGCNA.MM_B.csv")

#### BEGIN WGCNA #####

#### Maybe try this next time... ####
vst <- read.csv("vst.otu.WGCNA.MM_B.csv")
dim(vst)
names(vst)
otu_WGCNA <- as.data.frame(t(vst))
names(otu_WGCNA) = vst$X
rownames(otu_WGCNA) = names(vst)
dim(otu_WGCNA)
names(otu_WGCNA)


NMDS_META <- read.table("NUT_MM_B.txt", header = T, row.names = 1)


# Because we are using RStudio, we have to disable threads
# As consequence of this, maybe it would be better to do this step in regular ol' R
disableWGCNAThreads()
options(stringsAsFactors = FALSE)
## Identify beta to ensure scale free topology
powers = c(seq(4,10,by=1), seq(12,20, by=2))

pst <- pickSoftThreshold(otu_WGCNA, powerVector=powers, blockSize = 3767, verbose=2)

# Plot the results of soft thresholds if you wish:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(pst$fitIndices[,1], pst$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(pst$fitIndices[,1], pst$fitIndices[,5], labels=powers, cex=cex1,col="red")

#check the powerEstimate to choose a threshold
pst

# Create adjacency matrix by raising OTU matrix by beta and identify subnetworks (modules)
otu_WGCNA2 <- as.matrix(otu_WGCNA[2:13,1:3767])
mode(otu_WGCNA2)
class(otu_WGCNA2) <- "numeric"
mode(otu_WGCNA2)

# Check that the network ensures scale-free topology at that power
# R should be close to 1 (R > 0.8, I believe), should see a straight line.
##### scaleFreePlot #####
# here we define the adjacency matrix using soft thresholding with beta=12
ADJ1=abs(cor(otu_WGCNA2,use="p"))^12
# When you have relatively few genes (<5000) use the following code
#k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=otu_WGCNA2,power=12)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
scaleFreeFitIndex(k)

#R^2 of 0.84, this suggests we meet the assumption of scale-free topol.

# power of 12 chosen based on powerEstimate from 'pst'
net = blockwiseModules(otu_WGCNA2, power=12, minModuleSize=30, maxBlockSize = 3767,
                       corType = "pearson", saveTOMs = TRUE, 
                       saveTOMFileBase = "blockwiseTOM", pamStage=FALSE, verbose=5)
# Plot the dendrogram
moduleLabels = net$colors
moduleColors = net$colors
MEs = net$MEs
geneTree = net$dendrograms[[1]]
pdf("plotDendro_MM_B.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
# Save data

# Identify Eigenvalues for subnetworks by sample
nPop = ncol(otu_WGCNA2)
nSamples = nrow(otu_WGCNA2)
MEsO = moduleEigengenes(otu_WGCNA2, moduleColors)$eigengenes
MEs = orderMEs(MEsO)
save(MEs, moduleLabels, moduleColors, geneTree,file = "Module-networkConstruction-auto.RData")
# Save data
write.csv(file="Module_eigen_values_MM_B.csv",MEs)
write.csv(file="Module_composition_MM_B.csv",net$colors)

##Correlate Eigenvalues to metadata and create heatmap

moduleTraitCor = cor(MEs, NMDS_META[,c(5:14)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), " (",signif(moduleTraitPvalue, 1), ")"
                   , sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf("Correlation_MM_B.pdf",width=12,height=8)
par(mar = c(12, 12, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(NMDS_META[,c(5:14)]),yLabels = names(MEs)
               ,ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50)
               ,textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,
               cex.lab = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()


## Now make a plot for specific module <-> trait (metadata component) pairings
# This allows us to explore the structure of submodule OTU correlations with a given metadata component
# Here we will use "env variable" as trait and "color" as module
# First get the links between all modules and this trait

###DO-darkred
parameter<-"DO"
weight <- as.data.frame(NMDS_META[,parameter])
names(weight) = "DO"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(otu_WGCNA2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(otu_WGCNA2, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
# Then look at the specific module of interest
module<-"darkred"
column = match(module, modNames);
moduleGenes = moduleColors==module
pdf(paste(module,"-vs-",parameter,".pdf"))
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = paste("Population significance for ",parameter),main = paste("Module membership vs. population significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

##Incorporate OTU taxonomy info to identify specific OTUs within a module, coalesce correlative data by OTU for later use in network analysis
annot2 = read.csv("MM_B_seqnames.csv", colClasses = "character", header = TRUE)
colnames(annot2)[1] <- "OTU"
dim(annot2)
names(annot2)
probes2 = names(otu_WGCNA)
probes2annot2 = match(probes2, annot2$OTU)
# The following is the number or probes without annotation:
sum(is.na(probes2annot2))
# Should return 0.
# Create the starting data frame
geneInfo0 = data.frame(otu_orig = probes2,
                       OTU = annot2$OTU[probes2annot2],
                       Domain = annot2$Domain[probes2annot2],
                       Phylum = annot2$Phylum[probes2annot2],
                       Class = annot2$Class[probes2annot2],
                       Order = annot2$Order[probes2annot2],
                       Family = annot2$Family[probes2annot2],
                       Genus = annot2$Genus[probes2annot2],
                       names = annot2$OTU[probes2annot2],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance according to metadata component
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.DO)); # <- change geneInfo$'X'
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = paste("OTUInfo_MM_B_",module))

# Get the similarity between eigenvalues and weight
MET = orderMEs(cbind(MEs, weight))
pdf("Adjacency_trait_ME_MM_B_DO.pdf") # <- This changes with env. variable
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene dendrogram with ",parameter), marDendro = c(0,4,2,0),plotHeatmaps = FALSE)
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene adjacency heatmap with ",parameter), marHeatmap = c(3,4,2,2),plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

################################################################################
#           Moving to machine learning PLS and VIP scores                      #
################################################################################

library("pls")
source("VIP.R")
# metadata component (e.g. S, designated "weight" here as above) is the same as before, we just replicate the row names for pls
#the below r2 threshold changes with user input, good correlation above "0.x"
#this is important, because VIP scores don't really mean anything without a good correlation
th_r2<-0.3
subnetwork<-otu_WGCNA[2:13,moduleGenes]
row.names(weight)<-row.names(subnetwork)
weight <- as.matrix(weight)
subnetwork <- as.matrix(subnetwork)
class(weight) <- "numeric"
class(subnetwork) <- "numeric"
pls_result<-plsr(weight ~ subnetwork, validation="LOO",method="oscorespls")

r2_vector<-R2(pls_result)
max<-0
max_comp<--1
for (j in 1:length(r2_vector$val)){
  if(r2_vector$val[j]>th_r2){         # We will only look at the PLS if the correlation is better than 0.3
    if(r2_vector$val[j]>max){
      max<-r2_vector$val[j]
      max_comp<-r2_vector$comp[j]
    }
  }
}
print(paste(" the max r2 is ",max," corresponding to comp ",max_comp,sep="",".pdf"))

if(max==0){
  print ("No good correlation, we stop here")
} else{
  print("Good correlation, we check the VIP!")
}
# Checking the VIP
output<-paste("VIP_values_with_",parameter,sep="")
vip_result<-VIP(pls_result)
vip_components<-sort(vip_result[max_comp,],decreasing=TRUE)[1:41] # <- this value changes with subnetwork, check dim of 'subnetwork' file
for (i in 1:41){ # <- this value changes as the line above
  cat(paste("Rank ",i," we have ",names(vip_components[i])," with a VIP of ",vip_components[i],"\n",sep=""),file=output,append=TRUE)
}
weight_2 <- as.data.frame(weight[!is.na(weight)])
df<-data.frame(x=weight_2,y=pls_result$validation$pred[,,max_comp])
colnames(df)<-c("x","y")
pdf(paste("measured_vs_predicted_MM_B_",module,"-vs-",parameter,".pdf"))
ggplot(data=df) + geom_point(aes(x=x,y=y)) + geom_smooth(aes(x=x,y=y),method=lm) + xlab("Measured") + ylab("Predicted") + ggtitle(paste("Comparison of ",parameter," measured vs predicted for module ",module)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
dev.off()
# Establish the correlation between predicted and modeled
# This is the data to report with the figure (R2, CI, signif, etc.)
cor.test(df$x,df$y)

## Identify node centrality based on co-occurence data for each OTU in the module

TOM = TOMsimilarityFromExpr(otu_WGCNA2, power = 12, corType = "pearson");
# Select submodule of interest based on high correlation and signficance
module<-"darkred"; # <- this changes with module color being currently explored
# Select module probes
probes = names(otu_WGCNA)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
write.csv(modTOM, paste("nodeconnections_MM_B_",parameter,".csv"))

# Number of cells above 0.25 threshhold <-- this number is flexible and should change with your data
x<- as.data.frame(rowSums(modTOM > 0.1))
write.csv(x, paste("nodes_MM_B_",parameter,".csv"))


# make scatter hive plots
# You will need to make the Nodeworksheet by combining OTUinfo table for submodule of interest, VIP socres, and Nodeconnections. See workflow for more details.

hive_MM_B<- read.csv("DO_positive_darkred_MM_B_nodeworksheet.csv", header=T) #(Figure 5-8)
hive_MM_B$OTU<-factor(hive_MM_B$OTU, levels = unique(hive_MM_B$OTU))
MM_B_plot <- ggplot(hive_MM_B, aes(x= hive_MM_B$connectivity, y= hive_MM_B$GS.DO)) +
  geom_point(aes(size = hive_MM_B$VIP, colour = hive_MM_B$Phylum, alpha = 0.5)) +
  scale_size_area(max_size= 10) +
  scale_color_brewer(type = qual, palette = "Set1", direction = -1) +
  theme_bw() +
  scale_alpha(guide=FALSE) +
  geom_text_repel(aes(label = hive_MM_B$names), force = 3, size = 3) +
  labs(x="Node centrality", y="Correlation to DO", color="Phylum",
       size = "VIP")
MM_B_plot


## The darkred plot is also significantly related to Temp, so we'll try that too


###DO-darkred
parameter<-"Temp"
weight <- as.data.frame(NMDS_META[,parameter])
names(weight) = "Temp"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(otu_WGCNA2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(otu_WGCNA2, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
# Then look at the specific module of interest
module<-"darkred"
column = match(module, modNames);
moduleGenes = moduleColors==module
pdf(paste(module,"-vs-",parameter,".pdf"))
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = paste("Population significance for ",parameter),main = paste("Module membership vs. population significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

##Incorporate OTU taxonomy info to identify specific OTUs within a module, coalesce correlative data by OTU for later use in network analysis
annot2 = read.csv("MM_B_seqnames.csv", colClasses = "character", header = TRUE)
colnames(annot2)[1] <- "OTU"
dim(annot2)
names(annot2)
probes2 = names(otu_WGCNA)
probes2annot2 = match(probes2, annot2$OTU)
# The following is the number or probes without annotation:
sum(is.na(probes2annot2))
# Should return 0.
# Create the starting data frame
geneInfo0 = data.frame(otu_orig = probes2,
                       OTU = annot2$OTU[probes2annot2],
                       Domain = annot2$Domain[probes2annot2],
                       Phylum = annot2$Phylum[probes2annot2],
                       Class = annot2$Class[probes2annot2],
                       Order = annot2$Order[probes2annot2],
                       Family = annot2$Family[probes2annot2],
                       Genus = annot2$Genus[probes2annot2],
                       names = annot2$OTU[probes2annot2],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance according to metadata component
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Temp)); # <- change geneInfo$'X'
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = paste("OTUInfo_MM_B_",module))

# Get the similarity between eigenvalues and weight
MET = orderMEs(cbind(MEs, weight))
pdf("Adjacency_trait_ME_MM_B_DO.pdf") # <- This changes with env. variable
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene dendrogram with ",parameter), marDendro = c(0,4,2,0),plotHeatmaps = FALSE)
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene adjacency heatmap with ",parameter), marHeatmap = c(3,4,2,2),plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

################################################################################
#           Moving to machine learning PLS and VIP scores                      #
################################################################################

library("pls")
source("VIP.R")
# metadata component (e.g. S, designated "weight" here as above) is the same as before, we just replicate the row names for pls
#the below r2 threshold changes with user input, good correlation above "0.x"
#this is important, because VIP scores don't really mean anything without a good correlation
th_r2<-0.3
subnetwork<-otu_WGCNA[2:13,moduleGenes]
row.names(weight)<-row.names(subnetwork)
weight <- as.matrix(weight)
subnetwork <- as.matrix(subnetwork)
class(weight) <- "numeric"
class(subnetwork) <- "numeric"
pls_result<-plsr(weight ~ subnetwork, validation="LOO",method="oscorespls")

r2_vector<-R2(pls_result)
max<-0
max_comp<--1
for (j in 1:length(r2_vector$val)){
  if(r2_vector$val[j]>th_r2){         # We will only look at the PLS if the correlation is better than 0.3
    if(r2_vector$val[j]>max){
      max<-r2_vector$val[j]
      max_comp<-r2_vector$comp[j]
    }
  }
}
print(paste(" the max r2 is ",max," corresponding to comp ",max_comp,sep="",".pdf"))

if(max==0){
  print ("No good correlation, we stop here")
} else{
  print("Good correlation, we check the VIP!")
}
# Checking the VIP
output<-paste("VIP_values_with_",parameter,sep="")
vip_result<-VIP(pls_result)
vip_components<-sort(vip_result[max_comp,],decreasing=TRUE)[1:41] # <- this value changes with subnetwork, check dim of 'subnetwork' file
for (i in 1:41){ # <- this value changes as the line above
  cat(paste("Rank ",i," we have ",names(vip_components[i])," with a VIP of ",vip_components[i],"\n",sep=""),file=output,append=TRUE)
}
weight_2 <- as.data.frame(weight[!is.na(weight)])
df<-data.frame(x=weight_2,y=pls_result$validation$pred[,,max_comp])
colnames(df)<-c("x","y")
pdf(paste("measured_vs_predicted_MM_B_",module,"-vs-",parameter,".pdf"))
ggplot(data=df) + geom_point(aes(x=x,y=y)) + geom_smooth(aes(x=x,y=y),method=lm) + xlab("Measured") + ylab("Predicted") + ggtitle(paste("Comparison of ",parameter," measured vs predicted for module ",module)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
dev.off()
# Establish the correlation between predicted and modeled
# This is the data to report with the figure (R2, CI, signif, etc.)
cor.test(df$x,df$y)

## Identify node centrality based on co-occurence data for each OTU in the module

TOM = TOMsimilarityFromExpr(otu_WGCNA2, power = 12, corType = "pearson");
# Select submodule of interest based on high correlation and signficance
module<-"darkred"; # <- this changes with module color being currently explored
# Select module probes
probes = names(otu_WGCNA)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
write.csv(modTOM, paste("nodeconnections_MM_B_",parameter,".csv"))

# Number of cells above 0.25 threshhold <-- this number is flexible and should change with your data
x<- as.data.frame(rowSums(modTOM > 0.1))
write.csv(x, paste("nodes_MM_B_",parameter,".csv"))


# make scatter hive plots
# You will need to make the Nodeworksheet by combining OTUinfo table for submodule of interest, VIP socres, and Nodeconnections. See workflow for more details.

hive_MM_B_temp<- read.csv("Temp_darkred_MM_B_nodeworksheet.csv", header=T) #(Figure 5-8)
hive_MM_B_temp$OTU<-factor(hive_MM_B_temp$OTU, levels = unique(hive_MM_B_temp$OTU))
MM_B_temp_plot <- ggplot(hive_MM_B_temp, aes(x= hive_MM_B_temp$connectivity, y= hive_MM_B_temp$GS.Temp)) +
  geom_point(aes(size = hive_MM_B_temp$VIP, colour = hive_MM_B_temp$Phylum, alpha = 0.5)) +
  scale_size_area(max_size= 10) +
  scale_color_brewer(type = qual, palette = "Set1", direction = -1) +
  theme_bw() +
  scale_alpha(guide=FALSE) +
  geom_text_repel(aes(label = hive_MM_B_temp$names), force = 3, size = 3) +
  labs(x="Node centrality", y="Correlation to Temp", color="Phylum",
       size = "VIP")
MM_B_temp_plot

MM_B_legend <- g_legend(MM_B_plot)
grid.arrange(MM_B_plot+ theme(legend.position = 'none'),MM_B_temp_plot+ theme(legend.position = 'none'),MM_B_legend, nrow=1)

##################
### WGCNA SF_B ###
##################

NUT <- read.table("NUT_SF_B.txt", header = TRUE, row.names = 1, sep ="\t")
OTU <- read.table("OTU_SF_B.txt", header = TRUE, row.names = 1, sep="\t")
TAX <- read.table("diel_tax_singdoubrem_cols.txt", header = TRUE, row.names = 1, sep = "\t")

rownames(OTU) <- paste0("OTU", 1:nrow(OTU))
rownames(OTU)

TAX <- as.matrix(TAX, rownames.force = NA)
rownames(TAX) <- paste0("OTU", 1:nrow(TAX))
rownames(TAX)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)

physeq = phyloseq(OTU,TAX)

META = sample_data(NUT)
rownames(META) <- sample_names(physeq)

META = sample_data(META)

# get all the data in a phyloseq instance, or whatever
ALL = phyloseq(OTU,TAX,META)
ALL

### Prepping data to remove rare or erronious OTUs #####

# We need to de-noise the data by plotting the number of reads on a curve and look for the inflection point

at.least.n.in.m <- function(x, n, m){
  all(x[x>0]>=n)&length(x[x>0])>=m
}
counts<- rep(0,10)
for (i in 1:length(counts)){
  rows.to.keep<- apply(otu_table(ALL, taxa_are_rows = TRUE), 1,at.least.n.in.m, n=i, m=2)
  counts[i]<-sum(rows.to.keep)
}

plot(1:10, counts, xlab= 'Min sequences in 2 samples', ylab= 'Number of taxa remaining')

### Inflection point was 2 ###

# Filter taxa that arent seen more than twice in greater than 25% of the data.
CUT2<-filter_taxa(ALL, function(x) sum(x > 2) > (0.25*length(x)), TRUE)
CUT2
write.csv(tax_table(CUT2), "SF_B_seqnames.csv")
write.csv(otu_table(CUT2), "SF_B_seqOTU.csv")

### Deseq2 Normalization
### convert the counts to integer mode and group data to normalize together(e.g. region)
dds_WGCNA <- phyloseq_to_deseq2(CUT2, ~ Diel)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans_WGCNA = apply(counts(dds_WGCNA), 1, gm_mean)
dds_WGCNA = estimateSizeFactors(dds_WGCNA, geoMeans=geoMeans_WGCNA)
dds_WGCNA = estimateDispersions(dds_WGCNA)

### variance stabilizing transformation
#Make it a OTU again for phyloseq (no tax or sample data)
vst_WGCNA <- varianceStabilizingTransformation(dds_WGCNA, blind=FALSE)
vstMat_WGCNA <- assay(vst_WGCNA)
vstMat_WGCNA[vstMat_WGCNA<0]<-0
vst.otu.WGCNA <- otu_table(vstMat_WGCNA, taxa_are_rows=TRUE)
write.csv(otu_table(vst.otu.WGCNA), "vst.otu.WGCNA.SF_B.csv")

#### BEGIN WGCNA #####

#### Maybe try this next time... ####
vst <- read.csv("vst.otu.WGCNA.SF_B.csv")
dim(vst)
names(vst)
otu_WGCNA <- as.data.frame(t(vst))
names(otu_WGCNA) = vst$X
rownames(otu_WGCNA) = names(vst)
dim(otu_WGCNA)
names(otu_WGCNA)


NMDS_META <- read.table("NUT_SF_B.txt", header = T, row.names = 1)


# Because we are using RStudio, we have to disable threads
# As consequence of this, maybe it would be better to do this step in regular ol' R
disableWGCNAThreads()
options(stringsAsFactors = FALSE)
## Identify beta to ensure scale free topology
powers = c(seq(4,10,by=1), seq(12,20, by=2))

pst <- pickSoftThreshold(otu_WGCNA, powerVector=powers, blockSize = 3528, verbose=2)

# Plot the results of soft thresholds if you wish:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(pst$fitIndices[,1], pst$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(pst$fitIndices[,1], pst$fitIndices[,5], labels=powers, cex=cex1,col="red")

#check the powerEstimate to choose a threshold
pst

# Create adjacency matrix by raising OTU matrix by beta and identify subnetworks (modules)
otu_WGCNA2 <- as.matrix(otu_WGCNA[2:13,1:3528])
mode(otu_WGCNA2)
class(otu_WGCNA2) <- "numeric"
mode(otu_WGCNA2)

# Check that the network ensures scale-free topology at that power
# R should be close to 1 (R > 0.8, I believe), should see a straight line.
##### scaleFreePlot #####
# here we define the adjacency matrix using soft thresholding with beta=12
ADJ1=abs(cor(otu_WGCNA2,use="p"))^9
# When you have relatively few genes (<5000) use the following code
#k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=otu_WGCNA2,power=9)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
scaleFreeFitIndex(k)

#R^2 of 0.88, this suggests we meet the assumption of scale-free topol.

# power of 12 chosen based on powerEstimate from 'pst'
net = blockwiseModules(otu_WGCNA2, power=9, minModuleSize=30, maxBlockSize = 3528,
                       corType = "pearson", saveTOMs = TRUE, 
                       saveTOMFileBase = "blockwiseTOM", pamStage=FALSE, verbose=5)
# Plot the dendrogram
moduleLabels = net$colors
moduleColors = net$colors
MEs = net$MEs
geneTree = net$dendrograms[[1]]
pdf("plotDendro_SF_B.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
# Save data

# Identify Eigenvalues for subnetworks by sample
nPop = ncol(otu_WGCNA2)
nSamples = nrow(otu_WGCNA2)
MEsO = moduleEigengenes(otu_WGCNA2, moduleColors)$eigengenes
MEs = orderMEs(MEsO)
save(MEs, moduleLabels, moduleColors, geneTree,file = "Module-networkConstruction-auto.RData")
# Save data
write.csv(file="Module_eigen_values_SF_B.csv",MEs)
write.csv(file="Module_composition_SF_B.csv",net$colors)

##Correlate Eigenvalues to metadata and create heatmap

moduleTraitCor = cor(MEs, NMDS_META[,c(5:14)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), " (",signif(moduleTraitPvalue, 1), ")"
                   , sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf("Correlation_SF_B.pdf",width=12,height=8)
par(mar = c(12, 12, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(NMDS_META[,c(5:14)]),yLabels = names(MEs)
               ,ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50)
               ,textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,
               cex.lab = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()

### !!! No strong relationships between a module and DO, pH, or Temp !!! ###

##################
### WGCNA PM_B ###
##################

NUT <- read.table("NUT_PM_B.txt", header = TRUE, row.names = 1, sep ="\t")
OTU <- read.table("OTU_PM_B.txt", header = TRUE, row.names = 1, sep="\t")
TAX <- read.table("diel_tax_singdoubrem_cols.txt", header = TRUE, row.names = 1, sep = "\t")

rownames(OTU) <- paste0("OTU", 1:nrow(OTU))
rownames(OTU)

TAX <- as.matrix(TAX, rownames.force = NA)
rownames(TAX) <- paste0("OTU", 1:nrow(TAX))
rownames(TAX)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)

physeq = phyloseq(OTU,TAX)

META = sample_data(NUT)
rownames(META) <- sample_names(physeq)

META = sample_data(META)

# get all the data in a phyloseq instance, or whatever
ALL = phyloseq(OTU,TAX,META)
ALL

### Prepping data to remove rare or erronious OTUs #####

# We need to de-noise the data by plotting the number of reads on a curve and look for the inflection point

at.least.n.in.m <- function(x, n, m){
  all(x[x>0]>=n)&length(x[x>0])>=m
}
counts<- rep(0,10)
for (i in 1:length(counts)){
  rows.to.keep<- apply(otu_table(ALL, taxa_are_rows = TRUE), 1,at.least.n.in.m, n=i, m=2)
  counts[i]<-sum(rows.to.keep)
}

plot(1:10, counts, xlab= 'Min sequences in 2 samples', ylab= 'Number of taxa remaining')

### Inflection point was 3 ###

# Filter taxa that arent seen more than twice in greater than 25% of the data.
CUT2<-filter_taxa(ALL, function(x) sum(x > 3) > (0.25*length(x)), TRUE)
CUT2
write.csv(tax_table(CUT2), "PM_B_seqnames.csv")
write.csv(otu_table(CUT2), "PM_B_seqOTU.csv")

### Deseq2 Normalization
### convert the counts to integer mode and group data to normalize together(e.g. region)
dds_WGCNA <- phyloseq_to_deseq2(CUT2, ~ Diel)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans_WGCNA = apply(counts(dds_WGCNA), 1, gm_mean)
dds_WGCNA = estimateSizeFactors(dds_WGCNA, geoMeans=geoMeans_WGCNA)
dds_WGCNA = estimateDispersions(dds_WGCNA)

### variance stabilizing transformation
#Make it a OTU again for phyloseq (no tax or sample data)
vst_WGCNA <- varianceStabilizingTransformation(dds_WGCNA, blind=FALSE)
vstMat_WGCNA <- assay(vst_WGCNA)
vstMat_WGCNA[vstMat_WGCNA<0]<-0
vst.otu.WGCNA <- otu_table(vstMat_WGCNA, taxa_are_rows=TRUE)
write.csv(otu_table(vst.otu.WGCNA), "vst.otu.WGCNA.PM_B.csv")

#### BEGIN WGCNA #####

#### Maybe try this next time... ####
vst <- read.csv("vst.otu.WGCNA.PM_B.csv")
dim(vst)
names(vst)
otu_WGCNA <- as.data.frame(t(vst))
names(otu_WGCNA) = vst$X
rownames(otu_WGCNA) = names(vst)
dim(otu_WGCNA)
names(otu_WGCNA)


NMDS_META <- read.table("NUT_PM_B.txt", header = T, row.names = 1)


# Because we are using RStudio, we have to disable threads
# As consequence of this, maybe it would be better to do this step in regular ol' R
disableWGCNAThreads()
options(stringsAsFactors = FALSE)
## Identify beta to ensure scale free topology
powers = c(seq(4,10,by=1), seq(12,20, by=2))

pst <- pickSoftThreshold(otu_WGCNA, powerVector=powers, blockSize = 5123, verbose=2)

# Plot the results of soft thresholds if you wish:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(pst$fitIndices[,1], pst$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(pst$fitIndices[,1], pst$fitIndices[,5], labels=powers, cex=cex1,col="red")

#check the powerEstimate to choose a threshold
pst

# Create adjacency matrix by raising OTU matrix by beta and identify subnetworks (modules)
otu_WGCNA2 <- as.matrix(otu_WGCNA[2:12,1:5123])
mode(otu_WGCNA2)
class(otu_WGCNA2) <- "numeric"
mode(otu_WGCNA2)

# Check that the network ensures scale-free topology at that power
# R should be close to 1 (R > 0.8, I believe), should see a straight line.
##### scaleFreePlot #####
# here we define the adjacency matrix using soft thresholding with beta=12
ADJ1=abs(cor(otu_WGCNA2,use="p"))^12
# When you have relatively few genes (<5000) use the following code
#k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=otu_WGCNA2,power=12)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
scaleFreeFitIndex(k)

#R^2 of 0.84, this suggests we meet the assumption of scale-free topol.

# power of 12 chosen based on powerEstimate from 'pst'
net = blockwiseModules(otu_WGCNA2, power=12, minModuleSize=30, maxBlockSize = 5123,
                       corType = "pearson", saveTOMs = TRUE, 
                       saveTOMFileBase = "blockwiseTOM", pamStage=FALSE, verbose=5)
# Plot the dendrogram
moduleLabels = net$colors
moduleColors = net$colors
MEs = net$MEs
geneTree = net$dendrograms[[1]]
pdf("plotDendro_PM_B.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
# Save data

# Identify Eigenvalues for subnetworks by sample
nPop = ncol(otu_WGCNA2)
nSamples = nrow(otu_WGCNA2)
MEsO = moduleEigengenes(otu_WGCNA2, moduleColors)$eigengenes
MEs = orderMEs(MEsO)
save(MEs, moduleLabels, moduleColors, geneTree,file = "Module-networkConstruction-auto.RData")
# Save data
write.csv(file="Module_eigen_values_PM_B.csv",MEs)
write.csv(file="Module_composition_PM_B.csv",net$colors)

##Correlate Eigenvalues to metadata and create heatmap

moduleTraitCor = cor(MEs, NMDS_META[,c(5:14)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), " (",signif(moduleTraitPvalue, 1), ")"
                   , sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf("Correlation_PM_B.pdf",width=12,height=8)
par(mar = c(12, 12, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(NMDS_META[,c(5:14)]),yLabels = names(MEs)
               ,ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50)
               ,textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,
               cex.lab = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()

## Now make a plot for specific module <-> trait (metadata component) pairings
# This allows us to explore the structure of submodule OTU correlations with a given metadata component
# Here we will use "env variable" as trait and "color" as module
# First get the links between all modules and this trait

###DO-lightcyan
parameter<-"DO"
weight <- as.data.frame(NMDS_META[,parameter])
names(weight) = "DO"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(otu_WGCNA2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(otu_WGCNA2, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
# Then look at the specific module of interest
module<-"lightcyan"
column = match(module, modNames);
moduleGenes = moduleColors==module
pdf(paste(module,"-vs-",parameter,".pdf"))
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = paste("Population significance for ",parameter),main = paste("Module membership vs. population significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

##Incorporate OTU taxonomy info to identify specific OTUs within a module, coalesce correlative data by OTU for later use in network analysis
annot2 = read.csv("PM_B_seqnames.csv", colClasses = "character", header = TRUE)
colnames(annot2)[1] <- "OTU"
dim(annot2)
names(annot2)
probes2 = names(otu_WGCNA)
probes2annot2 = match(probes2, annot2$OTU)
# The following is the number or probes without annotation:
sum(is.na(probes2annot2))
# Should return 0.
# Create the starting data frame
geneInfo0 = data.frame(otu_orig = probes2,
                       OTU = annot2$OTU[probes2annot2],
                       Domain = annot2$Domain[probes2annot2],
                       Phylum = annot2$Phylum[probes2annot2],
                       Class = annot2$Class[probes2annot2],
                       Order = annot2$Order[probes2annot2],
                       Family = annot2$Family[probes2annot2],
                       Genus = annot2$Genus[probes2annot2],
                       names = annot2$OTU[probes2annot2],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance according to metadata component
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("PM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.DO)); # <- change geneInfo$'X'
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = paste("OTUInfo_PM_B_",module))

# Get the similarity between eigenvalues and weight
MET = orderMEs(cbind(MEs, weight))
pdf("Adjacency_trait_ME_PM_B_DO.pdf") # <- This changes with env. variable
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene dendrogram with ",parameter), marDendro = c(0,4,2,0),plotHeatmaps = FALSE)
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene adjacency heatmap with ",parameter), marHeatmap = c(3,4,2,2),plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

################################################################################
#           Moving to machine learning PLS and VIP scores                      #
################################################################################

library("pls")
source("VIP.R")
# metadata component (e.g. S, designated "weight" here as above) is the same as before, we just replicate the row names for pls
#the below r2 threshold changes with user input, good correlation above "0.x"
#this is important, because VIP scores don't really mean anything without a good correlation
th_r2<-0.3
subnetwork<-otu_WGCNA[2:12,moduleGenes]
row.names(weight)<-row.names(subnetwork)
weight <- as.matrix(weight)
subnetwork <- as.matrix(subnetwork)
class(weight) <- "numeric"
class(subnetwork) <- "numeric"
pls_result<-plsr(weight ~ subnetwork, validation="LOO",method="oscorespls")

r2_vector<-R2(pls_result)
max<-0
max_comp<--1
for (j in 1:length(r2_vector$val)){
  if(r2_vector$val[j]>th_r2){         # We will only look at the PLS if the correlation is better than 0.3
    if(r2_vector$val[j]>max){
      max<-r2_vector$val[j]
      max_comp<-r2_vector$comp[j]
    }
  }
}
print(paste(" the max r2 is ",max," corresponding to comp ",max_comp,sep="",".pdf"))

if(max==0){
  print ("No good correlation, we stop here")
} else{
  print("Good correlation, we check the VIP!")
}
# Checking the VIP
output<-paste("VIP_values_with_",parameter,sep="")
vip_result<-VIP(pls_result)
vip_components<-sort(vip_result[max_comp,],decreasing=TRUE)[1:68] # <- this value changes with subnetwork, check dim of 'subnetwork' file
for (i in 1:68){ # <- this value changes as the line above
  cat(paste("Rank ",i," we have ",names(vip_components[i])," with a VIP of ",vip_components[i],"\n",sep=""),file=output,append=TRUE)
}
weight_2 <- as.data.frame(weight[!is.na(weight)])
df<-data.frame(x=weight_2,y=pls_result$validation$pred[,,max_comp])
colnames(df)<-c("x","y")
pdf(paste("measured_vs_predicted_PM_B_",module,"-vs-",parameter,".pdf"))
ggplot(data=df) + geom_point(aes(x=x,y=y)) + geom_smooth(aes(x=x,y=y),method=lm) + xlab("Measured") + ylab("Predicted") + ggtitle(paste("Comparison of ",parameter," measured vs predicted for module ",module)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
dev.off()
# Establish the correlation between predicted and modeled
# This is the data to report with the figure (R2, CI, signif, etc.)
cor.test(df$x,df$y)

## Identify node centrality based on co-occurence data for each OTU in the module

TOM = TOMsimilarityFromExpr(otu_WGCNA2, power = 12, corType = "pearson");
# Select submodule of interest based on high correlation and signficance
module<-"lightcyan"; # <- this changes with module color being currently explored
# Select module probes
probes = names(otu_WGCNA)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
write.csv(modTOM, paste("nodeconnections_PM_B_",parameter,".csv"))

# Number of cells above 0.25 threshhold <-- this number is flexible and should change with your data
x<- as.data.frame(rowSums(modTOM > 0.1))
write.csv(x, paste("nodes_PM_B_",parameter,".csv"))


# make scatter hive plots
# You will need to make the Nodeworksheet by combining OTUinfo table for submodule of interest, VIP socres, and Nodeconnections. See workflow for more details.

hive_PM_B<- read.csv("DO_lightcyan_PM_B_nodeworksheet.csv", header=T) #(Figure 5-8)
hive_PM_B$OTU<-factor(hive_PM_B$OTU, levels = unique(hive_PM_B$OTU))
PM_B_plot <- ggplot(hive_PM_B, aes(x= hive_PM_B$connectivity, y= hive_PM_B$GS.DO)) +
  geom_point(aes(size = hive_PM_B$VIP, colour = hive_PM_B$Phylum, alpha = 0.5)) +
  scale_size_area(max_size= 10) +
  scale_color_brewer(type = qual, palette = "Set3", direction = -1) +
  theme_bw() +
  scale_alpha(guide=FALSE) +
  geom_text_repel(aes(label = hive_PM_B$names), force = 3, size = 3) +
  labs(x="Node centrality", y="Correlation to DO", color="Phylum",
       size = "VIP")
PM_B_plot

## The same subnetwork is most associated with pH. We'll try this now.

###pH-lightcyan
parameter<-"pH"
weight <- as.data.frame(NMDS_META[,parameter])
names(weight) = "pH"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(otu_WGCNA2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(otu_WGCNA2, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
# Then look at the specific module of interest
module<-"lightcyan"
column = match(module, modNames);
moduleGenes = moduleColors==module
pdf(paste(module,"-vs-",parameter,".pdf"))
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = paste("Population significance for ",parameter),main = paste("Module membership vs. population significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

##Incorporate OTU taxonomy info to identify specific OTUs within a module, coalesce correlative data by OTU for later use in network analysis
annot2 = read.csv("PM_B_seqnames.csv", colClasses = "character", header = TRUE)
colnames(annot2)[1] <- "OTU"
dim(annot2)
names(annot2)
probes2 = names(otu_WGCNA)
probes2annot2 = match(probes2, annot2$OTU)
# The following is the number or probes without annotation:
sum(is.na(probes2annot2))
# Should return 0.
# Create the starting data frame
geneInfo0 = data.frame(otu_orig = probes2,
                       OTU = annot2$OTU[probes2annot2],
                       Domain = annot2$Domain[probes2annot2],
                       Phylum = annot2$Phylum[probes2annot2],
                       Class = annot2$Class[probes2annot2],
                       Order = annot2$Order[probes2annot2],
                       Family = annot2$Family[probes2annot2],
                       Genus = annot2$Genus[probes2annot2],
                       names = annot2$OTU[probes2annot2],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance according to metadata component
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("PM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.pH)); # <- change geneInfo$'X'
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = paste("OTUInfo_PM_B_",module))

# Get the similarity between eigenvalues and weight
MET = orderMEs(cbind(MEs, weight))
pdf("Adjacency_trait_ME_PM_B_pH.pdf") # <- This changes with env. variable
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene dendrogram with ",parameter), marDendro = c(0,4,2,0),plotHeatmaps = FALSE)
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene adjacency heatmap with ",parameter), marHeatmap = c(3,4,2,2),plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

################################################################################
#           Moving to machine learning PLS and VIP scores                      #
################################################################################

library("pls")
source("VIP.R")
# metadata component (e.g. S, designated "weight" here as above) is the same as before, we just replicate the row names for pls
#the below r2 threshold changes with user input, good correlation above "0.x"
#this is important, because VIP scores don't really mean anything without a good correlation
th_r2<-0.3
subnetwork<-otu_WGCNA[2:12,moduleGenes]
row.names(weight)<-row.names(subnetwork)
weight <- as.matrix(weight)
subnetwork <- as.matrix(subnetwork)
class(weight) <- "numeric"
class(subnetwork) <- "numeric"
pls_result<-plsr(weight ~ subnetwork, validation="LOO",method="oscorespls")

r2_vector<-R2(pls_result)
max<-0
max_comp<--1
for (j in 1:length(r2_vector$val)){
  if(r2_vector$val[j]>th_r2){         # We will only look at the PLS if the correlation is better than 0.3
    if(r2_vector$val[j]>max){
      max<-r2_vector$val[j]
      max_comp<-r2_vector$comp[j]
    }
  }
}
print(paste(" the max r2 is ",max," corresponding to comp ",max_comp,sep="",".pdf"))

if(max==0){
  print ("No good correlation, we stop here")
} else{
  print("Good correlation, we check the VIP!")
}
# Checking the VIP
output<-paste("VIP_values_with_",parameter,sep="")
vip_result<-VIP(pls_result)
vip_components<-sort(vip_result[max_comp,],decreasing=TRUE)[1:68] # <- this value changes with subnetwork, check dim of 'subnetwork' file
for (i in 1:68){ # <- this value changes as the line above
  cat(paste("Rank ",i," we have ",names(vip_components[i])," with a VIP of ",vip_components[i],"\n",sep=""),file=output,append=TRUE)
}
weight_2 <- as.data.frame(weight[!is.na(weight)])
df<-data.frame(x=weight_2,y=pls_result$validation$pred[,,max_comp])
colnames(df)<-c("x","y")
pdf(paste("measured_vs_predicted_PM_B_",module,"-vs-",parameter,".pdf"))
ggplot(data=df) + geom_point(aes(x=x,y=y)) + geom_smooth(aes(x=x,y=y),method=lm) + xlab("Measured") + ylab("Predicted") + ggtitle(paste("Comparison of ",parameter," measured vs predicted for module ",module)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
dev.off()
# Establish the correlation between predicted and modeled
# This is the data to report with the figure (R2, CI, signif, etc.)
cor.test(df$x,df$y)

## Identify node centrality based on co-occurence data for each OTU in the module

TOM = TOMsimilarityFromExpr(otu_WGCNA2, power = 12, corType = "pearson");
# Select submodule of interest based on high correlation and signficance
module<-"lightcyan"; # <- this changes with module color being currently explored
# Select module probes
probes = names(otu_WGCNA)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
write.csv(modTOM, paste("nodeconnections_PM_B_",parameter,".csv"))

# Number of cells above 0.25 threshhold <-- this number is flexible and should change with your data
x<- as.data.frame(rowSums(modTOM > 0.1))
write.csv(x, paste("nodes_PM_B_",parameter,".csv"))


# make scatter hive plots
# You will need to make the Nodeworksheet by combining OTUinfo table for submodule of interest, VIP socres, and Nodeconnections. See workflow for more details.

hive_PM_B_pH<- read.csv("pH_lightcyan_PM_B_nodeworksheet.csv", header=T) #(Figure 5-8)
hive_PM_B_pH$OTU<-factor(hive_PM_B_pH$OTU, levels = unique(hive_PM_B_pH$OTU))
PM_B_pH_plot <- ggplot(hive_PM_B_pH, aes(x= hive_PM_B_pH$connectivity, y= hive_PM_B_pH$GS.pH)) +
  geom_point(aes(size = hive_PM_B_pH$VIP, colour = hive_PM_B_pH$Phylum, alpha = 0.5)) +
  scale_size_area(max_size= 10) +
  scale_color_brewer(type = qual, palette = "Set3", direction = -1) +
  theme_bw() +
  scale_alpha(guide=FALSE) +
  geom_text_repel(aes(label = hive_PM_B_pH$names), force = 3, size = 3) +
  labs(x="Node centrality", y="Correlation to pH", color="Phylum",
       size = "VIP")
PM_B_pH_plot

PM_B_legend <- g_legend(PM_B_plot)
grid.arrange(PM_B_plot+ theme(legend.position = 'none'),PM_B_pH_plot+ theme(legend.position = 'none'),PM_B_legend, nrow=1)

##################
### WGCNA CW_B ###
##################

NUT <- read.table("NUT_CW_B.txt", header = TRUE, row.names = 1, sep ="\t")
OTU <- read.table("OTU_CW_B.txt", header = TRUE, row.names = 1, sep="\t")
TAX <- read.table("diel_tax_singdoubrem_cols.txt", header = TRUE, row.names = 1, sep = "\t")

rownames(OTU) <- paste0("OTU", 1:nrow(OTU))
rownames(OTU)

TAX <- as.matrix(TAX, rownames.force = NA)
rownames(TAX) <- paste0("OTU", 1:nrow(TAX))
rownames(TAX)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)

physeq = phyloseq(OTU,TAX)

META = sample_data(NUT)
rownames(META) <- sample_names(physeq)

META = sample_data(META)

# get all the data in a phyloseq instance, or whatever
ALL = phyloseq(OTU,TAX,META)
ALL

### Prepping data to remove rare or erronious OTUs #####

# We need to de-noise the data by plotting the number of reads on a curve and look for the inflection point

at.least.n.in.m <- function(x, n, m){
  all(x[x>0]>=n)&length(x[x>0])>=m
}
counts<- rep(0,10)
for (i in 1:length(counts)){
  rows.to.keep<- apply(otu_table(ALL, taxa_are_rows = TRUE), 1,at.least.n.in.m, n=i, m=2)
  counts[i]<-sum(rows.to.keep)
}

plot(1:10, counts, xlab= 'Min sequences in 2 samples', ylab= 'Number of taxa remaining')

### Inflection point was 2 ###

# Filter taxa that arent seen more than twice in greater than 25% of the data.
CUT2<-filter_taxa(ALL, function(x) sum(x > 2) > (0.25*length(x)), TRUE)
CUT2
write.csv(tax_table(CUT2), "CW_B_seqnames.csv")
write.csv(otu_table(CUT2), "CW_B_seqOTU.csv")

### Deseq2 Normalization
### convert the counts to integer mode and group data to normalize together(e.g. region)
dds_WGCNA <- phyloseq_to_deseq2(CUT2, ~ Diel)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans_WGCNA = apply(counts(dds_WGCNA), 1, gm_mean)
dds_WGCNA = estimateSizeFactors(dds_WGCNA, geoMeans=geoMeans_WGCNA)
dds_WGCNA = estimateDispersions(dds_WGCNA)

### variance stabilizing transformation
#Make it a OTU again for phyloseq (no tax or sample data)
vst_WGCNA <- varianceStabilizingTransformation(dds_WGCNA, blind=FALSE)
vstMat_WGCNA <- assay(vst_WGCNA)
vstMat_WGCNA[vstMat_WGCNA<0]<-0
vst.otu.WGCNA <- otu_table(vstMat_WGCNA, taxa_are_rows=TRUE)
write.csv(otu_table(vst.otu.WGCNA), "vst.otu.WGCNA.CW_B.csv")

#### BEGIN WGCNA #####

#### Maybe try this next time... ####
vst <- read.csv("vst.otu.WGCNA.CW_B.csv")
dim(vst)
names(vst)
otu_WGCNA <- as.data.frame(t(vst))
names(otu_WGCNA) = vst$X
rownames(otu_WGCNA) = names(vst)
dim(otu_WGCNA)
names(otu_WGCNA)


NMDS_META <- read.table("NUT_CW_B.txt", header = T, row.names = 1)


# Because we are using RStudio, we have to disable threads
# As consequence of this, maybe it would be better to do this step in regular ol' R
disableWGCNAThreads()
options(stringsAsFactors = FALSE)
## Identify beta to ensure scale free topology
powers = c(seq(4,10,by=1), seq(12,20, by=2))

pst <- pickSoftThreshold(otu_WGCNA, powerVector=powers, blockSize = 3604, verbose=2)

# Plot the results of soft thresholds if you wish:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(pst$fitIndices[,1], pst$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(pst$fitIndices[,1], pst$fitIndices[,5], labels=powers, cex=cex1,col="red")

#check the powerEstimate to choose a threshold
pst

# Create adjacency matrix by raising OTU matrix by beta and identify subnetworks (modules)
otu_WGCNA2 <- as.matrix(otu_WGCNA[2:13,1:3604])
mode(otu_WGCNA2)
class(otu_WGCNA2) <- "numeric"
mode(otu_WGCNA2)

# Check that the network ensures scale-free topology at that power
# R should be close to 1 (R > 0.8, I believe), should see a straight line.
##### scaleFreePlot #####
# here we define the adjacency matrix using soft thresholding with beta=12
ADJ1=abs(cor(otu_WGCNA2,use="p"))^6
# When you have relatively few genes (<5000) use the following code
k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
#k=softConnectivity(datE=otu_WGCNA2,power=6)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
scaleFreeFitIndex(k)

#R^2 of 0.86, this suggests we meet the assumption of scale-free topol.

# power of 6 chosen based on powerEstimate from 'pst'
net = blockwiseModules(otu_WGCNA2, power=6, minModuleSize=30, maxBlockSize = 3604,
                       corType = "pearson", saveTOMs = TRUE, 
                       saveTOMFileBase = "blockwiseTOM", pamStage=FALSE, verbose=5)
# Plot the dendrogram
moduleLabels = net$colors
moduleColors = net$colors
MEs = net$MEs
geneTree = net$dendrograms[[1]]
pdf("plotDendro_cW_B.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
# Save data

# Identify Eigenvalues for subnetworks by sample
nPop = ncol(otu_WGCNA2)
nSamples = nrow(otu_WGCNA2)
MEsO = moduleEigengenes(otu_WGCNA2, moduleColors)$eigengenes
MEs = orderMEs(MEsO)
save(MEs, moduleLabels, moduleColors, geneTree,file = "Module-networkConstruction-auto.RData")
# Save data
write.csv(file="Module_eigen_values_CW_B.csv",MEs)
write.csv(file="Module_composition_CW_B.csv",net$colors)

##Correlate Eigenvalues to metadata and create heatmap

moduleTraitCor = cor(MEs, NMDS_META[,c(5:14)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), " (",signif(moduleTraitPvalue, 1), ")"
                   , sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf("Correlation_CW_B.pdf",width=12,height=8)
par(mar = c(12, 12, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(NMDS_META[,c(5:14)]),yLabels = names(MEs)
               ,ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50)
               ,textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,
               cex.lab = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()


### !!! No strong relationship between a module and DO !!! ###
### However, strong relationship with temperature ###

## Now make a plot for specific module <-> trait (metadata component) pairings
# This allows us to explore the structure of submodule OTU correlations with a given metadata component
# Here we will use "env variable" as trait and "color" as module
# First get the links between all modules and this trait

###Temp-tan
parameter<-"Temp"
weight <- as.data.frame(NMDS_META[,parameter])
names(weight) = "Temp"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(otu_WGCNA2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(otu_WGCNA2, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
# Then look at the specific module of interest
module<-"tan"
column = match(module, modNames);
moduleGenes = moduleColors==module
pdf(paste(module,"-vs-",parameter,".pdf"))
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = paste("Population significance for ",parameter),main = paste("Module membership vs. population significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

##Incorporate OTU taxonomy info to identify specific OTUs within a module, coalesce correlative data by OTU for later use in network analysis
annot2 = read.csv("CW_B_seqnames.csv", colClasses = "character", header = TRUE)
colnames(annot2)[1] <- "OTU"
dim(annot2)
names(annot2)
probes2 = names(otu_WGCNA)
probes2annot2 = match(probes2, annot2$OTU)
# The following is the number or probes without annotation:
sum(is.na(probes2annot2))
# Should return 0.
# Create the starting data frame
geneInfo0 = data.frame(otu_orig = probes2,
                       OTU = annot2$OTU[probes2annot2],
                       Domain = annot2$Domain[probes2annot2],
                       Phylum = annot2$Phylum[probes2annot2],
                       Class = annot2$Class[probes2annot2],
                       Order = annot2$Order[probes2annot2],
                       Family = annot2$Family[probes2annot2],
                       Genus = annot2$Genus[probes2annot2],
                       names = annot2$OTU[probes2annot2],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance according to metadata component
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Temp)); # <- change geneInfo$'X'
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = paste("OTUInfo_CW_B_",module))

# Get the similarity between eigenvalues and weight
MET = orderMEs(cbind(MEs, weight))
pdf("Adjacency_trait_ME_CW_B_Temp.pdf") # <- This changes with env. variable
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene dendrogram with ",parameter), marDendro = c(0,4,2,0),plotHeatmaps = FALSE)
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene adjacency heatmap with ",parameter), marHeatmap = c(3,4,2,2),plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

################################################################################
#           Moving to machine learning PLS and VIP scores                      #
################################################################################

library("pls")
source("VIP.R")
# metadata component (e.g. S, designated "weight" here as above) is the same as before, we just replicate the row names for pls
#the below r2 threshold changes with user input, good correlation above "0.x"
#this is important, because VIP scores don't really mean anything without a good correlation
th_r2<-0.3
subnetwork<-otu_WGCNA[2:13,moduleGenes]
row.names(weight)<-row.names(subnetwork)
weight <- as.matrix(weight)
subnetwork <- as.matrix(subnetwork)
class(weight) <- "numeric"
class(subnetwork) <- "numeric"
pls_result<-plsr(weight ~ subnetwork, validation="LOO",method="oscorespls")

r2_vector<-R2(pls_result)
max<-0
max_comp<--1
for (j in 1:length(r2_vector$val)){
  if(r2_vector$val[j]>th_r2){         # We will only look at the PLS if the correlation is better than 0.3
    if(r2_vector$val[j]>max){
      max<-r2_vector$val[j]
      max_comp<-r2_vector$comp[j]
    }
  }
}
print(paste(" the max r2 is ",max," corresponding to comp ",max_comp,sep="",".pdf"))

if(max==0){
  print ("No good correlation, we stop here")
} else{
  print("Good correlation, we check the VIP!")
}
# Checking the VIP
output<-paste("VIP_values_with_",parameter,sep="")
vip_result<-VIP(pls_result)
vip_components<-sort(vip_result[max_comp,],decreasing=TRUE)[1:72] # <- this value changes with subnetwork, check dim of 'subnetwork' file
for (i in 1:72){ # <- this value changes as the line above
  cat(paste("Rank ",i," we have ",names(vip_components[i])," with a VIP of ",vip_components[i],"\n",sep=""),file=output,append=TRUE)
}
weight_2 <- as.data.frame(weight[!is.na(weight)])
df<-data.frame(x=weight_2,y=pls_result$validation$pred[,,max_comp])
colnames(df)<-c("x","y")
pdf(paste("measured_vs_predicted_CW_B_",module,"-vs-",parameter,".pdf"))
ggplot(data=df) + geom_point(aes(x=x,y=y)) + geom_smooth(aes(x=x,y=y),method=lm) + xlab("Measured") + ylab("Predicted") + ggtitle(paste("Comparison of ",parameter," measured vs predicted for module ",module)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
dev.off()
# Establish the correlation between predicted and modeled
# This is the data to report with the figure (R2, CI, signif, etc.)
cor.test(df$x,df$y)

## Identify node centrality based on co-occurence data for each OTU in the module

TOM = TOMsimilarityFromExpr(otu_WGCNA2, power = 6, corType = "pearson");
# Select submodule of interest based on high correlation and signficance
module<-"tan"; # <- this changes with module color being currently explored
# Select module probes
probes = names(otu_WGCNA)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
write.csv(modTOM, paste("nodeconnections_CW_B_",parameter,".csv"))

# Number of cells above 0.25 threshhold <-- this number is flexible and should change with your data
x<- as.data.frame(rowSums(modTOM > 0.1))
write.csv(x, paste("nodes_CW_B_",parameter,".csv"))


# make scatter hive plots
# You will need to make the Nodeworksheet by combining OTUinfo table for submodule of interest, VIP socres, and Nodeconnections. See workflow for more details.

hive_CW_B<- read.csv("Temp_tan_CW_B_nodeworksheet.csv", header=T) #(Figure 5-8)
hive_CW_B$OTU<-factor(hive_CW_B$OTU, levels = unique(hive_CW_B$OTU))
CW_B_plot <- ggplot(hive_CW_B, aes(x= hive_CW_B$connectivity, y= hive_CW_B$GS.Temp)) +
  geom_point(aes(size = hive_CW_B$VIP, colour = hive_CW_B$Phylum, alpha = 0.5)) +
  scale_size_area(max_size= 10) +
  scale_color_brewer(type = qual, palette = "Set1", direction = -1) +
  theme_bw() +
  scale_alpha(guide=FALSE) +
  geom_text_repel(aes(label = hive_CW_B$names), force = 3, size = 4) +
  labs(x="Node centrality", y="Correlation to Temp", color="Phylum",
       size = "VIP")
CW_B_plot


##################
### WGCNA MM_A ###
##################

NUT <- read.table("NUT_MM_A.txt", header = TRUE, row.names = 1, sep ="\t")
OTU <- read.table("OTU_MM_A.txt", header = TRUE, row.names = 1, sep="\t")
TAX <- read.table("diel_tax_singdoubrem_cols.txt", header = TRUE, row.names = 1, sep = "\t")

rownames(OTU) <- paste0("OTU", 1:nrow(OTU))
rownames(OTU)

TAX <- as.matrix(TAX, rownames.force = NA)
rownames(TAX) <- paste0("OTU", 1:nrow(TAX))
rownames(TAX)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)

physeq = phyloseq(OTU,TAX)

META = sample_data(NUT)
rownames(META) <- sample_names(physeq)

META = sample_data(META)

# get all the data in a phyloseq instance, or whatever
ALL = phyloseq(OTU,TAX,META)
ALL

### Prepping data to remove rare or erronious OTUs #####

# We need to de-noise the data by plotting the number of reads on a curve and look for the inflection point

at.least.n.in.m <- function(x, n, m){
  all(x[x>0]>=n)&length(x[x>0])>=m
}
counts<- rep(0,10)
for (i in 1:length(counts)){
  rows.to.keep<- apply(otu_table(ALL, taxa_are_rows = TRUE), 1,at.least.n.in.m, n=i, m=2)
  counts[i]<-sum(rows.to.keep)
}

plot(1:10, counts, xlab= 'Min sequences in 2 samples', ylab= 'Number of taxa remaining')

### Inflection point was 2 ###

# Filter taxa that arent seen more than twice in greater than 25% of the data.
CUT2<-filter_taxa(ALL, function(x) sum(x > 2) > (0.25*length(x)), TRUE)
CUT2
write.csv(tax_table(CUT2), "MM_A_seqnames.csv")
write.csv(otu_table(CUT2), "MM_A_seqOTU.csv")

### Deseq2 Normalization
### convert the counts to integer mode and group data to normalize together(e.g. region)
dds_WGCNA <- phyloseq_to_deseq2(CUT2, ~ Diel)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans_WGCNA = apply(counts(dds_WGCNA), 1, gm_mean)
dds_WGCNA = estimateSizeFactors(dds_WGCNA, geoMeans=geoMeans_WGCNA)
dds_WGCNA = estimateDispersions(dds_WGCNA)

### variance stabilizing transformation
#Make it a OTU again for phyloseq (no tax or sample data)
vst_WGCNA <- varianceStabilizingTransformation(dds_WGCNA, blind=FALSE)
vstMat_WGCNA <- assay(vst_WGCNA)
vstMat_WGCNA[vstMat_WGCNA<0]<-0
vst.otu.WGCNA <- otu_table(vstMat_WGCNA, taxa_are_rows=TRUE)
write.csv(otu_table(vst.otu.WGCNA), "vst.otu.WGCNA.MM_A.csv")

#### BEGIN WGCNA #####

#### Maybe try this next time... ####
vst <- read.csv("vst.otu.WGCNA.MM_A.csv")
dim(vst)
names(vst)
otu_WGCNA <- as.data.frame(t(vst))
names(otu_WGCNA) = vst$X
rownames(otu_WGCNA) = names(vst)
dim(otu_WGCNA)
names(otu_WGCNA)


NMDS_META <- read.table("NUT_MM_A.txt", header = T, row.names = 1)


# Because we are using RStudio, we have to disable threads
# As consequence of this, maybe it would be better to do this step in regular ol' R
disableWGCNAThreads()
options(stringsAsFactors = FALSE)
## Identify beta to ensure scale free topology
powers = c(seq(4,10,by=1), seq(12,20, by=2))

pst <- pickSoftThreshold(otu_WGCNA, powerVector=powers, blockSize = 3708, verbose=2)

# Plot the results of soft thresholds if you wish:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(pst$fitIndices[,1], pst$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(pst$fitIndices[,1], pst$fitIndices[,5], labels=powers, cex=cex1,col="red")

#check the powerEstimate to choose a threshold
pst

# Create adjacency matrix by raising OTU matrix by beta and identify subnetworks (modules)
otu_WGCNA2 <- as.matrix(otu_WGCNA[2:13,1:3708])
mode(otu_WGCNA2)
class(otu_WGCNA2) <- "numeric"
mode(otu_WGCNA2)

# Check that the network ensures scale-free topology at that power
# R should be close to 1 (R > 0.8, I believe), should see a straight line.
##### scaleFreePlot #####
# here we define the adjacency matrix using soft thresholding with beta=12
ADJ1=abs(cor(otu_WGCNA2,use="p"))^12
# When you have relatively few genes (<5000) use the following code
#k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=otu_WGCNA2,power=12)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
scaleFreeFitIndex(k)

#R^2 of 0.87, this suggests we meet the assumption of scale-free topol.

# power of 12 chosen based on powerEstimate from 'pst'
net = blockwiseModules(otu_WGCNA2, power=12, minModuleSize=30, maxBlockSize = 3767,
                       corType = "pearson", saveTOMs = TRUE, 
                       saveTOMFileBase = "blockwiseTOM", pamStage=FALSE, verbose=5)
# Plot the dendrogram
moduleLabels = net$colors
moduleColors = net$colors
MEs = net$MEs
geneTree = net$dendrograms[[1]]
pdf("plotDendro_MM_A.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
# Save data

# Identify Eigenvalues for subnetworks by sample
nPop = ncol(otu_WGCNA2)
nSamples = nrow(otu_WGCNA2)
MEsO = moduleEigengenes(otu_WGCNA2, moduleColors)$eigengenes
MEs = orderMEs(MEsO)
save(MEs, moduleLabels, moduleColors, geneTree,file = "Module-networkConstruction-auto.RData")
# Save data
write.csv(file="Module_eigen_values_MM_A.csv",MEs)
write.csv(file="Module_composition_MM_A.csv",net$colors)

##Correlate Eigenvalues to metadata and create heatmap

moduleTraitCor = cor(MEs, NMDS_META[,c(5:14)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), " (",signif(moduleTraitPvalue, 1), ")"
                   , sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf("Correlation_MM_A.pdf",width=12,height=8)
par(mar = c(12, 12, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(NMDS_META[,c(5:14)]),yLabels = names(MEs)
               ,ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50)
               ,textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,
               cex.lab = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()

### !!! no strong relationship between a module and DO !!! ###


##################
### WGCNA SF_A ###
##################

NUT <- read.table("NUT_SF_A.txt", header = TRUE, row.names = 1, sep ="\t")
OTU <- read.table("OTU_SF_A.txt", header = TRUE, row.names = 1, sep="\t")
TAX <- read.table("diel_tax_singdoubrem_cols.txt", header = TRUE, row.names = 1, sep = "\t")

rownames(OTU) <- paste0("OTU", 1:nrow(OTU))
rownames(OTU)

TAX <- as.matrix(TAX, rownames.force = NA)
rownames(TAX) <- paste0("OTU", 1:nrow(TAX))
rownames(TAX)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)

physeq = phyloseq(OTU,TAX)

META = sample_data(NUT)
rownames(META) <- sample_names(physeq)

META = sample_data(META)

# get all the data in a phyloseq instance, or whatever
ALL = phyloseq(OTU,TAX,META)
ALL

### Prepping data to remove rare or erronious OTUs #####

# We need to de-noise the data by plotting the number of reads on a curve and look for the inflection point

at.least.n.in.m <- function(x, n, m){
  all(x[x>0]>=n)&length(x[x>0])>=m
}
counts<- rep(0,10)
for (i in 1:length(counts)){
  rows.to.keep<- apply(otu_table(ALL, taxa_are_rows = TRUE), 1,at.least.n.in.m, n=i, m=2)
  counts[i]<-sum(rows.to.keep)
}

plot(1:10, counts, xlab= 'Min sequences in 2 samples', ylab= 'Number of taxa remaining')

### Inflection point was 2 ###

# Filter taxa that arent seen more than twice in greater than 25% of the data.
CUT2<-filter_taxa(ALL, function(x) sum(x > 2) > (0.25*length(x)), TRUE)
CUT2
write.csv(tax_table(CUT2), "SF_A_seqnames.csv")
write.csv(otu_table(CUT2), "SF_A_seqOTU.csv")

### Deseq2 Normalization
### convert the counts to integer mode and group data to normalize together(e.g. region)
dds_WGCNA <- phyloseq_to_deseq2(CUT2, ~ Diel)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans_WGCNA = apply(counts(dds_WGCNA), 1, gm_mean)
dds_WGCNA = estimateSizeFactors(dds_WGCNA, geoMeans=geoMeans_WGCNA)
dds_WGCNA = estimateDispersions(dds_WGCNA)

### variance stabilizing transformation
#Make it a OTU again for phyloseq (no tax or sample data)
vst_WGCNA <- varianceStabilizingTransformation(dds_WGCNA, blind=FALSE)
vstMat_WGCNA <- assay(vst_WGCNA)
vstMat_WGCNA[vstMat_WGCNA<0]<-0
vst.otu.WGCNA <- otu_table(vstMat_WGCNA, taxa_are_rows=TRUE)
write.csv(otu_table(vst.otu.WGCNA), "vst.otu.WGCNA.SF_A.csv")

#### BEGIN WGCNA #####

#### Maybe try this next time... ####
vst <- read.csv("vst.otu.WGCNA.SF_A.csv")
dim(vst)
names(vst)
otu_WGCNA <- as.data.frame(t(vst))
names(otu_WGCNA) = vst$X
rownames(otu_WGCNA) = names(vst)
dim(otu_WGCNA)
names(otu_WGCNA)


NMDS_META <- read.table("NUT_SF_A.txt", header = T, row.names = 1)


# Because we are using RStudio, we have to disable threads
# As consequence of this, maybe it would be better to do this step in regular ol' R
disableWGCNAThreads()
options(stringsAsFactors = FALSE)
## Identify beta to ensure scale free topology
powers = c(seq(4,10,by=1), seq(12,20, by=2))

pst <- pickSoftThreshold(otu_WGCNA, powerVector=powers, blockSize = 3511, verbose=2)

# Plot the results of soft thresholds if you wish:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(pst$fitIndices[,1], pst$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(pst$fitIndices[,1], pst$fitIndices[,5], labels=powers, cex=cex1,col="red")

#check the powerEstimate to choose a threshold
pst

# Create adjacency matrix by raising OTU matrix by beta and identify subnetworks (modules)
otu_WGCNA2 <- as.matrix(otu_WGCNA[2:13,1:3511])
mode(otu_WGCNA2)
class(otu_WGCNA2) <- "numeric"
mode(otu_WGCNA2)

# Check that the network ensures scale-free topology at that power
# R should be close to 1 (R > 0.8, I believe), should see a straight line.
##### scaleFreePlot #####
# here we define the adjacency matrix using soft thresholding with beta=12
ADJ1=abs(cor(otu_WGCNA2,use="p"))^6
# When you have relatively few genes (<5000) use the following code
#k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=otu_WGCNA2,power=6)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
scaleFreeFitIndex(k)

#R^2 of 0.84, this suggests we meet the assumption of scale-free topol.

# power of 6 chosen based on powerEstimate from 'pst'
net = blockwiseModules(otu_WGCNA2, power=6, minModuleSize=30, maxBlockSize = 3511,
                       corType = "pearson", saveTOMs = TRUE, 
                       saveTOMFileBase = "blockwiseTOM", pamStage=FALSE, verbose=5)
# Plot the dendrogram
moduleLabels = net$colors
moduleColors = net$colors
MEs = net$MEs
geneTree = net$dendrograms[[1]]
pdf("plotDendro_SF_A.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
# Save data

# Identify Eigenvalues for subnetworks by sample
nPop = ncol(otu_WGCNA2)
nSamples = nrow(otu_WGCNA2)
MEsO = moduleEigengenes(otu_WGCNA2, moduleColors)$eigengenes
MEs = orderMEs(MEsO)
save(MEs, moduleLabels, moduleColors, geneTree,file = "Module-networkConstruction-auto.RData")
# Save data
write.csv(file="Module_eigen_values_SF_A.csv",MEs)
write.csv(file="Module_composition_SF_A.csv",net$colors)

##Correlate Eigenvalues to metadata and create heatmap

moduleTraitCor = cor(MEs, NMDS_META[,c(5:14)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), " (",signif(moduleTraitPvalue, 1), ")"
                   , sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf("Correlation_SF_A.pdf",width=12,height=8)
par(mar = c(12, 12, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(NMDS_META[,c(5:14)]),yLabels = names(MEs)
               ,ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50)
               ,textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,
               cex.lab = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()

### !!! No strong relationship between a module and DO !!! ###
### However, strong relationship between a module and pH ###

## Now make a plot for specific module <-> trait (metadata component) pairings
# This allows us to explore the structure of submodule OTU correlations with a given metadata component
# Here we will use "env variable" as trait and "color" as module
# First get the links between all modules and this trait

###pH-darkslateblue
parameter<-"pH"
weight <- as.data.frame(NMDS_META[,parameter])
names(weight) = "pH"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(otu_WGCNA2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(otu_WGCNA2, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
# Then look at the specific module of interest
module<-"darkslateblue"
column = match(module, modNames);
moduleGenes = moduleColors==module
pdf(paste(module,"-vs-",parameter,".pdf"))
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = paste("Population significance for ",parameter),main = paste("Module membership vs. population significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

##Incorporate OTU taxonomy info to identify specific OTUs within a module, coalesce correlative data by OTU for later use in network analysis
annot2 = read.csv("SF_A_seqnames.csv", colClasses = "character", header = TRUE)
colnames(annot2)[1] <- "OTU"
dim(annot2)
names(annot2)
probes2 = names(otu_WGCNA)
probes2annot2 = match(probes2, annot2$OTU)
# The following is the number or probes without annotation:
sum(is.na(probes2annot2))
# Should return 0.
# Create the starting data frame
geneInfo0 = data.frame(otu_orig = probes2,
                       OTU = annot2$OTU[probes2annot2],
                       Domain = annot2$Domain[probes2annot2],
                       Phylum = annot2$Phylum[probes2annot2],
                       Class = annot2$Class[probes2annot2],
                       Order = annot2$Order[probes2annot2],
                       Family = annot2$Family[probes2annot2],
                       Genus = annot2$Genus[probes2annot2],
                       names = annot2$OTU[probes2annot2],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance according to metadata component
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.pH)); # <- change geneInfo$'X'
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = paste("OTUInfo_SF_A_",module))

# Get the similarity between eigenvalues and weight
MET = orderMEs(cbind(MEs, weight))
pdf("Adjacency_trait_ME_SF_A_pH.pdf") # <- This changes with env. variable
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene dendrogram with ",parameter), marDendro = c(0,4,2,0),plotHeatmaps = FALSE)
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene adjacency heatmap with ",parameter), marHeatmap = c(3,4,2,2),plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

################################################################################
#           Moving to machine learning PLS and VIP scores                      #
################################################################################

library("pls")
source("VIP.R")
# metadata component (e.g. S, designated "weight" here as above) is the same as before, we just replicate the row names for pls
#the below r2 threshold changes with user input, good correlation above "0.x"
#this is important, because VIP scores don't really mean anything without a good correlation
th_r2<-0.3
subnetwork<-otu_WGCNA[2:13,moduleGenes]
row.names(weight)<-row.names(subnetwork)
weight <- as.matrix(weight)
subnetwork <- as.matrix(subnetwork)
class(weight) <- "numeric"
class(subnetwork) <- "numeric"
pls_result<-plsr(weight ~ subnetwork, validation="LOO",method="oscorespls")

r2_vector<-R2(pls_result)
max<-0
max_comp<--1
for (j in 1:length(r2_vector$val)){
  if(r2_vector$val[j]>th_r2){         # We will only look at the PLS if the correlation is better than 0.3
    if(r2_vector$val[j]>max){
      max<-r2_vector$val[j]
      max_comp<-r2_vector$comp[j]
    }
  }
}
print(paste(" the max r2 is ",max," corresponding to comp ",max_comp,sep="",".pdf"))

if(max==0){
  print ("No good correlation, we stop here")
} else{
  print("Good correlation, we check the VIP!")
}
# Checking the VIP
output<-paste("VIP_values_with_",parameter,sep="")
vip_result<-VIP(pls_result)
vip_components<-sort(vip_result[max_comp,],decreasing=TRUE)[1:30] # <- this value changes with subnetwork, check dim of 'subnetwork' file
for (i in 1:30){ # <- this value changes as the line above
  cat(paste("Rank ",i," we have ",names(vip_components[i])," with a VIP of ",vip_components[i],"\n",sep=""),file=output,append=TRUE)
}
weight_2 <- as.data.frame(weight[!is.na(weight)])
df<-data.frame(x=weight_2,y=pls_result$validation$pred[,,max_comp])
colnames(df)<-c("x","y")
pdf(paste("measured_vs_predicted_SF_A_",module,"-vs-",parameter,".pdf"))
ggplot(data=df) + geom_point(aes(x=x,y=y)) + geom_smooth(aes(x=x,y=y),method=lm) + xlab("Measured") + ylab("Predicted") + ggtitle(paste("Comparison of ",parameter," measured vs predicted for module ",module)) + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
dev.off()
# Establish the correlation between predicted and modeled
# This is the data to report with the figure (R2, CI, signif, etc.)
cor.test(df$x,df$y)

## Identify node centrality based on co-occurence data for each OTU in the module

TOM = TOMsimilarityFromExpr(otu_WGCNA2, power = 6, corType = "pearson");
# Select submodule of interest based on high correlation and signficance
module<-"darkslateblue"; # <- this changes with module color being currently explored
# Select module probes
probes = names(otu_WGCNA)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
write.csv(modTOM, paste("nodeconnections_SF_A_",parameter,".csv"))

# Number of cells above 0.25 threshhold <-- this number is flexible and should change with your data
x<- as.data.frame(rowSums(modTOM > 0.15))
write.csv(x, paste("nodes_SF_A_",parameter,".csv"))


# make scatter hive plots
# You will need to make the Nodeworksheet by combining OTUinfo table for submodule of interest, VIP socres, and Nodeconnections. See workflow for more details.

hive_SF_A<- read.csv("pH_darkslateblue_SF_A_nodeworksheet.csv", header=T) #(Figure 5-8)
hive_SF_A$OTU<-factor(hive_SF_A$OTU, levels = unique(hive_SF_A$OTU))
SF_A_plot <- ggplot(hive_SF_A, aes(x= hive_SF_A$connectivity, y= hive_SF_A$GS.pH)) +
  geom_point(aes(size = hive_SF_A$VIP, colour = hive_SF_A$Phylum, alpha = 0.5)) +
  scale_size_area(max_size= 10) +
  scale_color_brewer(type = qual, palette = "Set1", direction = -1) +
  theme_bw() +
  scale_alpha(guide=FALSE) +
  geom_text_repel(aes(label = hive_SF_A$names), force = 3, size = 4) +
  labs(x="Node centrality", y="Correlation to pH", color="Phylum",
       size = "VIP")
SF_A_plot


##################
### WGCNA CW_A ###
##################

NUT <- read.table("NUT_CW_A.txt", header = TRUE, row.names = 1, sep ="\t")
OTU <- read.table("OTU_CW_A.txt", header = TRUE, row.names = 1, sep="\t")
TAX <- read.table("diel_tax_singdoubrem_cols.txt", header = TRUE, row.names = 1, sep = "\t")

rownames(OTU) <- paste0("OTU", 1:nrow(OTU))
rownames(OTU)

TAX <- as.matrix(TAX, rownames.force = NA)
rownames(TAX) <- paste0("OTU", 1:nrow(TAX))
rownames(TAX)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)

physeq = phyloseq(OTU,TAX)

META = sample_data(NUT)
rownames(META) <- sample_names(physeq)

META = sample_data(META)

# get all the data in a phyloseq instance, or whatever
ALL = phyloseq(OTU,TAX,META)
ALL

### Prepping data to remove rare or erronious OTUs #####

# We need to de-noise the data by plotting the number of reads on a curve and look for the inflection point

at.least.n.in.m <- function(x, n, m){
  all(x[x>0]>=n)&length(x[x>0])>=m
}
counts<- rep(0,10)
for (i in 1:length(counts)){
  rows.to.keep<- apply(otu_table(ALL, taxa_are_rows = TRUE), 1,at.least.n.in.m, n=i, m=2)
  counts[i]<-sum(rows.to.keep)
}

plot(1:10, counts, xlab= 'Min sequences in 2 samples', ylab= 'Number of taxa remaining')

### Inflection point was 2 ###

# Filter taxa that arent seen more than twice in greater than 25% of the data.
CUT2<-filter_taxa(ALL, function(x) sum(x > 2) > (0.25*length(x)), TRUE)
CUT2
write.csv(tax_table(CUT2), "CW_A_seqnames.csv")
write.csv(otu_table(CUT2), "CW_A_seqOTU.csv")

### Deseq2 Normalization
### convert the counts to integer mode and group data to normalize together(e.g. region)
dds_WGCNA <- phyloseq_to_deseq2(CUT2, ~ Diel)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans_WGCNA = apply(counts(dds_WGCNA), 1, gm_mean)
dds_WGCNA = estimateSizeFactors(dds_WGCNA, geoMeans=geoMeans_WGCNA)
dds_WGCNA = estimateDispersions(dds_WGCNA)

### variance stabilizing transformation
#Make it a OTU again for phyloseq (no tax or sample data)
vst_WGCNA <- varianceStabilizingTransformation(dds_WGCNA, blind=FALSE)
vstMat_WGCNA <- assay(vst_WGCNA)
vstMat_WGCNA[vstMat_WGCNA<0]<-0
vst.otu.WGCNA <- otu_table(vstMat_WGCNA, taxa_are_rows=TRUE)
write.csv(otu_table(vst.otu.WGCNA), "vst.otu.WGCNA.CW_A.csv")

#### BEGIN WGCNA #####

#### Maybe try this next time... ####
vst <- read.csv("vst.otu.WGCNA.CW_A.csv")
dim(vst)
names(vst)
otu_WGCNA <- as.data.frame(t(vst))
names(otu_WGCNA) = vst$X
rownames(otu_WGCNA) = names(vst)
dim(otu_WGCNA)
names(otu_WGCNA)


NMDS_META <- read.table("NUT_CW_A.txt", header = T, row.names = 1)


# Because we are using RStudio, we have to disable threads
# As consequence of this, maybe it would be better to do this step in regular ol' R
disableWGCNAThreads()
options(stringsAsFactors = FALSE)
## Identify beta to ensure scale free topology
powers = c(seq(4,10,by=1), seq(12,20, by=2))

pst <- pickSoftThreshold(otu_WGCNA, powerVector=powers, blockSize = 4659, verbose=2)

# Plot the results of soft thresholds if you wish:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(pst$fitIndices[,1], pst$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(pst$fitIndices[,1], pst$fitIndices[,5], labels=powers, cex=cex1,col="red")

#check the powerEstimate to choose a threshold
pst

# Create adjacency matrix by raising OTU matrix by beta and identify subnetworks (modules)
otu_WGCNA2 <- as.matrix(otu_WGCNA[2:13,1:4659])
mode(otu_WGCNA2)
class(otu_WGCNA2) <- "numeric"
mode(otu_WGCNA2)

# Check that the network ensures scale-free topology at that power
# R should be close to 1 (R > 0.8, I believe), should see a straight line.
##### scaleFreePlot #####
# here we define the adjacency matrix using soft thresholding with beta=12
ADJ1=abs(cor(otu_WGCNA2,use="p"))^18
# When you have relatively few genes (<5000) use the following code
#k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=otu_WGCNA2,power=18)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
scaleFreeFitIndex(k)

#R^2 of 0.84, this suggests we meet the assumption of scale-free topol.

# power of 18 chosen based on powerEstimate from 'pst'
net = blockwiseModules(otu_WGCNA2, power=18, minModuleSize=30, maxBlockSize = 4659,
                       corType = "pearson", saveTOMs = TRUE, 
                       saveTOMFileBase = "blockwiseTOM", pamStage=FALSE, verbose=5)
# Plot the dendrogram
moduleLabels = net$colors
moduleColors = net$colors
MEs = net$MEs
geneTree = net$dendrograms[[1]]
pdf("plotDendro_CW_A.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
# Save data

# Identify Eigenvalues for subnetworks by sample
nPop = ncol(otu_WGCNA2)
nSamples = nrow(otu_WGCNA2)
MEsO = moduleEigengenes(otu_WGCNA2, moduleColors)$eigengenes
MEs = orderMEs(MEsO)
save(MEs, moduleLabels, moduleColors, geneTree,file = "Module-networkConstruction-auto.RData")
# Save data
write.csv(file="Module_eigen_values_CW_A.csv",MEs)
write.csv(file="Module_composition_CW_A.csv",net$colors)

##Correlate Eigenvalues to metadata and create heatmap

moduleTraitCor = cor(MEs, NMDS_META[,c(5:14)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), " (",signif(moduleTraitPvalue, 1), ")"
                   , sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf("Correlation_CW_A.pdf",width=12,height=8)
par(mar = c(12, 12, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(NMDS_META[,c(5:14)]),yLabels = names(MEs)
               ,ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50)
               ,textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,
               cex.lab = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()

### !!! No strong relationship between a module and DO !!! ###


##################
### WGCNA PM_A ###
##################

NUT <- read.table("NUT_PM_A.txt", header = TRUE, row.names = 1, sep ="\t")
OTU <- read.table("OTU_PM_A.txt", header = TRUE, row.names = 1, sep="\t")
TAX <- read.table("diel_tax_singdoubrem_cols.txt", header = TRUE, row.names = 1, sep = "\t")

rownames(OTU) <- paste0("OTU", 1:nrow(OTU))
rownames(OTU)

TAX <- as.matrix(TAX, rownames.force = NA)
rownames(TAX) <- paste0("OTU", 1:nrow(TAX))
rownames(TAX)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)

physeq = phyloseq(OTU,TAX)

META = sample_data(NUT)
rownames(META) <- sample_names(physeq)

META = sample_data(META)

# get all the data in a phyloseq instance, or whatever
ALL = phyloseq(OTU,TAX,META)
ALL

### Prepping data to remove rare or erronious OTUs #####

# We need to de-noise the data by plotting the number of reads on a curve and look for the inflection point

at.least.n.in.m <- function(x, n, m){
  all(x[x>0]>=n)&length(x[x>0])>=m
}
counts<- rep(0,10)
for (i in 1:length(counts)){
  rows.to.keep<- apply(otu_table(ALL, taxa_are_rows = TRUE), 1,at.least.n.in.m, n=i, m=2)
  counts[i]<-sum(rows.to.keep)
}

plot(1:10, counts, xlab= 'Min sequences in 2 samples', ylab= 'Number of taxa remaining')

### Inflection point was 2 ###

# Filter taxa that arent seen more than twice in greater than 25% of the data.
CUT2<-filter_taxa(ALL, function(x) sum(x > 2) > (0.25*length(x)), TRUE)
CUT2
write.csv(tax_table(CUT2), "PM_A_seqnames.csv")
write.csv(otu_table(CUT2), "PM_A_seqOTU.csv")

### Deseq2 Normalization
### convert the counts to integer mode and group data to normalize together(e.g. region)
dds_WGCNA <- phyloseq_to_deseq2(CUT2, ~ Diel)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans_WGCNA = apply(counts(dds_WGCNA), 1, gm_mean)
dds_WGCNA = estimateSizeFactors(dds_WGCNA, geoMeans=geoMeans_WGCNA)
dds_WGCNA = estimateDispersions(dds_WGCNA)

### variance stabilizing transformation
#Make it a OTU again for phyloseq (no tax or sample data)
vst_WGCNA <- varianceStabilizingTransformation(dds_WGCNA, blind=FALSE)
vstMat_WGCNA <- assay(vst_WGCNA)
vstMat_WGCNA[vstMat_WGCNA<0]<-0
vst.otu.WGCNA <- otu_table(vstMat_WGCNA, taxa_are_rows=TRUE)
write.csv(otu_table(vst.otu.WGCNA), "vst.otu.WGCNA.PM_A.csv")

#### BEGIN WGCNA #####

#### Maybe try this next time... ####
vst <- read.csv("vst.otu.WGCNA.PM_A.csv")
dim(vst)
names(vst)
otu_WGCNA <- as.data.frame(t(vst))
names(otu_WGCNA) = vst$X
rownames(otu_WGCNA) = names(vst)
dim(otu_WGCNA)
names(otu_WGCNA)


NMDS_META <- read.table("NUT_PM_A.txt", header = T, row.names = 1)


# Because we are using RStudio, we have to disable threads
# As consequence of this, maybe it would be better to do this step in regular ol' R
disableWGCNAThreads()
options(stringsAsFactors = FALSE)
## Identify beta to ensure scale free topology
powers = c(seq(4,10,by=1), seq(12,20, by=2))

pst <- pickSoftThreshold(otu_WGCNA, powerVector=powers, blockSize = 5022, verbose=2)

# Plot the results of soft thresholds if you wish:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(pst$fitIndices[,1], -sign(pst$fitIndices[,3])*pst$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(pst$fitIndices[,1], pst$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(pst$fitIndices[,1], pst$fitIndices[,5], labels=powers, cex=cex1,col="red")

#check the powerEstimate to choose a threshold
pst

# Create adjacency matrix by raising OTU matrix by beta and identify subnetworks (modules)
otu_WGCNA2 <- as.matrix(otu_WGCNA[2:13,1:5022])
mode(otu_WGCNA2)
class(otu_WGCNA2) <- "numeric"
mode(otu_WGCNA2)

# Check that the network ensures scale-free topology at that power
# R should be close to 1 (R > 0.8, I believe), should see a straight line.
##### scaleFreePlot #####
# here we define the adjacency matrix using soft thresholding with beta=12
ADJ1=abs(cor(otu_WGCNA2,use="p"))^6
# When you have relatively few genes (<5000) use the following code
#k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=otu_WGCNA2,power=6)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
scaleFreeFitIndex(k)

#R^2 of 0.83, this suggests we meet the assumption of scale-free topol.

# power of 18 chosen based on powerEstimate from 'pst'
net = blockwiseModules(otu_WGCNA2, power=6, minModuleSize=30, maxBlockSize = 5022,
                       corType = "pearson", saveTOMs = TRUE, 
                       saveTOMFileBase = "blockwiseTOM", pamStage=FALSE, verbose=5)
# Plot the dendrogram
moduleLabels = net$colors
moduleColors = net$colors
MEs = net$MEs
geneTree = net$dendrograms[[1]]
pdf("plotDendro_PM_A.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
# Save data

# Identify Eigenvalues for subnetworks by sample
nPop = ncol(otu_WGCNA2)
nSamples = nrow(otu_WGCNA2)
MEsO = moduleEigengenes(otu_WGCNA2, moduleColors)$eigengenes
MEs = orderMEs(MEsO)
save(MEs, moduleLabels, moduleColors, geneTree,file = "Module-networkConstruction-auto.RData")
# Save data
write.csv(file="Module_eigen_values_PM_A.csv",MEs)
write.csv(file="Module_composition_PM_A.csv",net$colors)

##Correlate Eigenvalues to metadata and create heatmap

moduleTraitCor = cor(MEs, NMDS_META[,c(5:14)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), " (",signif(moduleTraitPvalue, 1), ")"
                   , sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf("Correlation_PM_A.pdf",width=12,height=8)
par(mar = c(12, 12, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(NMDS_META[,c(5:14)]),yLabels = names(MEs)
               ,ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50)
               ,textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,
               cex.lab = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()

### No strong relationship between a module and DO !!! ###

# Double-plot for WGCNA figures
require(gridExtra)
grid.arrange(CW_B_plot,MM_B_plot,PM_B_plot,SF_A_plot,ncol = 2)
multiplot(CW_B_plot,MM_B_plot,PM_B_plot,SF_A_plot,cols = 2)
