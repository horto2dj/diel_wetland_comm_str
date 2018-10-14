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
library(Rmisc)

memory.size(max = TRUE)
memory.limit(size = 15000)
memory.limit()

setwd("~/Grad School/Dissertation/Diel/2015 DNA study/stats")


################
### Diel PCA ###
################

#chem data
chem <- read.table("chem_phys_dat_allsites_clean.txt",header=TRUE, row.names = 1, sep="\t")
dim(chem)

chem <- chem[,c(1:11,13:19)]

# change pH to [H+]
pH <- 10^-(chem[,8])

chem <- cbind(chem[,1:7],pH,chem[,9:18])


#pearson correlation analysis
#check to see if any variables are autocorrelated
rcorr(as.matrix(chem[,5:18]), type="pearson")

#running correlation on log-transformed data introduced zeros and -infs, would not run

#List (>0.7 and sig at 0.01 level)
#TDS and spCond correlate, represented by spCond
#NPOC and OrP correlate, represented by OrP
#TN and Ammonia correlate, represented by TN
#SRP and TP correlate, represented by TP

chem <- chem[,c(1:6,8:13,16:17)]

#apply PCA
chem.pca <- prcomp(chem[,5:14], center = TRUE, scale. = TRUE)

#print method
print(chem.pca)

#plot method
plot(chem.pca, type = "lines")

#summary method
summary(chem.pca)

#create PCA plot in ggplot

wetland <- as.factor(chem$Wetland)
diel <- as.factor(chem$Diel)
GLPCA <- ggbiplot(chem.pca, obs.scale = 1, var.scale = 1, 
                  groups = chem$Wetland, ellipse = FALSE, 
                  circle = FALSE, varname.adjust = 1) +
  #scale_color_hue(l = 0) +
  geom_point(aes(fill = wetland, shape = diel),size=3) +
  scale_shape_manual(name="Time",
                     labels=c("Dawn", "Dusk"),
                     values = c(21,22)) +
  scale_fill_hue(l = 40, name = "Wetland") +
  theme(legend.direction = 'vertical', 
        legend.position = 'right') +
  #geom_text_repel(aes(label = wetland)) +
  theme_classic() +
  guides(shape=guide_legend(title="Time"),
         fill=guide_legend(override.aes = list(shape = 21)),
         color=guide_legend(title="Wetland")) +
  stat_ellipse(aes(fill= wetland), show.legend = FALSE,
             geom="polygon", level=0.95, alpha=0.2)
GLPCA

#Finding point coordinates on biplot
names(chem.pca)
chem.pca$x

# MANOVA
PCA_manova <- manova(chem.pca$x ~ Wetland*Point*Diel, data = chem)
summary(PCA_manova)

PCA_perMANOVA <- adonis(chem.pca$x ~ Wetland*Point*Diel,
                        data = chem, method = "euclidean")
PCA_perMANOVA

pairwise.adonis.euc <- function(x,factors, sim.function = 'vegdist', sim.method = 'euclidean', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])],method = "euclidean" );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 *** 0.001 ** 0.01 * 0.05 . 0.1  1")
  return(pairw.res)
  
} 


chem_mat <- as.matrix(chem.pca$x)
pairwise.adonis.euc(chem_mat, factors = chem$Wetland)

chem_dist <- vegdist(chem.pca$x, method = "euclidean")
bd_chem <- betadisper(chem_dist, group = chem$Wetland)
plot(bd_chem)
summary(bd_chem)
anova(bd_chem)

# beta-dispersion is not significant, which means perMANOVA can be safely executed.

### Chemical Line Graphs ###

line.dat <- read.table("chem_phys_dat_allsites_clean_lines.txt",header=TRUE, row.names = 1, sep="\t")
head(line.dat)
line.dat <- line.dat[,c(1:6,8,12:13)]
head(line.dat)

#Get STDevs by groups - DO
ldstat <- summarySE(line.dat, measurevar = "DO", groupvars = c("Wetland","Time"))
ldstat

#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LDO <- ggplot(ldstat, aes(x=Time, y=DO, colour=Wetland, group = Wetland)) + 
  ylab("Temperature (°C)") +
  geom_errorbar(aes(ymin=DO-se, ymax=DO+se), colour = "black", width=.1, 
                position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Wetland", breaks = c("CW", "MM", "PM", "SF"),
                  labels = c("CW", "MM", "PM", "SF"), l = 40) +
  scale_y_continuous(breaks = 0:30*2) +
  theme_bw()

LDO

#Get STDevs by groups - Temp
ldstat_Temp <- summarySE(line.dat, measurevar = "Temp", groupvars = c("Wetland","Time"))
ldstat_Temp

#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LTemp <- ggplot(ldstat_Temp, aes(x=Time, y=Temp, colour=Wetland, group = Wetland)) + 
  ylab("Temperature (°C)") +
  geom_errorbar(aes(ymin=Temp-se, ymax=Temp+se), colour = "black", width=.1, 
                position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Wetland", breaks = c("CW", "MM", "PM", "SF"),
                  labels = c("CW", "MM", "PM", "SF"), l = 40) +
  scale_y_continuous(breaks = 0:30*2) +
  theme_bw()

#Get STDevs by groups - pH
ldstat_pH <- summarySE(line.dat, measurevar = "pH", groupvars = c("Wetland","Time"))
ldstat_pH

#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LpH <- ggplot(ldstat_pH, aes(x=Time, y=pH, colour=Wetland, group = Wetland)) + 
  ylab("pH") +
  geom_errorbar(aes(ymin=pH-se, ymax=pH+se), colour = "black", width=.1, 
                position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Wetland", breaks = c("CW", "MM", "PM", "SF"),
                  labels = c("CW", "MM", "PM", "SF"), l = 40) +
  #scale_color_brewer(type = "qual", palette = "Set1", name = "Wetland", l = 50) +
  #scale_y_continuous(breaks = 0:30*2) +
  theme_bw()

LpH
#Get STDevs by groups - Nitrate
ldstat_Nitrate <- summarySE(line.dat, measurevar = "Nitrate", groupvars = c("Wetland","Time"))
ldstat_Nitrate

#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LNitrate <- ggplot(ldstat_Nitrate, aes(x=Time, y=Nitrate, colour=Wetland, group = Wetland)) + 
  ylab("NO3 (mg/L)") +
  geom_errorbar(aes(ymin=Nitrate-se, ymax=Nitrate+se), colour = "black", width=.1, 
                position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Wetland", breaks = c("CW", "MM", "PM", "SF"),
                  labels = c("CW", "MM", "PM", "SF"), l = 40) +
  #scale_y_continuous(breaks = 0:30*2) +
  theme_bw()

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

LineCombo <- grid.arrange(arrangeGrob(LDO,LTemp,LpH,LNitrate,nrow=2))


### DO Line Graphs for all wetlands ###

line.dat.DO <- read.table("chem_phys_dat_allsites_clean_lines3.txt",header=TRUE, row.names = 1, sep="\t")
head(line.dat.DO)
line.dat.DO <- line.dat.DO[,c(2:5,12,16)]
head(line.dat.DO)

#Get STDevs by groups - DO
ldostat_CW <- summarySE(line.dat.DO[1:24,], measurevar = "DO", groupvars = c("Point","Time","Year"))
ldostat_CW


#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LDO_CW <- ggplot(ldostat_CW, aes(x=Time, y=DO, colour=Point, group = Point)) + 
  ylab("DO (mg/L)") +
  labs(x = "Time") +
  geom_errorbar(aes(ymin=DO-se, ymax=DO+se), colour = "black", width=.1, 
                position = pd) +
  geom_path(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Point", breaks = c("A", "B"),
                  labels = c("A - 2015", "B - 2015"), l = 40) +
  scale_y_continuous(breaks = 0:30*2, limits = c(0,14)) +
  #scale_x_discrete(labels=c("N1", "M1", "N2", "M2")) +
  geom_text(label="a. CW", x = 1.2, y = 14, size = 4, color = "black") +
  theme_bw()

LDO_CW

#Get STDevs by groups - DO
ldostat_MM <- summarySE(line.dat.DO[c(25:48,96:107),], measurevar = "DO", groupvars = c("Point","Time"))
ldostat_MM

#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LDO_MM <- ggplot(ldostat_MM, aes(x=Time, y=DO, 
                                 colour=Point, 
                                 group = Point)) + 
  ylab("DO (mg/L)") +
  geom_errorbar(aes(ymin=DO-se, ymax=DO+se), colour = "black", width=.1, 
                position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Point", breaks = c("A", "B_2015", "B_2016"),
                  labels = c("A - 2015", "B - 2015", "B - 2016"), l = 40) +
  scale_y_continuous(breaks = 0:30*2, limits = c(0,14)) +
  geom_text(label="b. MM", x = 1.2, y = 14, size = 4, color = "black") +
  theme_bw()

LDO_MM

#Get STDevs by groups - DO
ldostat_PM <- summarySE(line.dat.DO[49:71,], measurevar = "DO", groupvars = c("Point","Time"))
ldostat_PM

#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LDO_PM <- ggplot(ldostat_PM, aes(x=Time, y=DO, colour=Point, group = Point)) + 
  ylab("DO (mg/L)") +
  geom_errorbar(aes(ymin=DO-se, ymax=DO+se), colour = "black", width=.1, 
                position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Point", breaks = c("A", "B"),
                  labels = c("A - 2015", "B - 2015"), l = 40) +
  scale_y_continuous(breaks = 0:30*2, limits = c(0,14)) +
  geom_text(label="c. PM", x = 1.2, y = 14, size = 4, color = "black") +
  theme_bw()

LDO_PM

#Get STDevs by groups - DO
ldostat_SF <- summarySE(line.dat.DO[c(72:95,108:119),], measurevar = "DO", groupvars = c("Point","Time"))
ldostat_SF

#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LDO_SF <- ggplot(ldostat_SF, aes(x=Time, y=DO, colour=Point, group = Point)) + 
  ylab("DO (mg/L)") +
  geom_errorbar(aes(ymin=DO-se, ymax=DO+se), colour = "black", width=.1, 
                position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Point", breaks = c("A", "B_2015", "B_2016"),
                  labels = c("A - 2015", "B - 2015", "B - 2016"), l = 40) +
  scale_y_continuous(breaks = 0:30*2, limits = c(0,14)) +
  geom_text(label="d. SF", x = 1.2, y = 14, size = 4, color = "black") +
  theme_bw()

LDO_SF

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

LineCombo <- grid.arrange(arrangeGrob(LDO_CW,LDO_MM,LDO_PM,LDO_SF,nrow=2))

#pH line graphs
line.dat.pH <- read.table("chem_phys_dat_allsites_clean_lines3.txt",header=TRUE, row.names = 1, sep="\t")
head(line.dat.pH)
line.dat.pH <- line.dat.pH[,c(2:5,8,16)]
head(line.dat.pH)

#Get STDevs by groups - DO
lpHstat_CW <- summarySE(line.dat.pH[1:24,], measurevar = "pH", groupvars = c("Point","Time","Year"))
lpHstat_CW


#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LpH_CW <- ggplot(lpHstat_CW, aes(x=Time, y=pH, colour=Point, group = Point)) + 
  ylab("pH") +
  labs(x = "Time") +
  geom_errorbar(aes(ymin=pH-se, ymax=pH+se), colour = "black", width=.1, 
                position = pd) +
  geom_path(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Point", breaks = c("A", "B"),
                  labels = c("A - 2015", "B - 2015"), l = 40) +
  scale_y_continuous(breaks = 0:30*2, limits = c(5,9)) +
  #scale_x_discrete(labels=c("N1", "M1", "N2", "M2")) +
  #geom_text(label="a. CW", x = 1, y = 9, size = 4, color = "black") +
  theme_bw()

LpH_CW

#Get STDevs by groups - DO
lpHstat_MM <- summarySE(line.dat.pH[c(25:48,96:107),], measurevar = "pH", groupvars = c("Point","Time"))
lpHstat_MM

#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LpH_MM <- ggplot(lpHstat_MM, aes(x=Time, y=pH, 
                                 colour=Point, 
                                 group = Point)) + 
  ylab("pH") +
  geom_errorbar(aes(ymin=pH-se, ymax=pH+se), colour = "black", width=.1, 
                position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Point", breaks = c("A", "B_2015", "B_2016"),
                  labels = c("A - 2015", "B - 2015", "B - 2016"), l = 40) +
  scale_y_continuous(breaks = 0:30*2, limits = c(5,9)) +
  #geom_text(label="b. MM", x = 1, y = 9, size = 4, color = "black") +
  theme_bw()

LpH_MM

#Get STDevs by groups - DO
lpHstat_PM <- summarySE(line.dat.pH[49:71,], measurevar = "pH", groupvars = c("Point","Time"))
lpHstat_PM

#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LpH_PM <- ggplot(lpHstat_PM, aes(x=Time, y=pH, colour=Point, group = Point)) + 
  ylab("pH") +
  geom_errorbar(aes(ymin=pH-se, ymax=pH+se), colour = "black", width=.1, 
                position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Point", breaks = c("A", "B"),
                  labels = c("A - 2015", "B - 2015"), l = 40) +
  scale_y_continuous(breaks = 0:30*2, limits = c(5,9)) +
  #geom_text(label="c. PM", x = 1, y = 9, size = 4, color = "black") +
  theme_bw()

LpH_PM

#Get STDevs by groups - DO
lpHstat_SF <- summarySE(line.dat.pH[c(72:95,108:119),], measurevar = "pH", groupvars = c("Point","Time"))
lpHstat_SF

#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LpH_SF <- ggplot(lpHstat_SF, aes(x=Time, y=pH, colour=Point, group = Point)) + 
  ylab("pH") +
  geom_errorbar(aes(ymin=pH-se, ymax=pH+se), colour = "black", width=.1, 
                position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Point", breaks = c("A", "B_2015", "B_2016"),
                  labels = c("A - 2015", "B - 2015", "B - 2016"), l = 40) +
  scale_y_continuous(breaks = 0:30*2, limits = c(5,9)) +
  #geom_text(label="d. SF", x = 1, y = 9, size = 4, color = "black") +
  theme_bw()

LpH_SF

LineCombopH <- grid.arrange(arrangeGrob(LpH_CW,LpH_MM,LpH_PM,LpH_SF,nrow=2))

#Temp line graphs
line.dat.Temp <- read.table("chem_phys_dat_allsites_clean_lines3.txt",header=TRUE, row.names = 1, sep="\t")
head(line.dat.Temp)
line.dat.Temp <- line.dat.Temp[,c(2:6,16)]
head(line.dat.Temp)

#Get STDevs by groups - DO
lTempstat_CW <- summarySE(line.dat.Temp[1:24,], measurevar = "Temp", groupvars = c("Point","Time","Year"))
lTempstat_CW


#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LTemp_CW <- ggplot(lTempstat_CW, aes(x=Time, y=Temp, colour=Point, group = Point)) + 
  ylab("Temp (°C)") +
  labs(x = "Time") +
  geom_errorbar(aes(ymin=Temp-se, ymax=Temp+se), colour = "black", width=.1, 
                position = pd) +
  geom_path(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Point", breaks = c("A", "B"),
                  labels = c("A - 2015", "B - 2015"), l = 40) +
  scale_y_continuous(breaks = 0:30*2, limits = c(15,30)) +
  #scale_x_discrete(labels=c("N1", "M1", "N2", "M2")) +
  #geom_text(label="a. CW", x = 1, y = 30, size = 4, color = "black") +
  theme_bw()

LTemp_CW

#Get STDevs by groups - DO
lTempstat_MM <- summarySE(line.dat.Temp[c(25:48,96:107),], measurevar = "Temp", groupvars = c("Point","Time"))
lTempstat_MM

#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LTemp_MM <- ggplot(lTempstat_MM, aes(x=Time, y=Temp, 
                                 colour=Point, 
                                 group = Point)) + 
  ylab("Temp (°C)") +
  geom_errorbar(aes(ymin=Temp-se, ymax=Temp+se), colour = "black", width=.1, 
                position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Point", breaks = c("A", "B_2015", "B_2016"),
                  labels = c("A - 2015", "B - 2015", "B - 2016"), l = 40) +
  scale_y_continuous(breaks = 0:30*2, limits = c(15,30)) +
  #geom_text(label="b. MM", x = 1, y = 30, size = 4, color = "black") +
  theme_bw()

LTemp_MM

#Get STDevs by groups - DO
lTempstat_PM <- summarySE(line.dat.Temp[49:71,], measurevar = "Temp", groupvars = c("Point","Time"))
lTempstat_PM

#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LTemp_PM <- ggplot(lTempstat_PM, aes(x=Time, y=Temp, colour=Point, group = Point)) + 
  ylab("Temp (°C)") +
  geom_errorbar(aes(ymin=Temp-se, ymax=Temp+se), colour = "black", width=.1, 
                position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Point", breaks = c("A", "B"),
                  labels = c("A - 2015", "B - 2015"), l = 40) +
  scale_y_continuous(breaks = 0:30*2, limits = c(15,30)) +
  #geom_text(label="c. PM", x = 1, y = 30, size = 4, color = "black") +
  theme_bw()

LTemp_PM

#Get STDevs by groups - DO
lTempstat_SF <- summarySE(line.dat.Temp[c(72:95,108:119),], measurevar = "Temp", groupvars = c("Point","Time"))
lTempstat_SF

#plot
pd <- position_dodge(0.1) # move them .05 to the left and right to avoid overlaps
LTemp_SF <- ggplot(lTempstat_SF, aes(x=Time, y=Temp, colour=Point, group = Point)) + 
  ylab("Temp (°C)") +
  geom_errorbar(aes(ymin=Temp-se, ymax=Temp+se), colour = "black", width=.1, 
                position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") +
  scale_color_hue(name = "Point", breaks = c("A", "B_2015", "B_2016"),
                  labels = c("A - 2015", "B - 2015", "B - 2016"), l = 40) +
  scale_y_continuous(breaks = 0:30*2, limits = c(15,30)) +
  #geom_text(label="d. SF", x = 1, y = 30, size = 4, color = "black") +
  theme_bw() +
  theme(legend.position = c(0.7,0.235))

LTemp_SF

LineComboTemp <- grid.arrange(arrangeGrob(LTemp_CW,LTemp_MM,LTemp_PM,LTemp_SF,nrow=2))

SuperDuperGraph <- grid.arrange(LDO_CW+ theme(legend.position = 'none',axis.title.x=element_blank()),
                                LDO_MM+ theme(legend.position = 'none',axis.title.y=element_blank(),axis.title.x=element_blank()),
                                LDO_PM+ theme(legend.position = 'none',axis.title.y=element_blank(),axis.title.x=element_blank()),
                                LDO_SF+ theme(legend.position = 'none',axis.title.y=element_blank(),axis.title.x=element_blank()),
                                LpH_CW+ theme(legend.position = 'none',axis.title.x=element_blank()),
                                LpH_MM+ theme(legend.position = 'none',axis.title.y=element_blank(),axis.title.x=element_blank()),
                                LpH_PM+ theme(legend.position = 'none',axis.title.y=element_blank(),axis.title.x=element_blank()),
                                LpH_SF+ theme(legend.position = 'none',axis.title.y=element_blank(),axis.title.x=element_blank()),
                                LTemp_CW+ theme(legend.position = 'none'),
                                LTemp_MM+ theme(legend.position = 'none',axis.title.y=element_blank()),
                                LTemp_PM+ theme(legend.position = 'none',axis.title.y=element_blank()),
                                LTemp_SF+ theme(axis.title.y=element_blank()),
                                nrow=3)


#######################
### Alpha Diversity ###
#######################


### rarefaction curves
alpha_otu <- read.table("diel_t_singdoubrem.txt")
#alpha_otu2 <- alpha_otu[,-c(1,3)]
#talpha_OTU <- t(alpha_otu2)
#write.csv(talpha_OTU,"talpha_OTU.csv")
#alpha_otu <- read.csv("talpha_OTU.csv", header = TRUE, row.names = 1)
alpha_otut <- t(alpha_otu)
par(pty="s")
alpharare <- rarecurve(alpha_otut, step = 10000, col = 2:75, label = FALSE,
                       ylab = "OTUs")
#abline(v = 48226)
abline(a = 0, b = 1, lty = 2)


###########################
###########################
### NMDS Diel All Sites ###
###########################
###########################


##read in and transpose, 
#dieltransposed <- t(read.table("diel.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.abund.shared"))
#write.table(dieltransposed, "diel_t_singdoubrem.txt")

# PM_B1_N1 seems to be an outlier according to taxonomic phylotype breakdown.
# 45-50% Actinobacteria and Firmicutes. Removing from further analyses.

NUT <- read.table("chem_phys_dat_allsites_clean_deseq.txt", header = TRUE, row.names = 1, sep ="\t")
OTU <- read.table("diel_t_singdoubrem.txt", header = TRUE, row.names = 1, sep="\t")
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

#################################################
##### Normalizing and transforming the data #####
#################################################

### Deseq2 Normalization
### convert the counts to integer mode and group data to normalize together(e.g. region)
dds_wetland <- phyloseq_to_deseq2(ALL, ~ Wetland)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans_wetland = apply(counts(dds_wetland), 1, gm_mean)
dds_wetland = estimateSizeFactors(dds_wetland, geoMeans=geoMeans_wetland)
dds_wetland = estimateDispersions(dds_wetland)

### variance stabilizing transformation
#Make it a OTU again for phyloseq (no tax or sample data)
vst_wetland <- varianceStabilizingTransformation(dds_wetland, blind=FALSE)
vstMat_wetland <- assay(vst_wetland)
vstMat_wetland[vstMat_wetland<0]<-0
vst.otu.wetland <- otu_table(vstMat_wetland, taxa_are_rows=TRUE)
write.csv(otu_table(vst.otu.wetland), "vst.otu.wetland.csv")



###Create a distance matrix 
wetland_Dist <- phyloseq::distance(vst.otu.wetland, method= 'bray')

##NMDS of all samples
NMDS_wetland <- ordinate(vst.otu.wetland, method= 'NMDS', distance = wetland_Dist, formula = NULL)

## plot quick NMDS for comparisons
NMDS_wetland_plot <- plot_ordination(ALL, NMDS_wetland, color="Wetland", shape='Time')
NMDS_wetland_plot + geom_point(size=1) + theme_bw()


##################################################################
##### Publication-level NMDS figures and envfit correlations #####
##################################################################

##### Wetland NMDS
##### envfit depth
# NMDS_META <- read.table("GLCW_NMDS_metadata.txt", header = T, row.names = 1)
NMDS_META <- read.table("chem_phys_dat_allsites_clean_deseq.txt", header = T, row.names = 1)
ef_wetland <- envfit(NMDS_wetland, NMDS_META[,c(5:14)], permu=999, na.rm = TRUE)
ef_wetland

vectors_wetland<-as.data.frame(ef_wetland$vectors$arrows*sqrt(ef_wetland$vectors$r))
pvals_wetland=ef_wetland$vectors$pvals
r_wetland=ef_wetland$vectors$r
envfit_wetland<-cbind(vectors_wetland,pvals_wetland, r_wetland)
envfit_wetland <- envfit_wetland[c(1:5,7,9:10),]
envfit_wetland
#vectors are too long to fit on plot, transform with /2
envfit_wetland2 <- cbind(envfit_wetland[,1:2]/2, envfit_wetland[,3:4])
envfit_wetland2


NMDS_wetland_pub = data.frame(MDS1 = NMDS_wetland$points[,1], MDS2 = NMDS_wetland$points[,2])
NMDS_wetland_pub <- cbind(NMDS_wetland_pub,NMDS_META[,1:4])
NMDS_wetland_pub
NMDSplotWetland <- ggplot(NMDS_wetland_pub, aes(x=NMDS_wetland_pub$MDS1,
                                                y=NMDS_wetland_pub$MDS2),
                          groups = Wetland) +
  geom_point(aes(fill = factor(NMDS_wetland_pub$Wetland),
                 shape = NMDS_wetland_pub$Diel),
                 #color = NMDS_wetland_pub$Wetland),
             size=2) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  scale_shape_manual(name="Time",
                     labels=c("Dawn", "Dusk"),
                     values = c(21,22)) +
  scale_fill_hue(l = 40, name = "Wetland") +
  geom_text(label="Stress = 0.085", x = -0.4, y = -0.4, size = 4) +
  coord_equal(ylim = 0, xlim = 0) +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         color = guide_legend(title = "Wetland")) +
  theme_classic() +
  geom_segment(data = envfit_wetland2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_wetland2, aes(x=NMDS1*1.05, y=NMDS2*1.05, label=rownames(envfit_wetland2)),size=4)
NMDSplotWetland

### Test with perMANOVA

adonis_diel <- adonis(wetland_Dist ~ Wetland*Diel*Point, data = NMDS_META)
adonis_diel

#pairwise perMANOVA

pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 *** 0.001 ** 0.01 * 0.05 . 0.1  1")
  return(pairw.res)
  
} 

rdist_mat <- as.matrix(wetland_Dist)
pairwise.adonis(rdist_mat, factors = NMDS_META$Wetland)

memory.limit(size = 15000)
bd_dist <- vegdist(vst.otu.wetland, method = "bray")
bd_tax <- betadisper(bd_dist, group = NMDS_META$Wetland)
anova(bd_tax)

## remove incomplete cases
dat <- na.omit(rdist_mat)
## extract factor columns and drop redundant levels
fctr <- lapply(dat[sapply(dat, is.factor)], droplevels)
## count levels
sapply(fctr, nlevels)

pairwise.adonis(dat, factors = NMDS_META$Wetland)


###################
### Mantel Test ###
###################
### Testing for relationships between b-diversity and environ variables among wetlands
# create a distance matrix for nutrient data
num_nut <- NUT[5:14]
nut_dist <- vegdist(num_nut, method = "euclidean")

mantel_GL <- mantel(wetland_Dist,nut_dist, method = "spearman", permutations = 999)
mantel_GL


#######################
### Individual NMDS ###
#######################
# interaction term between wetland:point was significant, so we explore these wetlands
# on and individual basis.

### CW ###

###Create a distance matrix 
wetland_Dist_CW <- phyloseq::distance(vst.otu.wetland[,1:24], method= 'bray')

##NMDS of all samples
NMDS_wetland_CW <- ordinate(vst.otu.wetland[,1:24], method= 'NMDS', 
                            distance = wetland_Dist_CW, formula = NULL)
NMDS_wetland_CW
## plot quick NMDS for comparisons
NMDS_wetland_plot_CW <- plot_ordination(ALL, NMDS_wetland_CW, 
                                        color="Diel", label='Point')
NMDS_wetland_plot_CW + geom_point(size=1) + theme_bw()

##### envfit depth
ef_CW <- envfit(NMDS_wetland_CW, NMDS_META[1:24,5:14], permu=999, 
                na.rm = TRUE)
ef_CW

vectors_CW<-as.data.frame(ef_CW$vectors$arrows*sqrt(ef_CW$vectors$r))
pvals_CW=ef_CW$vectors$pvals
r_CW=ef_CW$vectors$r
envfit_CW<-cbind(vectors_CW,pvals_CW, r_CW)
envfit_CW
#keep only those vectors which significantly correlated to structure (p < 0.01)
envfit_CW <- envfit_CW[c(2:3,10),]
#vectors are too long to fit on plot, transform with /2
envfit_CW2 <- cbind(envfit_CW[,1:2]/3, envfit_CW[,3:4]/3)
envfit_CW2

NMDS_wetland_pub_CW = data.frame(MDS1 = NMDS_wetland_CW$points[,1],
                                 MDS2 = NMDS_wetland_CW$points[,2])
NMDS_wetland_pub_CW <- cbind(NMDS_wetland_pub_CW,NMDS_META[1:24,1:4])
NMDS_wetland_pub_CW
NMDSplotWetland_CW <- ggplot(NMDS_wetland_pub_CW, aes(x=NMDS_wetland_pub_CW$MDS1,
                                                      y=NMDS_wetland_pub_CW$MDS2)) +
  geom_point(aes(fill = factor(NMDS_wetland_pub_CW$Diel),
                 shape = NMDS_wetland_pub_CW$Point),
             size=5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  scale_shape_manual(name="Zone",
                     #labels=c("A", "B"),
                     values = c(21,22)) +
  scale_fill_hue(l = 40, name = "Time", labels = c("Dawn","Dusk")) +
  geom_text(label="Stress = 0.11", x = -0.35, y = -0.4, size = 4) +
  coord_equal(ylim = 0, xlim = 0) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  geom_text(label="a. CW", x = -0.4, y = 0.45, size = 4) +
  theme_classic() +
  geom_segment(data = envfit_CW2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_CW2, aes(x=NMDS1*1.05, y=NMDS2*1.05,
                                 label=rownames(envfit_CW2)),size=4)
NMDSplotWetland_CW

### MM ###

###Create a distance matrix 
wetland_Dist_MM <- phyloseq::distance(vst.otu.wetland[,25:48], method= 'bray')

##NMDS of all samples
NMDS_wetland_MM <- ordinate(vst.otu.wetland[,25:48], method= 'NMDS', 
                            distance = wetland_Dist_MM, formula = NULL)
NMDS_wetland_MM
## plot quick NMDS for comparisons
NMDS_wetland_plot_MM <- plot_ordination(ALL, NMDS_wetland_MM, 
                                        color="Diel", label='Point')
NMDS_wetland_plot_MM + geom_point(size=1) + theme_bw()

##### envfit depth
ef_MM <- envfit(NMDS_wetland_MM, NMDS_META[25:48,5:14], permu=999, 
                na.rm = TRUE)
ef_MM

vectors_MM<-as.data.frame(ef_MM$vectors$arrows*sqrt(ef_MM$vectors$r))
pvals_MM=ef_MM$vectors$pvals
r_MM=ef_MM$vectors$r
envfit_MM<-cbind(vectors_MM,pvals_MM, r_MM)
envfit_MM
#keep only those vectors which significantly correlated to structure
envfit_MM <- envfit_MM[10,]
#vectors are too long to fit on plot, transform with /2
envfit_MM2 <- cbind(envfit_MM[,1:2]/3, envfit_MM[,3:4]/3)
envfit_MM2

NMDS_wetland_pub_MM = data.frame(MDS1 = NMDS_wetland_MM$points[,1],
                                 MDS2 = NMDS_wetland_MM$points[,2])
NMDS_wetland_pub_MM <- cbind(NMDS_wetland_pub_MM,NMDS_META[25:48,1:4])
NMDS_wetland_pub_MM
NMDSplotWetland_MM <- ggplot(NMDS_wetland_pub_MM, aes(x=NMDS_wetland_pub_MM$MDS1,
                                                      y=NMDS_wetland_pub_MM$MDS2)) +
  geom_point(aes(fill = factor(NMDS_wetland_pub_MM$Diel),
                 shape = NMDS_wetland_pub_MM$Point), size=5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  scale_shape_manual(name="Zone",
                     values = c(21,22)) +
  scale_fill_hue(l = 40, name = "Time", labels = c("Dawn","Dusk")) +
  geom_text(label="Stress = 0.09", x = -0.35, y = -0.4, size = 4) +
  coord_equal(ylim = 0, xlim = 0) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  geom_text(label="b. MM", x = -0.4, y = 0.45, size = 4) +
  theme_classic() +
  geom_segment(data = envfit_MM2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_MM2, aes(x=NMDS1*1.05, y=NMDS2*1.05,
                                 label=rownames(envfit_MM2)),size=4)
NMDSplotWetland_MM


### PM ###

###Create a distance matrix 
wetland_Dist_PM <- phyloseq::distance(vst.otu.wetland[,49:71], method= 'bray')

##NMDS of all samples
NMDS_wetland_PM <- ordinate(vst.otu.wetland[,49:71], method= 'NMDS', 
                            distance = wetland_Dist_PM, formula = NULL)
NMDS_wetland_PM
## plot quick NMDS for comparisons
NMDS_wetland_plot_PM <- plot_ordination(ALL, NMDS_wetland_PM, 
                                        color="Diel", label='Point')
NMDS_wetland_plot_PM + geom_point(size=1) + theme_bw()

##### envfit depth
ef_PM <- envfit(NMDS_wetland_PM, NMDS_META[49:71,5:14], permu=999, 
                na.rm = TRUE)
ef_PM

vectors_PM<-as.data.frame(ef_PM$vectors$arrows*sqrt(ef_PM$vectors$r))
pvals_PM=ef_PM$vectors$pvals
r_PM=ef_PM$vectors$r
envfit_PM<-cbind(vectors_PM,pvals_PM, r_PM)
envfit_PM
#keep only those vectors which significantly correlated to structure
envfit_PM <- envfit_PM[c(2:3,9:10),]
#vectors are too long to fit on plot, transform with /2
envfit_PM2 <- cbind(envfit_PM[,1:2]/3, envfit_PM[,3:4]/3)
envfit_PM2

NMDS_wetland_pub_PM = data.frame(MDS1 = NMDS_wetland_PM$points[,1],
                                 MDS2 = NMDS_wetland_PM$points[,2])
NMDS_wetland_pub_PM <- cbind(NMDS_wetland_pub_PM,NMDS_META[49:71,1:4])
NMDS_wetland_pub_PM
NMDSplotWetland_PM <- ggplot(NMDS_wetland_pub_PM, aes(x=NMDS_wetland_pub_PM$MDS1,
                                                      y=NMDS_wetland_pub_PM$MDS2)) +
  geom_point(aes(fill = factor(NMDS_wetland_pub_PM$Diel),
                 shape = NMDS_wetland_pub_PM$Point),
             size=5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  scale_shape_manual(name="Zone",
                     values = c(21,22)) +
  scale_fill_hue(l = 40, name = "Time", labels = c("Dawn","Dusk")) +
  geom_text(label="Stress = 0.07", x = -0.35, y = -0.4, size = 4) +
  coord_equal(ylim = 0, xlim = 0) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  geom_text(label="c. PM", x = -0.4, y = 0.45, size = 4) +
  theme_classic() +
  geom_segment(data = envfit_PM2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_PM2, aes(x=NMDS1*1.05, y=NMDS2*1.05,
                                 label=rownames(envfit_PM2)),size=4)
NMDSplotWetland_PM

### SF ###

###Create a distance matrix 
wetland_Dist_SF <- phyloseq::distance(vst.otu.wetland[,72:95], method= 'bray')

##NMDS of all samples
NMDS_wetland_SF <- ordinate(vst.otu.wetland[,72:95], method= 'NMDS', 
                            distance = wetland_Dist_SF, formula = NULL)
NMDS_wetland_SF
## plot quick NMDS for comparisons
NMDS_wetland_plot_SF <- plot_ordination(ALL, NMDS_wetland_SF, 
                                        color="Diel", label='Point')
NMDS_wetland_plot_SF + geom_point(size=1) + theme_bw()

##### envfit depth
ef_SF <- envfit(NMDS_wetland_SF, NMDS_META[72:95,5:14], permu=999, 
                na.rm = TRUE)
ef_SF

vectors_SF<-as.data.frame(ef_SF$vectors$arrows*sqrt(ef_SF$vectors$r))
pvals_SF=ef_SF$vectors$pvals
r_SF=ef_SF$vectors$r
envfit_SF<-cbind(vectors_SF,pvals_SF, r_SF)
envfit_SF
#keep only those vectors which significantly correlated to structure
envfit_SF <- envfit_SF[9,]
#vectors are too long to fit on plot, transform with /3
envfit_SF2 <- cbind(envfit_SF[,1:2]/3, envfit_SF[,3:4]/3)
envfit_SF2

NMDS_wetland_pub_SF = data.frame(MDS1 = NMDS_wetland_SF$points[,1],
                                 MDS2 = NMDS_wetland_SF$points[,2])
NMDS_wetland_pub_SF <- cbind(NMDS_wetland_pub_SF,NMDS_META[72:95,1:4])
NMDS_wetland_pub_SF
NMDSplotWetland_SF <- ggplot(NMDS_wetland_pub_SF, aes(x=NMDS_wetland_pub_SF$MDS1,
                                                      y=NMDS_wetland_pub_SF$MDS2)) +
  geom_point(aes(fill = factor(NMDS_wetland_pub_SF$Diel),
                 shape = NMDS_wetland_pub_SF$Point),
             size=5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  scale_shape_manual(name="Zone",
                     values = c(21,22)) +
  scale_fill_hue(l = 40, name = "Time", labels = c("Dawn","Dusk")) +
  geom_text(label="Stress = 0.15", x = -0.35, y = -0.4, size = 4) +
  coord_equal(ylim = 0, xlim = 0) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  geom_text(label="d. SF", x = -0.4, y = 0.45, size = 4) +
  theme_classic() +
  geom_segment(data = envfit_SF2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_SF2, aes(x=NMDS1*1.05, y=NMDS2*1.05,
                                 label=rownames(envfit_SF2)),size=4)
NMDSplotWetland_SF


## Multiplot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
multiplot(NMDSplotWetland_CW,NMDSplotWetland_MM,NMDSplotWetland_PM,NMDSplotWetland_SF, cols = 2)


####################################
### Within Wetland Location NMDS ###
####################################

### CW_A ###

###Create a distance matrix 
wetland_Dist_CW_A <- phyloseq::distance(vst.otu.wetland[,1:12], method= 'bray')

##NMDS of all samples
NMDS_wetland_CW_A <- ordinate(vst.otu.wetland[,1:12], method= 'NMDS', 
                            distance = wetland_Dist_CW_A, formula = NULL)
NMDS_wetland_CW_A
## plot quick NMDS for comparisons
NMDS_wetland_plot_CW_A <- plot_ordination(ALL, NMDS_wetland_CW_A, 
                                        color="Diel", label='Point')
NMDS_wetland_plot_CW_A + geom_point(size=1) + theme_bw()

##### envfit depth
ef_CW_A <- envfit(NMDS_wetland_CW_A, NMDS_META[1:12,5:14], permu=999, 
                na.rm = TRUE)
ef_CW_A

vectors_CW_A<-as.data.frame(ef_CW_A$vectors$arrows*sqrt(ef_CW_A$vectors$r))
pvals_CW_A=ef_CW_A$vectors$pvals
r_CW_A=ef_CW_A$vectors$r
envfit_CW_A<-cbind(vectors_CW_A,pvals_CW_A, r_CW_A)
envfit_CW_A

## No relationships between NMDS and envfit


### CW_B ###

###Create a distance matrix 
wetland_Dist_CW_B <- phyloseq::distance(vst.otu.wetland[,13:24], method= 'bray')

##NMDS of all samples
NMDS_wetland_CW_B <- ordinate(vst.otu.wetland[,13:24], method= 'NMDS', 
                              distance = wetland_Dist_CW_B, formula = NULL)
NMDS_wetland_CW_B
## plot quick NMDS for comparisons
NMDS_wetland_plot_CW_B <- plot_ordination(ALL, NMDS_wetland_CW_B, 
                                          color="Diel", label='Point')
NMDS_wetland_plot_CW_B + geom_point(size=1) + theme_bw()

##### envfit depth
ef_CW_B <- envfit(NMDS_wetland_CW_B, NMDS_META[13:24,5:14], permu=999, 
                  na.rm = TRUE)
ef_CW_B

vectors_CW_B<-as.data.frame(ef_CW_B$vectors$arrows*sqrt(ef_CW_B$vectors$r))
pvals_CW_B=ef_CW_B$vectors$pvals
r_CW_B=ef_CW_B$vectors$r
envfit_CW_B<-cbind(vectors_CW_B,pvals_CW_B, r_CW_B)
envfit_CW_B
#keep only those vectors which significantly correlated to structure (p < 0.01)
envfit_CW_B <- envfit_CW_B[c(6),]
#vectors are too long to fit on plot, transform with /2
envfit_CW_B2 <- cbind(envfit_CW_B[,1:2]/3, envfit_CW_B[,3:4]/3)
envfit_CW_B2

NMDS_wetland_pub_CW_B = data.frame(MDS1 = NMDS_wetland_CW_B$points[,1],
                                   MDS2 = NMDS_wetland_CW_B$points[,2])
NMDS_wetland_pub_CW_B <- cbind(NMDS_wetland_pub_CW_B,NMDS_META[13:24,1:4])
NMDS_wetland_pub_CW_B
NMDSplotWetland_CW_B <- ggplot(NMDS_wetland_pub_CW_B, aes(x=NMDS_wetland_pub_CW_B$MDS1,
                                                          y=NMDS_wetland_pub_CW_B$MDS2)) +
  geom_point(aes(fill = factor(NMDS_wetland_pub_CW_B$Diel),
                 color = NMDS_wetland_pub_CW_B$Diel),
             size=5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Time") +
  scale_fill_brewer(type = "qual", palette = "Set1", guide = FALSE) +
  geom_text(label="Stress = ", x = -0.2, y = -0.4, size = 4) +
  coord_equal(ylim = 0, xlim = 0) +
  guides(fill = FALSE) +
  geom_text(label="a. CW_B", x = -0.36, y = 0.45, size = 4) +
  theme_classic()# +  
#guides(fill = guide_legend(override.aes= list(colour = c("white","grey","black")))) +
#theme(legend.key = element_rect(fill = "grey90"))
#NMDSplotWetland_CW_B <- NMDSplotWetland_CW_B + geom_text_repel(aes(label = NMDS_wetland_pub_CW_B$Diel))
NMDSplotWetland_CW_B <- NMDSplotWetland_CW_B + 
  geom_segment(data = envfit_CW_B, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_CW_B2, aes(x=NMDS1*1.05, y=NMDS2*1.05,
                                   label=rownames(envfit_CW_B2)),size=4)
NMDSplotWetland_CW_B


### MM_A ###

###Create a distance matrix 
wetland_Dist_MM_A <- phyloseq::distance(vst.otu.wetland[,25:36], method= 'bray')

##NMDS of all samples
NMDS_wetland_MM_A <- ordinate(vst.otu.wetland[,25:36], method= 'NMDS', 
                              distance = wetland_Dist_MM_A, formula = NULL)
NMDS_wetland_MM_A
## plot quick NMDS for comparisons
NMDS_wetland_plot_MM_A <- plot_ordination(ALL, NMDS_wetland_MM_A, 
                                          color="Diel", label='Point')
NMDS_wetland_plot_MM_A + geom_point(size=1) + theme_bw()

##### envfit depth
ef_MM_A <- envfit(NMDS_wetland_MM_A, NMDS_META[25:36,5:14], permu=999, 
                  na.rm = TRUE)
ef_MM_A

# No relationships between NMDS and envfit

### MM_B ###

###Create a distance matrix 
wetland_Dist_MM_B <- phyloseq::distance(vst.otu.wetland[,37:48], method= 'bray')

##NMDS of all samples
NMDS_wetland_MM_B <- ordinate(vst.otu.wetland[,37:48], method= 'NMDS', 
                              distance = wetland_Dist_MM_B, formula = NULL)
NMDS_wetland_MM_B
## plot quick NMDS for comparisons
NMDS_wetland_plot_MM_B <- plot_ordination(ALL, NMDS_wetland_MM_B, 
                                          color="Diel", label='Point')
NMDS_wetland_plot_MM_B + geom_point(size=1) + theme_bw()

##### envfit depth
ef_MM_B <- envfit(NMDS_wetland_MM_B, NMDS_META[37:48,5:14], permu=999, 
                  na.rm = TRUE)
ef_MM_B

vectors_MM_B<-as.data.frame(ef_MM_B$vectors$arrows*sqrt(ef_MM_B$vectors$r))
pvals_MM_B=ef_MM_B$vectors$pvals
r_MM_B=ef_MM_B$vectors$r
envfit_MM_B<-cbind(vectors_MM_B,pvals_MM_B, r_MM_B)
envfit_MM_B
#keep only those vectors which significantly correlated to structure (p < 0.01)
envfit_MM_B <- envfit_MM_B[c(1:4,7),]
#vectors are too long to fit on plot, transform with /2
envfit_MM_B2 <- cbind(envfit_MM_B[,1:2]/3, envfit_MM_B[,3:4]/3)
envfit_MM_B2

NMDS_wetland_pub_MM_B = data.frame(MDS1 = NMDS_wetland_MM_B$points[,1],
                                   MDS2 = NMDS_wetland_MM_B$points[,2])
NMDS_wetland_pub_MM_B <- cbind(NMDS_wetland_pub_MM_B,NMDS_META[37:48,1:4])
NMDS_wetland_pub_MM_B
NMDSplotWetland_MM_B <- ggplot(NMDS_wetland_pub_MM_B, aes(x=NMDS_wetland_pub_MM_B$MDS1,
                                                          y=NMDS_wetland_pub_MM_B$MDS2)) +
  geom_point(aes(fill = factor(NMDS_wetland_pub_MM_B$Diel),
                 color = NMDS_wetland_pub_MM_B$Diel),
             size=5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Time") +
  scale_fill_brewer(type = "qual", palette = "Set1", guide = FALSE) +
  geom_text(label="Stress = ", x = -0.2, y = -0.4, size = 4) +
  coord_equal(ylim = 0, xlim = 0) +
  guides(fill = FALSE) +
  geom_text(label="a. MM_B", x = -0.36, y = 0.45, size = 4) +
  theme_classic()# +  
#guides(fill = guide_legend(override.aes= list(colour = c("white","grey","black")))) +
#theme(legend.key = element_rect(fill = "grey90"))
#NMDSplotWetland_MM_B <- NMDSplotWetland_MM_B + geom_text_repel(aes(label = NMDS_wetland_pub_MM_B$Diel))
NMDSplotWetland_MM_B <- NMDSplotWetland_MM_B + 
  geom_segment(data = envfit_MM_B, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_MM_B2, aes(x=NMDS1*1.05, y=NMDS2*1.05,
                                   label=rownames(envfit_MM_B2)),size=4)
NMDSplotWetland_MM_B

### Significance found with DO, pH, Temp, and OrP all headed the same direction


### PM_A ###

###Create a distance matrix 
wetland_Dist_PM_A <- phyloseq::distance(vst.otu.wetland[,49:60], method= 'bray')

##NMDS of all samples
NMDS_wetland_PM_A <- ordinate(vst.otu.wetland[,49:60], method= 'NMDS', 
                              distance = wetland_Dist_PM_A, formula = NULL)
NMDS_wetland_PM_A
## plot quick NMDS for comparisons
NMDS_wetland_plot_PM_A <- plot_ordination(ALL, NMDS_wetland_PM_A, 
                                          color="Diel", label='Point')
NMDS_wetland_plot_PM_A + geom_point(size=1) + theme_bw()

##### envfit depth
ef_PM_A <- envfit(NMDS_wetland_PM_A, NMDS_META[49:60,5:14], permu=999, 
                  na.rm = TRUE)
ef_PM_A

### No relationships between NMDS and envfit


### PM_B ###

###Create a distance matrix 
wetland_Dist_PM_B <- phyloseq::distance(vst.otu.wetland[,61:71], method= 'bray')

##NMDS of all samples
NMDS_wetland_PM_B <- ordinate(vst.otu.wetland[,61:71], method= 'NMDS', 
                              distance = wetland_Dist_PM_B, formula = NULL)
NMDS_wetland_PM_B
## plot quick NMDS for comparisons
NMDS_wetland_plot_PM_B <- plot_ordination(ALL, NMDS_wetland_PM_B, 
                                          color="Diel", label='Point')
NMDS_wetland_plot_PM_B + geom_point(size=1) + theme_bw()

##### envfit depth
ef_PM_B <- envfit(NMDS_wetland_PM_B, NMDS_META[61:71,5:14], permu=999, 
                  na.rm = TRUE)
ef_PM_B

### No relationships between NMDS and envfit


### SF_A ###

###Create a distance matrix 
wetland_Dist_SF_A <- phyloseq::distance(vst.otu.wetland[,72:83], method= 'bray')

##NMDS of all samples
NMDS_wetland_SF_A <- ordinate(vst.otu.wetland[,72:83], method= 'NMDS', 
                              distance = wetland_Dist_SF_A, formula = NULL)
NMDS_wetland_SF_A
## plot quick NMDS for comparisons
NMDS_wetland_plot_SF_A <- plot_ordination(ALL, NMDS_wetland_SF_A, 
                                          color="Diel", label='Point')
NMDS_wetland_plot_SF_A + geom_point(size=1) + theme_bw()

##### envfit depth
ef_SF_A <- envfit(NMDS_wetland_SF_A, NMDS_META[72:83,5:14], permu=999, 
                  na.rm = TRUE)
ef_SF_A

### No relationships between NMDS and envfit


### SF_B ###

###Create a distance matrix 
wetland_Dist_SF_B <- phyloseq::distance(vst.otu.wetland[,84:95], method= 'bray')

##NMDS of all samples
NMDS_wetland_SF_B <- ordinate(vst.otu.wetland[,84:95], method= 'NMDS', 
                              distance = wetland_Dist_SF_B, formula = NULL)
NMDS_wetland_SF_B
## plot quick NMDS for comparisons
NMDS_wetland_plot_SF_B <- plot_ordination(ALL, NMDS_wetland_SF_B, 
                                          color="Diel", label='Point')
NMDS_wetland_plot_SF_B + geom_point(size=1) + theme_bw()

##### envfit depth
ef_SF_B <- envfit(NMDS_wetland_SF_B, NMDS_META[84:95,5:14], permu=999, 
                  na.rm = TRUE)
ef_SF_B

### No relationships between NMDS and envfit


##################
### perMANOVAs ###
##################

##############
### ADONIS ###
##############

#CW
#perMANOVA
NMDS_META_CW <- NMDS_META[1:24,]
CW_adonis <- adonis(wetland_Dist_CW ~ Diel*Point*Time, data = NMDS_META_CW)
CW_adonis

# Point significant

#CW_A
#perMANOVA
NMDS_META_CW_A <- NMDS_META[1:12,]
CW_A_adonis <- adonis(wetland_Dist_CW_A ~ Diel*Time, data = NMDS_META_CW_A)
CW_A_adonis

#No significance

#CW_B
#perMANOVA
NMDS_META_CW_B <- NMDS_META[13:24,]
CW_B_adonis <- adonis(wetland_Dist_CW_B ~ Diel*Time, data = NMDS_META_CW_B)
CW_B_adonis

#no significance

#MM
#perMANOVA
NMDS_META_MM <- NMDS_META[25:48,]
MM_adonis <- adonis(wetland_Dist_MM ~ Diel*Point*Time, data = NMDS_META_MM)
MM_adonis

# Point Significance (p < 0.001) Diel significance (p <0.05)

#MM_A
#perMANOVA
NMDS_META_MM_A <- NMDS_META[25:36,]
MM_A_adonis <- adonis(wetland_Dist_MM_A ~ Diel*Time, data = NMDS_META_MM_A)
MM_A_adonis

# Time significance (p < 0.05)

#MM_B
#perMANOVA
NMDS_META_MM_B <- NMDS_META[37:48,]
MM_B_adonis <- adonis(wetland_Dist_MM_B ~ Diel*Time, data = NMDS_META_MM_B)
MM_B_adonis

# Diel significance (p < 0.01)

#PM
#perMANOVA
NMDS_META_PM <- NMDS_META[49:71,]
PM_adonis <- adonis(wetland_Dist_PM ~ Diel*Point*Time, data = NMDS_META_PM)
PM_adonis

# Point Significance (p < 0.001)

#PM_A
#perMANOVA
NMDS_META_PM_A <- NMDS_META[49:60,]
PM_A_adonis <- adonis(wetland_Dist_PM_A ~ Diel*Time, data = NMDS_META_PM_A)
PM_A_adonis

#no significance

#PM_B
#perMANOVA
NMDS_META_PM_B <- NMDS_META[61:71,]
PM_B_adonis <- adonis(wetland_Dist_PM_B ~ Diel*Time, data = NMDS_META_PM_B)
PM_B_adonis

#no significance

#SF
#perMANOVA
NMDS_META_SF <- NMDS_META[72:95,]
SF_adonis <- adonis(wetland_Dist_SF ~ Diel*Point*Time, data = NMDS_META_SF)
SF_adonis

# Point Significance (p < 0.01)

#SF_A
#perMANOVA
NMDS_META_SF_A <- NMDS_META[72:83,]
SF_A_adonis <- adonis(wetland_Dist_SF_A ~ Diel*Time, data = NMDS_META_SF_A)
SF_A_adonis

#no significance

#SF_B
#perMANOVA
NMDS_META_SF_B <- NMDS_META[84:95,]
SF_B_adonis <- adonis(wetland_Dist_SF_B ~ Diel*Time, data = NMDS_META_SF_B)
SF_B_adonis

#no significance


############
### DGGE ###
############

MM_DGGE <- read.csv("MM_DGGE_interpret.csv", header=TRUE)
NM_DGGE <- read.csv("NM_DGGE_interpret.csv", header=TRUE)

MM_DGGE_dist <- vegdist(MM_DGGE[,3:34], method="jaccard", binary = TRUE)
NM_DGGE_dist <- vegdist(NM_DGGE[,3:35], method="jaccard", binary = TRUE)

adonis(MM_DGGE_dist ~ Time, data = MM_DGGE)
adonis(NM_DGGE_dist ~ Time, data = NM_DGGE)
