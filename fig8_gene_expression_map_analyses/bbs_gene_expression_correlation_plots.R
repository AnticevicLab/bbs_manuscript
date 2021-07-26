rm(list = ls())
setwd("../fig8_gene_expression_map_analyses")
library(ggplot2)
library(stats)
library(PerformanceAnalytics)
library(grDevices)

# -- Set ggplot basic theme
theme1<- theme(panel.border = element_rect(size = 2, colour = 'NA', fill='NA'), 
               axis.line = element_line(colour = "black", size = 1),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               plot.background = element_rect(fill = "transparent", colour = NA),
               panel.background = element_rect(fill = "white", colour = NA),
               legend.position="none",
               axis.ticks = element_line(colour = "black", size=1),
               axis.ticks.length = unit(.4, "cm"),
               axis.text.x = element_text(size = 30, color="black"), 
               axis.title.x = element_text(size =30,margin=margin(-20,0,0,0)),
               axis.text.y = element_text(size = 30, color="black",angle=90), 
               axis.title.y = element_text(size =35,margin=margin(0,-20,0,0)))

##### Read data -- table of correlations between PC maps and gene expression maps
genecor = read.table("geneSweepResults_parcellated_cortex.csv", sep=",", header=TRUE)

# -- Format table
signs <- c(1,1,1,1,1)
pccor <- data.frame(cbind(data.frame("Gene"=genecor[,1]), t(signs*t(genecor[,2:6]))))
names(pccor) <- c("Gene","PC1","PC2","PC3","PC4","PC5")

##### Main text figure
# -- List genes of interest (GOI)
GOI1="^SST$"
GOI2="^PVALB$"
GOI3="^GABRA1$"
GOI4="^GABRA5$"
GOI5="^HTR2A$"
GOI6="^HTR1E$" 
GOI7="^HTR2C$"

##### Compute the percentile for each GOI relative to the entire distribution of correlation values for that map
# -- Compute for the 5 PC maps
for (pcno in 1:5){
  pc <- paste0("pccor$PC",pcno)
  for (gno in 1:7){
    assign(paste0("PC",pcno,"GOIr", gno), c(eval(parse(text=pc)))[grep(get(paste0("GOI", gno)),pccor$Gene)])
    print(sum(c(eval(parse(text=pc))) <= get(paste0("PC",pcno,"GOIr", gno)), na.rm=TRUE) / nrow(pccor))
  }
}

# -- Plot distribution of genes with GOI highlighted
pdf("GeneExpressionCorrelation_PC3.pdf", height=4, width=12)
q <- ggplot(pccor, aes(x=get("PC3"))) +
  geom_histogram(colour="grey50", fill="grey50",binwidth=0.01) +
  geom_vline(aes(xintercept=mean(get("PC3"), na.rm=T)),color="black", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC3GOIr",1))),color="#e66101", linetype="dashed", size=3) +
  geom_vline(aes(xintercept=get(paste0("PC3GOIr",2))),color="#fdb863", linetype="dashed", size=3) +
  geom_vline(aes(xintercept=get(paste0("PC3GOIr",3))),color="#8073ac", linetype="dashed", size=3) +
  geom_vline(aes(xintercept=get(paste0("PC3GOIr",4))),color="#542788", linetype="dashed", size=3) +
  geom_vline(aes(xintercept=get(paste0("PC3GOIr",5))),color="#006837", linetype="dashed", size=3) +
  geom_vline(aes(xintercept=get(paste0("PC3GOIr",6))),color="#66bd63", linetype="dashed", size=3) +
  geom_vline(aes(xintercept=get(paste0("PC3GOIr",7))),color="#1a9850", linetype="dashed", size=3) +
  ylab("") +
  xlab("Pearson's r") +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-0.45,0.45), breaks=c(-0.4,0,0.4), labels=c(-0.4,0,0.4)) +
  scale_y_continuous(expand = c(0.001, 0), limits=c(0,930), breaks=c(0,400,800), labels=c(0,400,800)) +
  theme1 +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.x = element_text(size = 20, color="black"),
        axis.title.x = element_text(size =30,margin=margin(0,0,0,0)))
print(q)
dev.off()


##### Plot correlation scatterplots for GOI
# -- Read data
pcmaps <- read.table("../fig3_univariate_mapping_analyses/univariate_maps/GBC_PC3_ztstat.pscalar.txt") # -- Get values of 718 parcels
names(pcmaps) <- c("PC3")
pcmapssymm <- (pcmaps[1:180,] + pcmaps[181:360,])/2 # -- Symmetrize the cortex only

# -- Get values of gene expression for each GOI, using gene expression data from Burt et al. (2018)
genemaps <- vector()
for (genemap in c("SST","PVALB","GABRA1","GABRA5","HTR2A","HTR1E","HTR2C")){
  genename <- paste0(genemap, "map")
  assign(genename, as.data.frame(read.table(paste0("Burt_et_al_2018_Hierarchy_Kx5n/",genemap,".txt")))) # -- 
  genemaps <- cbind(genemaps,get(genename)$V1)
}

# -- Organize into dataframe
genemaps <- as.data.frame(genemaps)
names(genemaps) <- c("SST","PVALB","GABRA1","GABRA5","HTR2A","HTR1E","HTR2C")
pcgenedf <- data.frame("PC3"=pcmapssymm, genemaps)

# -- Plot the correlations between PC3 map and each GOI expression map
pdf("GeneExpression_Correlation_Scatterplot_SST.pdf", width=6, height=6)
ggplot(dat=pcgenedf, aes(x=PC3,y=SST)) +
  geom_point(col="#DE6A10", size=4, alpha=0.8) +
  geom_smooth(method='lm',formula=y~x, se=FALSE, lwd=5, col="black") +
  scale_y_continuous(expand=c(0.001,0), limits=c(-3.4,2.8), breaks=seq(-3,2,length=2),labels=seq(-3,2,length=2)) +
  scale_x_continuous(expand=c(0.001,0), limits=c(-4.4,4.4), labels=format(seq(-4,4,length=2)), breaks=seq(-4,4,length=2)) +
  xlab("PC3") +
  ylab("Gene") +
  theme(aspect.ratio=1) +
  theme1
dev.off()
  
pdf("GeneExpression_Correlation_Scatterplot_PVALB.pdf", width=6, height=6)
ggplot(dat=pcgenedf, aes(x=PC3,y=PVALB)) +
  geom_point(col="#FBB76A", size=4, alpha=0.8) +
  geom_smooth(method='lm',formula=y~x, se=FALSE, lwd=5, col="black") +
  scale_y_continuous(expand=c(0.001,0), limits=c(-3.4,2.8), breaks=seq(-3,2,length=2),labels=c("","")) +
  scale_x_continuous(expand=c(0.001,0), limits=c(-4.4,4.4), labels=format(seq(-4,4,length=2)), breaks=seq(-4,4,length=2)) +
  xlab("PC3") +
  ylab("Gene") +
  theme(aspect.ratio=1) +
  theme1
dev.off()

pdf("GeneExpression_Correlation_Scatterplot_GABRA1.pdf", width=6, height=6)
ggplot(dat=pcgenedf, aes(x=PC3,y=GABRA1)) +
  geom_point(col="#8074AA", size=4, alpha=0.8) +
  geom_smooth(method='lm',formula=y~x, se=FALSE, lwd=5, col="black") +
  scale_y_continuous(expand=c(0.001,0), limits=c(-3.4,2.8), breaks=seq(-3,2,length=2),labels=seq(-3,2,length=2)) +
  scale_x_continuous(expand=c(0.001,0), limits=c(-4.4,4.4), labels=format(seq(-4,4,length=2)), breaks=seq(-4,4,length=2)) +
  xlab("PC3") +
  ylab("Gene") +
  theme(aspect.ratio=1) +
  theme1
dev.off()

pdf("GeneExpression_Correlation_Scatterplot_GABRA5.pdf", width=6, height=6)
ggplot(dat=pcgenedf, aes(x=PC3,y=GABRA5)) +
  geom_point(col="#542B86", size=4, alpha=0.8) +
  geom_smooth(method='lm',formula=y~x, se=FALSE, lwd=5, col="black") +
  scale_y_continuous(expand=c(0.001,0), limits=c(-3.4,2.8), breaks=seq(-3,2,length=2),labels=c("","")) +
  scale_x_continuous(expand=c(0.001,0), limits=c(-4.4,4.4), labels=format(seq(-4,4,length=2)), breaks=seq(-4,4,length=2)) +
  xlab("PC3") +
  ylab("Gene") +
  theme(aspect.ratio=1) +
  theme1
dev.off()

pdf("GeneExpression_Correlation_Scatterplot_HTR2A.pdf", width=6, height=6)
ggplot(dat=pcgenedf, aes(x=PC3,y=HTR2A)) +
  geom_point(col="#0A6739", size=4, alpha=0.8) +
  geom_smooth(method='lm',formula=y~x, se=FALSE, lwd=5, col="black") +
  scale_y_continuous(expand=c(0.001,0), limits=c(-4.4,3.4), breaks=seq(-4,3,length=2),labels=c("","")) +
  scale_x_continuous(expand=c(0.001,0), limits=c(-4.4,4.4), labels=format(seq(-4,4,length=2)), breaks=seq(-4,4,length=2)) +
  xlab("PC3") +
  ylab("Gene") +
  theme(aspect.ratio=1) +
  theme1
dev.off()

pdf("GeneExpression_Correlation_Scatterplot_HTR2C.pdf", width=6, height=6)
ggplot(dat=pcgenedf, aes(x=PC3,y=HTR2C)) +
  geom_point(col="#00882B", size=4, alpha=0.8) +
  geom_smooth(method='lm',formula=y~x, se=FALSE, lwd=5, col="black") +
  scale_y_continuous(expand=c(0.001,0), limits=c(-4.4,3.4), breaks=seq(-4,3,length=2),labels=c("","")) +
  scale_x_continuous(expand=c(0.001,0), limits=c(-4.4,4.4), labels=format(seq(-4,4,length=2)), breaks=seq(-4,4,length=2)) +
  xlab("PC3") +
  ylab("Gene") +
  theme(aspect.ratio=1) +
  theme1
dev.off()

pdf("GeneExpression_Correlation_Scatterplot_HTR1E.pdf", width=6, height=6)
ggplot(dat=pcgenedf, aes(x=PC3,y=HTR1E)) +
  geom_point(col="#69BC67", size=4, alpha=0.8) +
  geom_smooth(method='lm',formula=y~x, se=FALSE, lwd=5, col="black") +
  scale_y_continuous(expand=c(0.001,0), limits=c(-4.4,3.4), breaks=seq(-4,3,length=2),labels=seq(-4,3,length=2)) +
  scale_x_continuous(expand=c(0.001,0), limits=c(-4.4,4.4), labels=format(seq(-4,4,length=2)), breaks=seq(-4,4,length=2)) +
  xlab("PC3") +
  ylab("Gene") +
  theme(aspect.ratio=1) +
  theme1
dev.off()

##### Supplementary figure -- distributions for all PCs
# -- Select genes of interest (GOI) for each PC
# -- PC1
pcno=1
GOI1="^SST$"
GOI2="^PVALB$"
GOI3="^GABRA3$"
GOI4="^GABRA1$"
GOI5="^HTR2A$"
GOI6="^HTR1E$" #pos
GOI7="^HTR1F$" #neg

# -- Compute percentiles
pc <- paste0("pccor$PC",pcno)
for (gno in 1:7){
  assign(paste0("PC",pcno,"GOIr", gno), c(eval(parse(text=pc)))[grep(get(paste0("GOI", gno)),pccor$Gene)])
  print(sum(c(eval(parse(text=pc))) <= get(paste0("PC",pcno,"GOIr", gno)), na.rm=TRUE) / nrow(pccor))
}

# -- Plot distributions with GOI highlighted
pdf("GeneExpressionCorrelation_PC1.pdf", height=4, width=12)
q <- ggplot(pccor, aes(x=eval(parse(text=paste0("pccor$PC",pcno))))) +
  geom_histogram(colour="grey50", fill="grey50",binwidth=0.01) +
  geom_vline(aes(xintercept=mean(get(paste0("PC",pcno)), na.rm=T)),color="black", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",2))),color="#fdb863", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",3))),color="#8073ac", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",4))),color="#542788", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",5))),color="#006837", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",6))),color="#66bd63", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",7))),color="#1a9850", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",1))),color="#e66101", linetype="dashed", size=2) +
  ylab("") +
  xlab("") +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-0.8,0.8), breaks=c(-0.5,0,0.5), labels=c(-0.5,0,0.5)) +
  scale_y_continuous(expand = c(0.001, 0), limits=c(0,930), breaks=c(0,400,800), labels=c(0,400,800)) +
  theme1 +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.x = element_text(size = 20, color="black"))
print(q)
dev.off()

# -- PC2
pcno=2
GOI1="^SST$"
GOI2="^PVALB$"
GOI3="^GABRA1$" #pos
GOI4="^GABRA2$" #neg
GOI5="^HTR2A$"
GOI6="^HTR1F$" #pos
GOI7="^HTR5A$" #neg

# -- Compute percentiles
pc <- paste0("pccor$PC",pcno)
for (gno in 1:7){
  assign(paste0("PC",pcno,"GOIr", gno), c(eval(parse(text=pc)))[grep(get(paste0("GOI", gno)),pccor$Gene)])
  print(sum(c(eval(parse(text=pc))) <= get(paste0("PC",pcno,"GOIr", gno)), na.rm=TRUE) / nrow(pccor))
}

# -- Plot distributions with GOI highlighted
pdf("GeneExpressionCorrelation_PC2.pdf", height=4, width=12)
q <- ggplot(pccor, aes(x=eval(parse(text=paste0("pccor$PC",pcno))))) +
  geom_histogram(colour="grey50", fill="grey50",binwidth=0.01) +
  geom_vline(aes(xintercept=mean(get(paste0("PC",pcno)), na.rm=T)),color="black", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",2))),color="#fdb863", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",3))),color="#8073ac", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",4))),color="#542788", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",5))),color="#006837", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",6))),color="#66bd63", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",7))),color="#1a9850", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",1))),color="#e66101", linetype="dashed", size=2) +
  ylab("") +
  xlab("") +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-0.8,0.8), breaks=c(-0.5,0,0.5), labels=c(-0.5,0,0.5)) +
  scale_y_continuous(expand = c(0.001, 0), limits=c(0,930), breaks=c(0,400,800), labels=c(0,400,800)) +
  theme1 +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.x = element_text(size = 20, color="black"))
print(q)
dev.off()

# -- PC4
pcno=4
GOI1="^SST$"
GOI2="^PVALB$"
GOI3="^GABRA3$"
GOI4="^GABRA5$"
GOI5="^HTR2A$"
GOI6="^HTR7$" #pos
GOI7="^HTR3B$" #neg

# -- Compute percentiles
pc <- paste0("pccor$PC",pcno)
for (gno in 1:7){
  assign(paste0("PC",pcno,"GOIr", gno), c(eval(parse(text=pc)))[grep(get(paste0("GOI", gno)),pccor$Gene)])
  print(sum(c(eval(parse(text=pc))) <= get(paste0("PC",pcno,"GOIr", gno)), na.rm=TRUE) / nrow(pccor))
}

# -- Plot distributions with GOI highlighted
pdf("GeneExpressionCorrelation_PC4.pdf", height=4, width=12)
q <- ggplot(pccor, aes(x=eval(parse(text=paste0("pccor$PC",pcno))))) +
  geom_histogram(colour="grey50", fill="grey50",binwidth=0.01) +
  geom_vline(aes(xintercept=mean(get(paste0("PC",pcno)), na.rm=T)),color="black", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",2))),color="#fdb863", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",3))),color="#8073ac", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",4))),color="#542788", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",5))),color="#006837", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",6))),color="#66bd63", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",7))),color="#1a9850", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",1))),color="#e66101", linetype="dashed", size=2) +
  ylab("") +
  xlab("") +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-0.8,0.8), breaks=c(-0.5,0,0.5), labels=c(-0.5,0,0.5)) +
  scale_y_continuous(expand = c(0.001, 0), limits=c(0,930), breaks=c(0,400,800), labels=c(0,400,800)) +
  theme1 +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.x = element_text(size = 20, color="black"))
print(q)
dev.off()

# -- PC5
pcno=5
GOI1="^SST$"
GOI2="^PVALB$"
GOI3="^GABRA1$"
GOI4="^GABRA2$"
GOI5="^HTR2A$"
GOI6="^HTR5A$" #pos
GOI7="^HTR2C$" #neg

# -- Compute percentiles
pc <- paste0("pccor$PC",pcno)
for (gno in 1:7){
  assign(paste0("PC",pcno,"GOIr", gno), c(eval(parse(text=pc)))[grep(get(paste0("GOI", gno)),pccor$Gene)])
  print(sum(c(eval(parse(text=pc))) <= get(paste0("PC",pcno,"GOIr", gno)), na.rm=TRUE) / nrow(pccor))
}

# -- Plot distributions with GOI highlighted
pdf("GeneExpressionCorrelation_PC5.pdf", height=4, width=12)
q <- ggplot(pccor, aes(x=eval(parse(text=paste0("pccor$PC",pcno))))) +
  geom_histogram(colour="grey50", fill="grey50",binwidth=0.01) +
  geom_vline(aes(xintercept=mean(get(paste0("PC",pcno)), na.rm=T)),color="black", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",2))),color="#fdb863", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",3))),color="#8073ac", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",4))),color="#542788", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",5))),color="#006837", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",6))),color="#66bd63", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",7))),color="#1a9850", linetype="dashed", size=2) +
  geom_vline(aes(xintercept=get(paste0("PC",pcno,"GOIr",1))),color="#e66101", linetype="dashed", size=2) +
  ylab("") +
  xlab("") +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-0.8,0.8), breaks=c(-0.5,0,0.5), labels=c(-0.5,0,0.5)) +
  scale_y_continuous(expand = c(0.001, 0), limits=c(0,930), breaks=c(0,400,800), labels=c(0,400,800)) +
  theme1 +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.x = element_text(size = 20, color="black"))
print(q)
dev.off()