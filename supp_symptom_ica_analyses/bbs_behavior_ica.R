rm(list = ls())
setwd("bbs_manuscript/supp_symptom_ica_analyses")
library("ica")
library("Hmisc")
library("gplots")
library("mgcv")
library("scales")

# -- Functions
theme1<- theme(panel.border = element_rect(size = 2, colour = 'NA', fill='NA'), 
               axis.line = element_line(colour = "black", size = 2),
               axis.text.x = element_text(size = 25, color="black"), 
               axis.title.x = element_text(size = 30,margin=margin(15,0,0,0)),
               axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
               axis.title.y = element_text(size = 30,margin=margin(0,15,0,0)),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               plot.background = element_rect(fill = "transparent", colour = NA),
               panel.background = element_rect(fill = "white", colour = NA),
               legend.position="none",
               axis.ticks = element_line(colour = "black", size=2),
               axis.ticks.length = unit(.4, "cm"))

# -- Read data
beh_raw = read.table("../data_symptoms/BSNIP_PSD_N436_Behavior.dat", sep="\t", header=TRUE)
beh=beh_raw[1:436,2:37]
icaA = icafast(beh,nc=5)
beh_raw$Group <- as.factor(beh_raw$Group)

pcaA = prcomp(beh, scale. = TRUE)

# -- Correlation matrix between PCA and ICA component scores
pcaica5 <- data.frame(cbind(pcaA$x[,1:5],icaA$S))
names(pcaica5)[6:10] <- c("IC1","IC2","IC3","IC4","IC5")
pcaica5r <- rcorr(as.matrix(pcaica5))
pcaica5rplot <- pcaica5r$r[6:10,1:5]
colours=colorRampPalette(c("blue","white","red"))(99)
heatmap.2(pcaica5rplot, 
          Rowv=FALSE, 
          Colv=FALSE,
          dendrogram="none", 
          trace="none", 
          density.info="none", 
          srtCol=0,
          offsetCol=1, 
          col=colours,
          breaks=seq(-1,1,length.out=100),
          cexRow = 2,
          cexCol = 2,
          cellnote=format(round(pcaica5rplot, 2), nsmall = 2),
          notecol="black",
          notecex = 3,
          adjCol=c(NA,1),
          key.par=list(cex.axis=1.5,mar=c(5.4, 0, 0.5, 5)),
          key.xlab=NA,
          key.ylab=NA,
          key.title = NA,
          lmat=rbind(c(7,8,9), c(5, 1, 2),  c(6, 4, 3)), lhei=c(0.1, 7, 0.1), lwid=c(0.1,7, 0.1))
dev.off()

# -- Plot screeplot of variance explained
prop_var5 <- data.frame("IC"=1:5, "Variance"=icaA$vafs)
ggplot(data=prop_var5, aes(x=IC, y=Variance)) +
  geom_line(lwd=1.2)+
  geom_point(size=4,col="black",pch=21,bg="grey35") +
  geom_point(data=prop_var5, aes(size=Variance),col="black",pch=21,bg="green") +
  scale_size_area(name=prop_var5, max_size=15) +
  ylab('% Variance Accounted') +
  xlab('Independent\nComponent') +
  scale_x_continuous(limits=c(0.5,5.2),labels=1:5, breaks=1:5) +
  scale_y_continuous(limits=c(0,0.2),labels=c(0,5,10,15,20), breaks=c(0,0.05,0.10,0.15,0.20)) +
  theme1  +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
       axis.title.y = element_text(size = 30),
       axis.text.y =  element_text(size = 30),
       axis.title.x = element_text(size = 30,margin=margin(10,0,0,0)),
       axis.text.x =  element_text(size = 30, margin=margin(10,10,0,0)))
dev.off()

# -- Plot piechart variance explained
pveIC_prop <- data.frame(ICs=c("1-5","Unaccounted"),"Total Variance"=c(sum(prop_var5$Variance[1:5]),1-sum(prop_var5$Variance[1:5])))
ggplot(data=pveIC_prop,aes(x=factor(1),y=Total.Variance,fill=factor(ICs))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta="y") +
  scale_fill_manual(values = c("1-5"="green","Unaccounted"="grey35")) +
  geom_text(aes(y = c(0.75,0.22), label = percent(Total.Variance)), size=11, col=c("black","white")) +
  theme1 +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line=element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA))
dev.off()