rm(list = ls())
setwd("../fig2_symptom_pca_stability_analyses")
library(ggplot2)
library(nFactors)
library(ggbiplot)
library(scales)
library(ggridges)
library(plot3D)
library(rgl)
library(magick)
library(caret)
library(dmm)
library(plyr)

# -- Themes
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

# -- Functions
sign.pc<-function(x,R=5000,s=36,...){
  # run PCA
  pc.out<-prcomp(x, scale=TRUE,...)
  # the proportion of variance of each PC
  pve=(pc.out$sdev^2/sum(pc.out$sdev^2))[1:s]
  # a matrix with R rows and s columns that contains the proportion of variance explained by each pc for each randomization replicate.
  pve.perm<-matrix(NA,ncol=s,nrow=R)
  for(i in 1:R){
    # permutation each column
    x.perm<-apply(x,2,sample)
    # run PCA
    pc.perm.out<-prcomp(x.perm,scale=TRUE,...)
    # the proportion of variance of each PC.perm
    pve.perm[i,]=(pc.perm.out$sdev^2/sum(pc.perm.out$sdev^2))[1:s]
  }
  # calcalute the p-values
  pval<-apply(t(pve.perm)>pve,1,sum)/R
  return(list(pve=pve,pval=pval, pve.perm=pve.perm))
}

# -- leave-site-out validation
site.loo <- function(behS,R=5000,s=36,...){
  pcaS = prcomp(behS, scale. = TRUE)
  loadings <- pcaS$rotation
  loadingsm <- as.matrix(loadings)
  prop_var <- data.frame(cbind((1:s), pcaS$sdev^2/sum(pcaS$sdev^2)))
  names(prop_var) <- c("PC", "Variance")
  signif <- sign.pc(behS,R=R,s=s)
  pve.perm <- signif$pve.perm
  pve <- signif$pve
  pval <- signif$pval
  crit95 <- matrix(0,nrow=1,ncol=s) # Compute 0.95 chance
  crit95i <- ceiling(0.05 * R)
  for(i in 1:s){
    crit95[i] <- sort(pve.perm[1:R, i])[1]
  }
  crit=data.frame(PC=c(1:s),crit=c(crit95))
  return(list(pve=pve,pval=pval, pve.perm=pve.perm, crit=crit, prop_var=prop_var, loadings=loadingsm))
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

###############################
##### LOO Site Validation #####
###############################
# -- Read data -- get site data from BSNIP
beh = read.table("../data_symptoms/BSNIP_PSD_N436_Behavior.dat", sep="\t", header=TRUE)
beh_Site<-cbind(Site=c("GP","GT","CT","GP","JS","JS","GP","JS","JS","MK","GT","GP","CT","GT","GT","MK","MK","GT","GP","GT","GT","CT","GT","GP","GT","CT","CT","GT","JS","CT","JS","GP","CT","MK","JS","JS","GP","CT","MK","CT","MK","JS","JS","CT","GT","CT","CT","CT","GP","GP","GT","GT","GT","GP","MK","GT","JS","GP","MK","GP","CT","JS","GT","CT","GP","JS","GP","GT","CT","JS","JS","GT","JS","GP","JS","CT","CT","CT","JS","GP","MK","GT","GP","GP","JS","CT","GT","JS","JS","CT","GP","MK","GP","CT","GP","GT","GP","GT","CT","GP","GT","GT","CT","MK","GP","GT","JS","GP","GT","GT","CT","GP","JS","GP","GP","GT","GT","GT","GT","CT","GT","JS","JS","JS","GT","JS","JS","MK","JS","CT","GT","CT","GP","GT","CT","GT","GP","MK","GP","GP","CT","GT","GP","GP","JS","GP","MB","MK","CT","GT","GP","JS","GT","CT","GT","MB","GT","JS","CT","GP","GT","GT","JS","GT","JS","GP","GP","GT","GT","GT","GP","GP","GP","GT","MK","CT","GT","CT","JS","GP","JS","JS","JS","JS","GT","GP","GP","CT","MB","GP","GP","JS","MB","GP","JS","CT","JS","GP","JS","GT","GP","MK","GT","CT","CT","GT","JS","GT","GP","JS","CT","JS","JS","GT","GT","JS","GT","GP","GT","GT","JS","MK","CT","JS","CT","CT","GT","MB","MK","JS","MK","MK","CT","GP","CT","GT","CT","GP","CT","MK","JS","GP","GT","GP","JS","GP","GP","JS","GP","GT","GP","CT","GT","GT","GT","GP","GP","CT","GP","MK","MK","GT","MK","CT","CT","MK","GT","JS","JS","GT","CT","MK","GT","GT","GT","MK","GP","JS","JS","GP","JS","GP","GT","GP","GP","MK","JS","MK","GT","JS","GT","GT","JS","MB","GP","GT","GT","GT","MK","MK","GT","JS","CT","GT","JS","JS","MK","MK","GT","JS","GP","GT","JS","JS","GP","GP","MK","GT","CT","GT","GT","JS","JS","GP","JS","CT","JS","JS","CT","JS","GP","GT","CT","GP","CT","JS","CT","CT","GP","MK","GT","JS","MB","GT","GP","CT","CT","CT","CT","CT","JS","GT","GP","GT","JS","MK","CT","GT","GP","CT","MK","GT","CT","CT","GP","JS","JS","GP","GT","JS","GP","GP","GP","GT","GT","MK","MB","GT","GP","CT","GP","CT","JS","CT","GT","JS","MK","GP","CT","JS","CT","CT","GT","MK","JS","JS","GP","GT","JS","JS","MK","GT","GT","MB","GP","MK","JS","JS","JS","GT","JS","CT","GP","MB","GP","CT","GP","GP","GT","GP","JS","GP","GT","GT","GT","CT","MK","JS","GP","JS","JS","CT","MK","MK","MK","JS"),beh)
beh_Site$Group <- as.factor(beh_Site$Group)
beh_Site$Site  <- as.factor(beh_Site$Site)

# -- Run PCA with subjects from all sites except one
R=5000
s=36
for (site in c("CT", "GP", "GT", "JS", "MB", "MK")){
  nam <- paste("LOO_", site, sep = "")
  behS <- subset(beh_Site, Site != site) # select all subs not in one site
  assign(nam, site.loo(behS[,-1:-2], R=R, s=s))
}

# -- Compute proportion variance explained for each LOO PCA
for (site in c("CT", "GP", "GT", "JS", "MB", "MK")){
  nam <- paste("pve_prop", site, sep = "")
  assign(nam, data.frame(PCs=c("1-5","6-36"),"Total Variance"=c(sum(get(paste0("LOO_", site))$pve[1:5]),sum(get(paste0("LOO_", site))$pve[6:36]))))
}

# -- Plot the screeplots
for (site in c("CT", "GP", "GT", "JS", "MB", "MK")){
  pdf(paste0("Symptom_PCA_LOSO_",site,"_Screeplot.pdf"), height=6, width=10)
  q <- ggplot(data=get(paste0("LOO_", site))$prop_var, aes(x=PC, y=Variance)) +
    geom_line(data=get(paste0("LOO_", site))$crit, aes(x=PC,y=crit), lty=2, col="dark red",lwd=1.2) +
    geom_line(lwd=1.2)+
    geom_point(size=4,col="black",pch=21,bg="grey35") +
    geom_point(data=get(paste0("LOO_", site))$prop_var[1:5,1:2], aes(size=1/PC),col="black",pch=21,bg="green") +
    scale_size_area(name=get(paste0("LOO_", site))$prop_var[1:5,1:2], max_size=15) +
    ylab('% Variance Explained') +
    xlab('Principal Component') +
    scale_x_continuous(labels=c(0,36), breaks=c(0,36)) +
    scale_y_continuous(labels=c(0,5,10,15,20), breaks=c(0,0.05,0.10,0.15,0.20)) +
    theme1 + 
    theme(panel.background = element_rect(fill = "transparent", colour = NA),
          axis.title.y = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.x = element_text(size = 35,margin=margin(-20,-30,0,0)),
          axis.text.x = element_text(size = 30, margin=margin(10,0,0,0)))
  print(q)
  dev.off()
}

# -- Pie charts of proportion variance explained by significant vs. non-sig. PCs
for (site in c("CT", "GP", "GT", "JS", "MB", "MK")){
  pdf(paste0("Symptom_PCA_LOSO_",site,"_Piechart.pdf"), height=6, width=10)
  q <- ggplot(data=get(paste0("pve_prop", site)),aes(x=factor(1),y=Total.Variance,fill=factor(PCs))) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar(theta="y") +
    scale_fill_manual(values = c("1-5"="green","6-36"="grey35")) +
    geom_text(aes(y = Total.Variance/2 + c(cumsum(Total.Variance)[-length(Total.Variance)],0), 
                  label = percent(as.numeric(format(round(get(paste0("pve_prop", site))$Total.Variance, 4),nsmall=4)))), size=8, col=c("black","white")) +
    theme1 +
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line=element_blank(),
          legend.position="none",
          panel.background = element_rect(fill = "transparent", colour = NA)) +
    guides(fill=guide_legend(keywidth=0.7,keyheight=0.7, default.unit="inch",title="PCs"))
  print(q)
  dev.off()
}

# -- Computing the PC scores on holdout subjects in the left-out site
pcaA <- prcomp(beh_Site[-1:-2], scale. = TRUE)
scoresA <- cbind(beh_Site[,1:2], pcaA$x)
pcVar <- pcaA$x %*% pcaA$rotation

# -- Normalize by site mean and sd
for (site in c("CT", "GP", "GT", "JS", "MB", "MK")){
  nam <- paste0("beh", site)
  assign(nam, as.matrix(scale(subset(beh_Site, Site == site, select= -(Group:Site)))[,])) # columns are Sx, rows are subjects
  assign(paste0("loadings", site), as.matrix(get(paste0("LOO_", site))$loadings))
  assign(paste0("pred", site), as.data.frame(as.matrix(get(nam)) %*% as.matrix(get(paste0("loadings", site)))))
  assign(paste0("obs", site), subset(scoresA, Site == site, select= -(Group:Site)))
}

# -- Scatterplots showing correlation between computed and observed PC1 scores, for each site
cols <- subset(beh_Site, Site == "CT", select=Group)
pdf(paste0("Symptom_PCA_LOSO_Predicted_Site1.pdf"), height=6, width=10)
ggplot(predCT, aes(x=PC1, y=obsCT$PC1, col=cols$Group)) +
  stat_smooth(method = "lm", col = "dark blue", lwd=2) +
  geom_point(size=5, alpha=0.7) +
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_continuous(limits=c(-7.5,7.5),breaks=seq(-8,8,16),labels=seq(-8,8,16)) +
  scale_x_continuous(limits=c(-7.5,7.5),breaks=seq(-8,8,16),labels=seq(-8,8,16)) +
  ylab("Observed\nPC1 Scores") +
  xlab("Predicted\nPC1 Scores") +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme1 + 
  theme(axis.text.x = element_text(size = 25, color="black"), 
    axis.title.x = element_text(size = 25, margin=margin(-10,0,0,0)),
    axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
    axis.title.y = element_text(size = 25,margin=margin(0,-15,0,0)))
dev.off()

cols <- subset(beh_Site, Site == "GP", select=Group)
pdf(paste0("Symptom_PCA_LOSO_Predicted_Site2.pdf"), height=6, width=10)
ggplot(predGP, aes(x=PC1, y=obsGP$PC1, col=cols$Group)) +
  stat_smooth(method = "lm", col = "dark blue", lwd=2) +
  geom_point(size=5, alpha=0.7)+
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_continuous(limits=c(-10,5),breaks=seq(-10,5,15),labels=seq(-10,5,15)) +
  scale_x_continuous(limits=c(-10,5),breaks=seq(-10,5,15),labels=seq(-10,5,15)) +
  ylab("Observed\nPC1 Scores") +
  xlab("Predicted\nPC1 Scores") +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme1 + 
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(-10,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,-15,0,0)))
dev.off()

cols <- subset(beh_Site, Site == "GT", select=Group)
pdf(paste0("Symptom_PCA_LOSO_Predicted_Site3.pdf"), height=6, width=10)
ggplot(predGT, aes(x=PC1, y=obsGT$PC1, col=cols$Group)) +
  stat_smooth(method = "lm", col = "dark blue", lwd=2) +
  geom_point(size=5, alpha=0.7)+
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_continuous(limits=c(-8,6),breaks=seq(-8,6,14),labels=seq(-8,6,14)) +
  scale_x_continuous(limits=c(-8,6),breaks=seq(-8,6,14),labels=seq(-8,6,14)) +
  ylab("Observed\nPC1 Scores") +
  xlab("Predicted\nPC1 Scores") +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme1 + 
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(-10,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,-15,0,0)))
dev.off()

cols <- subset(beh_Site, Site == "JS", select=Group)
pdf(paste0("Symptom_PCA_LOSO_Predicted_Site4.pdf"), height=6, width=10)
ggplot(predJS, aes(x=PC1, y=obsJS$PC1, col=cols$Group)) +
  stat_smooth(method = "lm", col = "dark blue", lwd=2) +
  geom_point(size=5, alpha=0.7)+
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_continuous(limits=c(-8,6),breaks=seq(-8,6,14),labels=seq(-8,6,14)) +
  scale_x_continuous(limits=c(-8,6),breaks=seq(-8,6,14),labels=seq(-8,6,14)) +
  ylab("Observed\nPC1 Scores") +
  xlab("Predicted\nPC1 Scores") +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme1 + 
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(-10,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,-15,0,0)))
dev.off()

cols <- subset(beh_Site, Site == "MK", select=Group)
pdf(paste0("Symptom_PCA_LOSO_Predicted_Site5.pdf"), height=6, width=10)
ggplot(predMK, aes(x=PC1, y=obsMK$PC1, col=cols$Group)) +
  stat_smooth(method = "lm", col = "dark blue", lwd=2) +
  geom_point(size=5, alpha=0.7)+
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_continuous(limits=c(-7,7),breaks=seq(-7,7,14),labels=seq(-7,7,14)) +
  scale_x_continuous(limits=c(-7,7),breaks=seq(-7,7,14),labels=seq(-7,7,14)) +
  ylab("Observed\nPC1 Scores") +
  xlab("Predicted\nPC1 Scores") +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme1 + 
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(-10,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,-15,0,0)))
dev.off()

cols <- subset(beh_Site, Site == "MB", select=Group)
pdf(paste0("Symptom_PCA_LOSO_Predicted_Site6.pdf"), height=6, width=10)
ggplot(predMB, aes(x=PC1, y=obsMB$PC1, col=cols$Group)) +
  stat_smooth(method = "lm", col = "dark blue", lwd=2) +
  geom_point(size=5, alpha=0.7)+
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_continuous(limits=c(-5,5),breaks=seq(-5,5,10),labels=seq(-5,5,10)) +
  scale_x_continuous(limits=c(-5,5),breaks=seq(-5,5,10),labels=seq(-5,5,10)) +
  ylab("Observed\nPC1 Scores") +
  xlab("Predicted\nPC1 Scores") +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme1 + 
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(-10,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,-15,0,0)))
dev.off()

# -- Compute the correlation values
for (site in c("CT", "GP", "GT", "JS", "MB", "MK")){
  for (pc in c(1:5)){
    print(paste0("pred",site, pc))
    print(cor(get(paste0("pred",site))[,pc],get(paste0("obs",site))[,pc]))
  }
}

# -- Summarize 
propvarsumm <- data.frame(cbind(c("CT","GP","GT","JS","MB","MK"),5,rbind(pve_propCT$Total.Variance[1],pve_propGP$Total.Variance[1],pve_propGT$Total.Variance[1],pve_propJS$Total.Variance[1],pve_propMB$Total.Variance[1],pve_propMK$Total.Variance[1])))
names(propvarsumm) <- c("Site","PCSig","PropVar")
propvarsumm$PropVar <- unfactor(propvarsumm$PropVar)
propvarsumm$PCSig <- unfactor(propvarsumm$PCSig)

# -- Plot summary of no. sig. PCs and total var explained
pdf(paste0("Symptom_PCA_LOSO_PropVarExplained.pdf"), height=6, width=6)
ggplot(propvarsumm, aes(x=Site, y=PropVar)) +
  geom_col(fill="white", col="black", lwd=1) +
  theme1 +
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 30, margin=margin(15,0,0,0)),
        axis.text.y = element_text(size = 25, color="black", angle=90, hjust=0.5, margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,10,0,0)),
        axis.title.y.right = element_text(size = 25, color="green3", angle=270),
        #axis.line.y.right = element_line(color = "red"),
        axis.text.y.right = element_text(size = 25, color="green3", angle=270, hjust=0.5, margin=margin(0,5,0,5))) +
  geom_line(aes(y=PCSig/12,group = 1), col="green3", lwd=2) +
  geom_point(aes(y=PCSig/12), col="green3",size=6) +
  ylab("Total % Var Explained") +
  xlab("") +
  scale_x_discrete(labels=seq(1,6,1)) +
  scale_y_continuous(expand = c(0.001, 0), limits=c(0,0.55), breaks=seq(0,0.5,0.1),labels=seq(0,50,10),
                     sec.axis = sec_axis(~.*12, name = "No. Significant PCs", breaks=seq(0,6,1),labels=seq(0,6,1)))
dev.off()

###############################
###### 5-fold Validation ######
###############################
# -- Break up the data into k folds
k=5
R=5000
s=36
flds <- createFolds(1:nrow(beh_Site), k = k, list = TRUE, returnTrain = FALSE)

# -- Compute PCAs, each time leaving out a fold
for (fold in 1:k){
  tnam <- paste("t",fold, sep = "")
  nam <- paste("CV_t",fold, sep = "")
  assign(tnam, beh_Site[-c(eval(parse(text=paste0("flds$Fold",fold)))),])
  assign(nam, site.loo(get(tnam)[,-1:-2], R=R, s=s))
}

# -- Compute prop variance explained, each time leaving out a fold
for (fold in 1:k){
  nam <- paste0("pve_prop", fold)
  assign(nam, data.frame(PCs=c("1-5","6-36"),"Total Variance"=c(sum(get(paste0("CV_t", fold))$pve[1:5]),sum(get(paste0("CV_t",fold))$pve[6:36]))))
}

# -- Plot the screeplots
for (fold in 1:5){
  pdf(paste0("Symptom_PCA_kfold_CV_Fold",fold,"_Screeplot.pdf"), height=6, width=10)
  q <- ggplot(data=get(paste0("CV_t", fold))$prop_var, aes(x=PC, y=Variance)) +
    geom_line(data=get(paste0("CV_t", fold))$crit, aes(x=PC,y=crit), lty=2, col="dark red",lwd=1.2) +
    geom_line(lwd=1.2)+
    geom_point(size=4,col="black",pch=21,bg="grey35") +
    geom_point(data=get(paste0("CV_t", fold))$prop_var[1:5,1:2], aes(size=1/PC),col="black",pch=21,bg="green") +
    scale_size_area(name=get(paste0("CV_t", fold))$prop_var[1:5,1:2], max_size=15) +
    ylab('% Variance Explained') +
    xlab('Principal Component') +
    scale_x_continuous(labels=c(0,36), breaks=c(0,36)) +
    scale_y_continuous(labels=c(0,5,10,15,20), breaks=c(0,0.05,0.10,0.15,0.20)) +
    theme1 + 
    theme(panel.background = element_rect(fill = "transparent", colour = NA),
          axis.title.y = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.x = element_text(size = 35,margin=margin(-20,-30,0,0)),
          axis.text.x = element_text(size = 30, margin=margin(10,0,0,0)))
  print(q)
  dev.off()
}

# -- Plot the pie charts
for (fold in 1:5){
  pdf(paste0("Symptom_PCA_kfold_CV_Fold",fold,"_Piechart.pdf"), height=6, width=10)
  q <- ggplot(data=get(paste0("pve_prop", fold)),aes(x=factor(1),y=Total.Variance,fill=factor(PCs))) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar(theta="y") +
    scale_fill_manual(values = c("1-5"="green","6-36"="grey35")) +
    geom_text(aes(y = Total.Variance/2 + c(cumsum(Total.Variance)[-length(Total.Variance)],0), 
                  label = percent(as.numeric(format(round(get(paste0("pve_prop", fold))$Total.Variance, 4),nsmall=4)))), size=8, col=c("black","white")) +
    theme1 +
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line=element_blank(),
          legend.position="none",
          panel.background = element_rect(fill = "transparent", colour = NA)) +
    guides(fill=guide_legend(keywidth=0.7,keyheight=0.7, default.unit="inch",title="PCs"))
  print(q)
  dev.off()
}

# -- Compute PC scores for the left out subjects
pcaA <- prcomp(beh_Site[-1:-2], scale. = TRUE)
scoresA <- cbind(beh_Site[,1:2], pcaA$x)
pcVar <- pcaA$x %*% pcaA$rotation

# -- Normalize by mean and sd
for (fold in 1:k){
  namm <- paste0("CVtest", fold)
  assign(namm, as.matrix(scale(beh_Site[c(eval(parse(text=paste0("flds$Fold",fold)))),][-1:-2]))[,]) # columns are Sx, rows are subjects
  assign(paste0("loadings_train", fold), as.matrix(get(paste0("CV_t", fold))$loadings))
  assign(paste0("pred_test", fold), as.data.frame(as.matrix(get(namm)) %*% as.matrix(get(paste0("loadings_train", fold)))))
  assign(paste0("obs_test", fold), scoresA[c(eval(parse(text=paste0("flds$Fold",fold)))),-1:-2])
}

# -- Plot the scatterplots of pred vs obs
pdf(paste0("Symptom_PCA_PredvsObs_Fold1.pdf"), height=6, width=10)
cols <- beh_Site[c(flds$Fold1),]$Group
ggplot(pred_test1, aes(x=PC1, y=obs_test1$PC1, col=cols)) +
  stat_smooth(method = "lm", col = "dark blue", lwd=2) +
  geom_point(size=5, alpha=0.7) +
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_continuous(limits=c(-8,6),breaks=seq(-6,6,3),labels=seq(-6,6,3)) +
  scale_x_continuous(limits=c(-8,6),breaks=seq(-6,6,3),labels=seq(-6,6,3)) +
  ylab("Observed") +
  xlab("Predicted") +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme1 + 
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 30, margin=margin(15,0,0,0)),
        axis.text.y = element_text(size = 25, color="black", angle=90, hjust=0.5, margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 30,margin=margin(0,10,0,0)),
        axis.title.y.right = element_text(size = 25, color="green3", angle=270),
        axis.text.y.right = element_text(size = 25, color="green3", angle=270, hjust=0.5, margin=margin(0,5,0,5)))
dev.off()

pdf(paste0("Symptom_PCA_PredvsObs_Fold2.pdf"), height=6, width=10)
cols <- beh_Site[c(flds$Fold2),]$Group
ggplot(pred_test2, aes(x=PC1, y=obs_test2$PC1, col=cols)) +
  stat_smooth(method = "lm", col = "dark blue", lwd=2) +
  geom_point(size=5, alpha=0.7) +
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_continuous(limits=c(-7,6),breaks=seq(-7,6,13),labels=seq(-7,6,13)) +
  scale_x_continuous(limits=c(-7,6),breaks=seq(-7,6,13),labels=seq(-7,6,13)) +
  ylab("Observed\nPC1 Scores") +
  xlab("Predicted\nPC1 Scores") +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme1 + 
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(-10,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,-15,0,0)))
dev.off()

pdf(paste0("Symptom_PCA_PredvsObs_Fold3.pdf"), height=6, width=10)
cols <- beh_Site[c(flds$Fold3),]$Group
ggplot(pred_test3, aes(x=PC1, y=obs_test3$PC1, col=cols)) +
  stat_smooth(method = "lm", col = "dark blue", lwd=2) +
  geom_point(size=5, alpha=0.7) +
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_continuous(limits=c(-6.2,6),breaks=seq(-6,6,12),labels=seq(-6,6,12)) +
  scale_x_continuous(limits=c(-6.2,6),breaks=seq(-6,6,12),labels=seq(-6,6,12)) +
  ylab("Observed\nPC1 Scores") +
  xlab("Predicted\nPC1 Scores") +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme1 + 
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(-10,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,-15,0,0)))
dev.off()

pdf(paste0("Symptom_PCA_PredvsObs_Fold4.pdf"), height=6, width=10)
cols <- beh_Site[c(flds$Fold4),]$Group
ggplot(pred_test4, aes(x=PC1, y=obs_test4$PC1, col=cols)) +
  stat_smooth(method = "lm", col = "dark blue", lwd=2) +
  geom_point(size=5, alpha=0.7) +
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_continuous(limits=c(-10,6),breaks=seq(-10,6,16),labels=seq(-10,6,16)) +
  scale_x_continuous(limits=c(-10,6),breaks=seq(-10,6,16),labels=seq(-10,6,16)) +
  ylab("Observed\nPC1 Scores") +
  xlab("Predicted\nPC1 Scores") +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme1 + 
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(-10,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,-15,0,0)))
dev.off()

pdf(paste0("Symptom_PCA_PredvsObs_Fold5.pdf"), height=6, width=10)
cols <- beh_Site[c(flds$Fold5),]$Group
ggplot(pred_test5, aes(x=PC1, y=obs_test5$PC1, col=cols)) +
  stat_smooth(method = "lm", col = "dark blue", lwd=2) +
  geom_point(size=5, alpha=0.7) +
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_continuous(limits=c(-7.2,5),breaks=seq(-7,5,12),labels=seq(-7,5,12)) +
  scale_x_continuous(limits=c(-7.2,5),breaks=seq(-7,5,12),labels=seq(-7,5,12)) +
  ylab("Observed\nPC1 Scores") +
  xlab("Predicted\nPC1 Scores") +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme1 + 
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(-10,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,-15,0,0)))
dev.off()

# -- Compute the correlation values
for (fold in 1:5){
  for (pc in c(1:5)){
    print(paste0("pred_test",fold,"_PC", pc))
    print(cor(get(paste0("pred_test",fold))[,pc],get(paste0("obs_test",fold))[,pc]))
  }
}

# -- Summarize 
propvarsummkfold <- data.frame(cbind(1:5,5,rbind(pve_prop1$Total.Variance[1],pve_prop2$Total.Variance[1],pve_prop3$Total.Variance[1],pve_prop4$Total.Variance[1],pve_prop5$Total.Variance[1])))
names(propvarsummkfold) <- c("Fold","PCSig","PropVar")

# -- Plot summary of no. sig. PCs and total var explained
pdf(paste0("Symptom_PCA_Summary_KFoldCV.pdf"), height=6, width=10)
ggplot(propvarsummkfold, aes(x=Fold, y=PropVar)) +
  geom_col(fill="white", col="black",lwd=1) +
  theme1 +
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 30, margin=margin(15,0,0,0)),
        axis.text.y = element_text(size = 25, color="black", angle=90, hjust=0.5, margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,10,0,0)),
        axis.title.y.right = element_text(size = 25, color="green3", angle=270),
        #axis.line.y.right = element_line(color = "green3"),
        axis.text.y.right = element_text(size = 25, color="green3", angle=270, hjust=0.5, margin=margin(0,5,0,5))) +
  geom_line(aes(y=PCSig/12,group = 1), col="green3", lwd=2) +
  geom_point(aes(y=PCSig/12), col="green3",size=6) +
  ylab("Total % Var Explained") +
  xlab("") +
  scale_x_continuous(labels=seq(1,5,1),breaks=seq(1,5,1)) +
  scale_y_continuous(expand = c(0.001, 0), limits=c(0,0.55), breaks=seq(0,0.5,0.1),labels=seq(0,50,10),
                     sec.axis = sec_axis(~.*12, name = "No. Significant PCs", breaks=seq(0,6,1),labels=seq(0,6,1)))
dev.off()

###############################
###### k-fold validation ######
###############################
# -- Split into k folds, where k = 2 to 20
k=2:20
R=10
s=36
for (kfold in k){
  assign(paste0("k",kfold,"fold"), createFolds(1:nrow(beh_Site), k = kfold, list = TRUE, returnTrain = FALSE))
}

# -- Compute the PCAs and prop var explained for each value of k
for (kfold in k){
  knames=NULL
  for (fold in 1:kfold){
    knames <- c(knames, paste0("Fold",fold))
  }
  x <- get(paste0("k",kfold,"fold"))
  names(x) <- knames
  assign(paste0("k",kfold,"fold"), x)
}
for (kfold in k){
  for (fold in 1:kfold){
    tnam <- paste("k",kfold,"beh_train",fold, sep = "")
    nam <- paste("k",kfold,"CV_train",fold, sep = "")
    assign(tnam, beh_Site[-c(eval(parse(text=paste0("k",kfold,"fold$Fold",fold)))),])
    assign(nam, site.loo(get(tnam)[,-1:-2], R=R, s=s))
  }
}
for (kfold in k){
  for (fold in 1:kfold){
  nam <- paste0("k",kfold,"pve_prop", fold)
  assign(nam, data.frame(PCs=c("1-5","6-36"),"Total Variance"=c(sum(get(paste0("k",kfold,"CV_train", fold))$pve[1:5]),sum(get(paste0("k",kfold,"CV_train",fold))$pve[6:36]))))
  }
}

## Predicting left out subjects, normalize by site mean and sd
pcaA <- prcomp(beh_Site[-1:-2], scale. = TRUE)
scoresA <- cbind(beh_Site[,1:2], pcaA$x)
pcVar <- pcaA$x %*% pcaA$rotation
for (kfold in k){
  for (fold in 1:kfold){
    namm <- paste0("k",kfold,"CV_test", fold)
    assign(namm, as.matrix(scale(beh_Site[c(eval(parse(text=paste0("k",kfold,"fold$Fold",fold)))),][-1:-2]))[,]) # columns are Sx, rows are subjects
    assign(paste0("k",kfold,"loadings_train", fold), as.matrix(get(paste0("k",kfold,"CV_train", fold))$loadings))
    assign(paste0("k",kfold,"pred_test", fold), as.data.frame(as.matrix(get(namm)) %*% as.matrix(get(paste0("k",kfold,"loadings_train", fold)))))
    assign(paste0("k",kfold,"obs_test", fold), scoresA[c(eval(parse(text=paste0("k",kfold,"fold$Fold",fold)))),-1:-2])
  }
}

# -- Compute the correlation between pred and obs
for (kfold in k){
  assign(paste0("k",kfold,"test"), matrix(data=NA,nrow=5,ncol=kfold))
  x <- matrix(data=NA,nrow=5,ncol=kfold)
  for (fold in c(1:kfold)){
    for (pc in c(1:5)){
      x[pc,fold] <- cor(get(paste0("k",kfold,"pred_test",fold))[,pc],get(paste0("k",kfold,"obs_test",fold))[,pc])
    }
  }
  assign(paste0("k",kfold,"test"), x)
}

# -- Make the data frame with the pertinent information for all values of k, all folds, all PCs, correlation values for each
kcat=NULL
for (k in 2:20){kcat=c(kcat,replicate(5*k,k))}
fcat=NULL
for (f in 2:20){for (c in 1:f){fcat=c(fcat,replicate(5,c))}}
dat<-NULL
for (f in 2:20){
  dat <- c(dat,get(paste0("k",f,"test")))
}
DF <- data.frame("k"=kcat, "Fold"=fcat,"PC"=1:5, "R"=abs(dat))
# -- Compute summary statistics across k values and PCs
DFsumm <- summarySE(DF, measurevar="R", groupvars=c("k","PC"))

pdf(paste0("Symptom_PCA_Summmary_kfoldCV.pdf"), height=6, width=10)
ggplot(DFsumm, aes(x=k, y=R, alpha=as.factor(PC))) +
  geom_hline(yintercept=1, col="red", lwd=2) +
  stat_summary(position="dodge", geom = "bar", fun="mean",fill=c("#485C77")) +
  geom_errorbar(aes(x=k, ymin=R-se, ymax=R+se), width=.5, position=position_dodge(0.9), lwd=1.5) +
  scale_alpha_manual(values = c("1"=1,"2"=0.9,"3"=0.8, "4"=0.7, "5"=0.6)) +
  scale_y_continuous(expand=c(0.0001,0), limits=c(0.0, 1.1),breaks=seq(0,1,0.5),labels=format(seq(0,1,0.5), nsmall=1)) +
  scale_x_continuous(expand=c(0.000001,0),limits=c(1.45,10.5),breaks=seq(2,10,by=1),labels=seq(2,10,1)) +
  ylab("r-value") + 
  xlab("") +
  theme1 +
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 30, margin=margin(15,0,0,0)),
        axis.text.y = element_text(size = 25, color="black", angle=90, hjust=0.5, margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,10,0,0)))
dev.off()
