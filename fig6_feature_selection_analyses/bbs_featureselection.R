rm(list = ls())
setwd("bbs_manuscript/fig6_feature_selection_analyses")

library(ggplot2)
library(emojifont)
library(plyr)

# -- Functions
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

# -- Themes
theme1<- theme(panel.border = element_rect(size = 2, colour = 'NA', fill='NA'), 
               axis.line = element_line(colour = "black", size = 2),
               axis.text.x = element_text(size = 40, color="black"), 
               axis.title.x = element_text(size = 40,margin=margin(15,0,0,0)),
               axis.text.y = element_text(size = 40, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
               axis.title.y = element_text(size = 40,margin=margin(0,15,0,0)),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               plot.background = element_rect(fill = "transparent", colour = NA),
               panel.background = element_rect(fill = "white", colour = NA),
               legend.position="none",
               axis.ticks = element_line(colour = "black", size=2),
               axis.ticks.length = unit(.4, "cm"))     
            
# -- Read data (outputs from featureselection.m)
dpGBC_PredObsv_PC1 <- as.data.frame(read.table("BSNIP_PC1_dpGBCPred_dpGBCObsv_r.txt", header=FALSE)[1:717,])
dpGBCObsv_PC1ScorePred <- as.data.frame(read.table("BSNIP_PC1_dpGBCObs_PC1ScorePred_r.txt", header=FALSE)[1:717,])
names(dpGBCObsv_PC1ScorePred)=c("V1")
names(dpGBC_PredObsv_PC1)=c("V1")

dpGBC_PredObsv_PC2 <- as.data.frame(read.table("BSNIP_PC2_dpGBCPred_dpGBCObsv_r.txt", header=FALSE)[1:717,])
dpGBCObsv_PC2ScorePred <- as.data.frame(read.table("BSNIP_PC2_dpGBCObs_PC2ScorePred_r.txt", header=FALSE)[1:717,])
names(dpGBCObsv_PC2ScorePred)=c("V1")
names(dpGBC_PredObsv_PC2)=c("V1")

dpGBC_PredObsv_PC3 <- as.data.frame(read.table("BSNIP_PC3_dpGBCPred_dpGBCObsv_r.txt", header=FALSE)[1:717,])
dpGBCObsv_PC3ScorePred <- as.data.frame(read.table("BSNIP_PC3_dpGBCObs_PC3ScorePred_r.txt", header=FALSE)[1:717,])
names(dpGBCObsv_PC3ScorePred)=c("V1")
names(dpGBC_PredObsv_PC3)=c("V1")

dpGBC_PredObsv_PC4 <- as.data.frame(read.table("BSNIP_PC4_dpGBCPred_dpGBCObsv_r.txt", header=FALSE)[1:717,])
dpGBCObsv_PC4ScorePred <- as.data.frame(read.table("BSNIP_PC4_dpGBCObs_PC4ScorePred_r.txt", header=FALSE)[1:717,])
names(dpGBCObsv_PC4ScorePred)=c("V1")
names(dpGBC_PredObsv_PC4)=c("V1")

dpGBC_PredObsv_PC5 <- as.data.frame(read.table("BSNIP_PC5_dpGBCPred_dpGBCObsv_r.txt", header=FALSE)[1:717,])
dpGBCObsv_PC5ScorePred <- as.data.frame(read.table("BSNIP_PC5_dpGBCObs_PC5ScorePred_r.txt", header=FALSE)[1:717,])
names(dpGBCObsv_PC5ScorePred)=c("V1")
names(dpGBC_PredObsv_PC5)=c("V1")

# -- Plot correlations between dpGBCobs and predicted symptom PC scores across step-down parcel selection (718 to 1)
for (pcno in 1:5){
	q <- ggplot(data=eval(parse(text=paste0("dpGBCObsv_PC",pcno,"ScorePred"))),aes(x=2:718, y=V1)) + 
	geom_point() +
	scale_x_continuous(limits=c(0,740),breaks=c(1,718),labels=c(718,1)) +
	scale_y_continuous(limits=c(0,0.38),breaks=c(0,0.35),labels=c(0,0.35)) +
	ylab(expression(r(dpGBC^"obs",PC~Score^pred))) +
	xlab("# Parcels\nin Model") +
	theme1 + 
	theme(axis.text.x = element_text(size = 25, color="black",margin=margin(5,0,0,0)), 
			axis.title.x = element_text(size = 25, margin=margin(-25,0,0,0)),
			axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,5,0,0)), 
			axis.title.y = element_text(size = 18, hjust=0.3,margin=margin(0,-25,0,0)),
			axis.line= element_line(colour = "black", size = 2),
			axis.ticks = element_line(colour = "black", size=2),
			axis.ticks.length = unit(.4, "cm"))
	print(q)
	dev.off()
}

# -- Plot correlations between dpGBCobs and dpGBCpred across step-down parcel selection (718 to 1)
for (pcno in 1:5){
	q <- ggplot(data=eval(parse(text=paste0("s",pcno))),aes(x=2:718, y=V1)) + 
	geom_point(col="grey40") +
	scale_x_continuous(limits=c(0,740),breaks=c(1,718),labels=c(718,1)) +
	scale_y_continuous(limits=c(0,0.38),breaks=c(0,0.35),labels=c(0,0.35)) +
	ylab(expression(r(dpGBC^"obs",dpGBC^"pred"))) +
	xlab("# Parcels\nin Model") +
	theme1 + 
	theme(axis.text.x = element_text(size = 25, color="black",margin=margin(5,0,0,0)), 
			axis.title.x = element_text(size = 25, margin=margin(-25,0,0,0)),
			axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,5,0,0)), 
			axis.title.y = element_text(size = 18, hjust=0.2,margin=margin(0,-25,0,0)),
			axis.line= element_line(colour = "black", size = 2),
			axis.ticks = element_line(colour = "black", size=2),
			axis.ticks.length = unit(.4, "cm"))
	print(q)
	dev.off()
}

pc1beta <- read.table("/Users/jielisaji/Dropbox/N-BRIDGE/BSNIP_Analyses/fcMRI/GBC_Parcellated718/GBC_PC1_2k_ColeGlasserGSRParcels718_dat_ztstat.txt", header=FALSE)
pc2beta <- read.table("/Users/jielisaji/Dropbox/N-BRIDGE/BSNIP_Analyses/fcMRI/GBC_Parcellated718/GBC_PC2_2k_ColeGlasserGSRParcels718_dat_ztstat.txt", header=FALSE)
pc3beta <- read.table("/Users/jielisaji/Dropbox/N-BRIDGE/BSNIP_Analyses/fcMRI/GBC_Parcellated718/GBC_PC3_2k_ColeGlasserGSRParcels718_dat_ztstat.txt", header=FALSE)
pc4beta <- read.table("/Users/jielisaji/Dropbox/N-BRIDGE/BSNIP_Analyses/fcMRI/GBC_Parcellated718/GBC_PC4_2k_ColeGlasserGSRParcels718_dat_ztstat.txt", header=FALSE)
pc5beta <- read.table("/Users/jielisaji/Dropbox/N-BRIDGE/BSNIP_Analyses/fcMRI/GBC_Parcellated718/GBC_PC5_2k_ColeGlasserGSRParcels718_dat_ztstat.txt", header=FALSE)

# -- Get P_select for which metrics are maximal
for (pcno in 1:4){
	# -- Get P_select maximal (get local maximum)
	for (i in 1:717){
		Pselect = which(eval(parse(text=paste0("pc",pcno,"ob_dotdGBCpc",pcno,"b_r$V1")))==sort(abs(eval(parse(text=paste0("pc",pcno,"ob_dotdGBCpc",pcno,"b_r$V1")))), decreasing=TRUE)[i])
		if ((Pselect[1]> 630) & (Pselect[1] < 687)){
			Pselect=Pselect-1
			break
		}
	}

	# -- Assign the index for which metrics are maximal
	assign(paste0("pc",pcno,"_maxind"), 717-Pselect)
	assign("pc5_maxind",31) #local max
}

# -- For PC1
pc1beta_thresh <- pc1beta$V1
for (parcel in 1:718){
	if (abs(pc1beta$V1[parcel]) < sort(abs(pc1beta$V1), decreasing=TRUE)[pc1_maxind]){
		pc1beta_thresh[parcel] <- 0
	} else {
		pc1beta_thresh[parcel] <-  pc1beta$V1[parcel]
	}
}
write.table(pc1beta_thresh, "GBC_PC1_Top79.txt", col.names=FALSE, row.names=FALSE)

pc2beta_thresh <- pc2beta$V1
for (parcel in 1:718){
	if (abs(pc2beta$V1[parcel]) < sort(abs(pc2beta$V1), decreasing=TRUE)[pc2_maxind]){
		pc2beta_thresh[parcel] <- 0
	} else {
		pc2beta_thresh[parcel] <-  pc2beta$V1[parcel]
	}
}
write.table(pc2beta_thresh, "GBC_PC2_Top68.txt", col.names=FALSE, row.names=FALSE)

pc3beta_thresh <- pc3beta$V1
for (parcel in 1:718){
	if (abs(pc3beta$V1[parcel]) < sort(abs(pc3beta$V1), decreasing=TRUE)[pc3_maxind]){
		pc3beta_thresh[parcel] <- 0
	} else {
		pc3beta_thresh[parcel] <-  pc3beta$V1[parcel]
	}
}
write.table(pc3beta_thresh, "GBC_PC3_Top39.txt", col.names=FALSE, row.names=FALSE)

pc4beta_thresh <- pc4beta$V1
for (parcel in 1:718){
	if (abs(pc4beta$V1[parcel]) < sort(abs(pc4beta$V1), decreasing=TRUE)[pc4_maxind]){
		pc4beta_thresh[parcel] <- 0
	} else {
		pc4beta_thresh[parcel] <-  pc4beta$V1[parcel]
	}
}
write.table(pc4beta_thresh, "GBC_PC4_Top50.txt", col.names=FALSE, row.names=FALSE)

pc5beta_thresh <- pc5beta$V1
for (parcel in 1:718){
	if (abs(pc5beta$V1[parcel]) < sort(abs(pc5beta$V1), decreasing=TRUE)[pc5_maxind]){
		pc5beta_thresh[parcel] <- 0
	} else {
		pc5beta_thresh[parcel] <-  pc5beta$V1[parcel]
	}
}#colSums(matrix(pc5beta_thresh) != 0)
write.table(pc5beta_thresh, "GBC_PC5_Top31.txt", col.names=FALSE, row.names=FALSE)

# -- Predicted vs Observed dpGBC calculated for Pselect_max (for PC3)
bsnippred = read.table("BSNIP_N436_dpGBC_PreddGBCDotPC3beta_Top39Parcels.txt", sep="\t")
bsnipobsv = read.table("BSNIP_N436_dpGBC_ObsvdGBCDotPC3beta_Top39Parcels.txt", sep="\t")
Group=c("SADP","SCZP","BPP","SCZP","SADP","SCZP","SCZP","BPP","BPP","BPP","SCZP","BPP","SCZP","BPP","SCZP","SCZP","BPP","BPP","BPP","SADP","SCZP","SCZP","SCZP","BPP","BPP","SADP","BPP","SADP","BPP","SCZP","SADP","BPP","SADP","SCZP","BPP","BPP","SADP","SADP","BPP","SCZP","SCZP","SADP","BPP","SADP","SCZP","SADP","SADP","SADP","SADP","SCZP","SCZP","SADP","SADP","SCZP","BPP","SCZP","BPP","SCZP","SCZP","SADP","BPP","BPP","BPP","BPP","SADP","SCZP","BPP","BPP","SADP","BPP","BPP","BPP","BPP","SCZP","BPP","BPP","SADP","SCZP","SADP","SADP","SCZP","SADP","SCZP","SADP","BPP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","SADP","SCZP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","BPP","SADP","SADP","SCZP","SCZP","BPP","BPP","SCZP","SADP","SCZP","SCZP","SCZP","BPP","SCZP","SADP","SADP","BPP","BPP","SCZP","SCZP","BPP","SCZP","BPP","SADP","SCZP","SADP","BPP","SADP","SADP","BPP","BPP","SADP","SADP","BPP","BPP","SADP","BPP","BPP","SCZP","SCZP","SCZP","BPP","BPP","BPP","BPP","SADP","SADP","SADP","SCZP","SCZP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SCZP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SADP","SCZP","BPP","BPP","SCZP","SADP","SCZP","SADP","SCZP","BPP","SCZP","BPP","BPP","SADP","BPP","BPP","BPP","BPP","BPP","BPP","BPP","BPP","SADP","SCZP","SCZP","BPP","SADP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SADP","SCZP","SCZP","BPP","BPP","BPP","BPP","BPP","BPP","SADP","SCZP","SADP","SADP","SCZP","SCZP","BPP","SCZP","BPP","BPP","SADP","BPP","SADP","SCZP","SADP","SCZP","SADP","SCZP","BPP","SCZP","SCZP","BPP","SADP","SADP","SADP","BPP","SADP","SCZP","SADP","SADP","BPP","SADP","SCZP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SCZP","SADP","BPP","SCZP","BPP","BPP","SCZP","BPP","BPP","SCZP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SADP","SCZP","BPP","SCZP","SCZP","SADP","SADP","BPP","BPP","SADP","SADP","BPP","BPP","SCZP","SCZP","SADP","SADP","BPP","BPP","BPP","BPP","SCZP","SCZP","SCZP","SCZP","SCZP","BPP","SADP","SCZP","SCZP","SCZP","BPP","SADP","SCZP","SADP","SCZP","SCZP","BPP","BPP","SCZP","SCZP","SADP","BPP","SADP","BPP","SADP","SADP","SADP","SADP","SCZP","BPP","SADP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SADP","SADP","SADP","SCZP","BPP","SCZP","SCZP","SCZP","SADP","BPP","BPP","SCZP","SCZP","BPP","SADP","SADP","SCZP","SADP","SADP","SADP","SCZP","SADP","SCZP","SADP","SADP","SADP","SADP","SCZP","BPP","SCZP","BPP","SADP","SCZP","SCZP","BPP","SADP","SCZP","BPP","SCZP","BPP","SCZP","SADP","BPP","SCZP","SCZP","BPP","SCZP","BPP","BPP","BPP","SADP","SCZP","BPP","SCZP","BPP","SCZP","BPP","SCZP","SADP","SCZP","SADP","BPP","SADP","SCZP","SCZP","SCZP","SADP","BPP","SADP","BPP","SCZP","BPP","BPP","BPP","BPP","SADP","SADP","SADP","SCZP","SADP","SCZP","BPP","SADP","SCZP","BPP","SCZP","BPP","BPP","BPP","SCZP")

bsnippred=data.frame("Group"=Group,"Value"=c(bsnippred))
bsnipobsv=data.frame("Group"=Group,"Value"=c(bsnipobsv))
names(bsnippred) = c("Group","Value")
names(bsnipobsv) = c("Group","Value")

# -- Compute correlations between obs and pred dpGBC
bsnipsimdf <- data.frame("Dx"=c("SCZP", "SADP", "BPP", "PSD"),
"r"=c(cor(bsnipobsv$Value[which(bsnipobsv$Group=="SCZP")],bsnippred$Value[which(bsnippred$Group=="SCZP")]),
  cor(bsnipobsv$Value[which(bsnipobsv$Group=="SADP")],bsnippred$Value[which(bsnippred$Group=="SADP")]),
  cor(bsnipobsv$Value[which(bsnipobsv$Group=="BPP")],bsnippred$Value[which(bsnippred$Group=="BPP")]),
  cor(bsnipobsv$Value,bsnippred$Value)))

# -- Plot summary broken down by Dx group
ggplot(data=bsnipsimdf, aes(x=Dx,y=r, fill=Dx)) +
  geom_col() +
  geom_text(aes(label=c("n=167","n=119","n=150","n=436")), size=5, vjust=-1) +
  scale_fill_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026", "OCD"="cornflowerblue","PSD"="black")) +
  theme1 +
  scale_x_discrete(limits=c("PSD","BPP","SADP","SCZP"),labels=c("PSD","BPP", "SADP","SCZP")) +
  scale_y_continuous(expand = c(0.001, 0), limits=c(0,0.45), breaks=seq(0,0.4,0.2),labels=format(seq(0,0.4,0.2), nsmall=1)) +
  xlab("") +
  ylab("r(Pred,Obsv)\nDot Products") +
  theme(axis.text.x = element_text(size = 25, color="black", margin=margin(5,0,0,0), angle=60,hjust=1, vjust=1), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 25, color="black", angle=90, hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,5,0,0)),
        axis.line= element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size=2),
        axis.ticks.length = unit(.4, "cm"))
dev.off()

# -- Plot the single subject example
preddgbc_t39 <- as.data.frame(t(read.table("BSNIP_N436_PredGBC_Top39Parcels.txt")))
obsvdgbc_t39 <- as.data.frame(t(read.table("BSNIP_N436_ObsvGBC_Top39Parcels.txt")))
bsnipsub1 <- data.frame("Predicted"=preddgbc_t39$V380, "Observed"=obsvdgbc_t39$V380)
ggplot(data=bsnipsub1, aes(x=Predicted,y=Observed, col=Predicted)) +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  geom_point(size=4, alpha=0.7, col="grey30") +
  theme1 +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-0.0026,0.0026), breaks=seq(-0.002,0.002,length=3),labels=seq(-0.002,0.002,length=3)) +
  scale_y_continuous(expand = c(0.01, 0), limits=c(-0.025,0.025), breaks=seq(-0.02,0.02,length=3),labels=seq(-0.02,0.02,length=3)) +
  xlab(expression("Predicted "*Delta*"GBC")) +
  ylab(expression("Observed "*Delta*"GBC")) +
  theme(axis.text.x = element_text(size = 25, color="black",margin=margin(5,0,0,0)), 
        axis.title.x = element_text(size = 25, margin=margin(5,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,5,0,0)),
        axis.line= element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size=2),
        axis.ticks.length = unit(.4, "cm"))
dev.off()