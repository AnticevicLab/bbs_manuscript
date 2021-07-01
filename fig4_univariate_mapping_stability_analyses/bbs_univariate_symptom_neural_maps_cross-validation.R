### -- This script is a continuation of the "../fig2_symptom_pca_stability_analyses/bbs_behavior_pca_cross-validation.R" script

#setwd("/Users/jielisaji/Dropbox/bbs_manuscript/fig4_univariate_mapping_stability_analyses")

###############################
##### LOO Site Validation #####
###############################

# -- Organizing the data for plotting purposes 
## 1-CT ; 2-GP ; 3-GT ; 4-JS ; 5-MB ; 6-MK
sdat<-NULL
for (site in c("CT","GP","GT","JS","MB","MK")){
  for (pc in 1:5){
    sdat <- data.frame(rbind(sdat, cbind(site, pc, abs(cor(eval(parse(text=paste0("pred",site,"$PC",pc))),eval(parse(text=paste0("obs",site,"$PC",pc))))))))
}}
names(sdat) <- c("Site","PC","r")
sdat$r <- unfactor(sdat$r)
sdat$Site <- as.character(sdat$Site)
sdat$Site[sdat$Site == "CT"] <- 1
sdat$Site[sdat$Site == "GP"] <- 2
sdat$Site[sdat$Site == "GT"] <- 3
sdat$Site[sdat$Site == "JS"] <- 4
sdat$Site[sdat$Site == "MB"] <- 5
sdat$Site[sdat$Site == "MK"] <- 6
sdat$Site <- as.factor(sdat$Site)

# -- Plot the exemplar scatterplot of betas for one run
ggplot(data=PC3all, aes(x=V1,y=PC3S3$V1, col=PC3S3$V1)) +
  geom_point(size=3, alpha=0.8) +
  geom_smooth(method = lm, se=FALSE, col="forest green",lwd=4) +
  scale_color_gradientn(colours=c("cyan","blue","black","red","yellow")) +
  theme1 +
  theme(axis.text.x = element_text(size = 40, color="black"), 
        axis.title.x = element_text(size = 40,margin=margin(15,0,0,0)),
        axis.text.y = element_text(size = 40, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
        axis.title.y = element_text(size = 40,margin=margin(0,15,0,0))) +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-5,5.2), breaks=seq(-5,5,5),labels=seq(-5,5,5)) +
  scale_y_continuous(expand = c(0.01, 0), limits=c(-5,5.2), breaks=seq(-5,5,5),labels=seq(-5,5,5)) +
  xlab("Full model betas") +
  ylab("W/o Site 3 betas") +
  coord_fixed() +
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(15,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,10,0,0)))
dev.off()
cor.test(PC3S3$V1,PC3all$V1)

### Site effects
# -- Read in the data from the beta_PC_GBC maps from LOSO
for (map in c("Positive","General","Negative","Cognitive","PC1","PC2","PC3","PC4","PC5")){
  assign(paste0(map,"S1") ,read.table(paste0("leave_one_site_out_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_",map,"_SITE1trainN357_beta.txt")))
  assign(paste0(map,"S2") ,read.table(paste0("leave_one_site_out_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_",map,"_SITE2trainN339_beta.txt")))
  assign(paste0(map,"S3") ,read.table(paste0("leave_one_site_out_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_",map,"_SITE3trainN329_beta.txt")))
  assign(paste0(map,"S4") ,read.table(paste0("leave_one_site_out_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_",map,"_SITE4trainN339_beta.txt")))
  assign(paste0(map,"S5") ,read.table(paste0("leave_one_site_out_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_",map,"_SITE5trainN426_beta.txt")))
  assign(paste0(map,"S6") ,read.table(paste0("leave_one_site_out_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_",map,"_SITE6trainN390_beta.txt")))
}

PC1all       <-read.table("../fig3_univariate_mapping_analyses/univariate_maps/GBC_PC1_ztstat.pscalar.txt")
PC2all       <-read.table("../fig3_univariate_mapping_analyses/univariate_maps/GBC_PC2_ztstat.pscalar.txt")
PC3all       <-read.table("../fig3_univariate_mapping_analyses/univariate_maps/GBC_PC3_ztstat.pscalar.txt")
PC4all       <-read.table("../fig3_univariate_mapping_analyses/univariate_maps/GBC_PC4_ztstat.pscalar.txt")
PC5all       <-read.table("../fig3_univariate_mapping_analyses/univariate_maps/GBC_PC5_ztstat.pscalar.txt")
Positiveall  <-read.table("../fig3_univariate_mapping_analyses/univariate_maps/GBC_Pos_ztstat.pscalar.txt")
Negativeall  <-read.table("../fig3_univariate_mapping_analyses/univariate_maps/GBC_Neg_ztstat.pscalar.txt")
Generalall   <-read.table("../fig3_univariate_mapping_analyses/univariate_maps/GBC_Gen_ztstat.pscalar.txt")
Cognitiveall <-read.table("../fig3_univariate_mapping_analyses/univariate_maps/GBC_Cog_ztstat.pscalar.txt")

# -- Organise the dataframe for plotting
corSDF <- NULL
for (pc in c("Cognitive","Negative","Positive","General","PC1","PC2","PC3","PC4","PC5")){
  for (kval in 1:6){
    corSDF <- rbind(corSDF, cbind(pc, kval, as.numeric(cor.test(eval(parse(text=paste0(pc,"S",kval,"$V1"))),eval(parse(text=paste0(pc,"all$V1"))))$estimate[[1]])))
  }}
corSDF <- as.data.frame(corSDF)
names(corSDF) <- c("Map", "Site", "r")
corSDF$r<-dmm::unfactor(corSDF$r)
corSDF$Map <- factor(corSDF$Map, levels = c("Cognitive","Negative","Positive","General","PC1","PC2","PC3","PC4","PC5"))

# -- Plot the summary of LOSO neural mapping for all maps
cols <- rev(viridisLite::plasma(10))
cols <- cols[-1]
cols2 <- viridisLite::viridis(10)
ggplot(corSDF, aes(x=Map, y=r, col=Map,fill=factor(Map),shape=Site)) + 
  geom_point(size=10,alpha=0.8) +
  theme1 + 
  scale_color_manual(values=c("Cognitive"="#1FE20C","Negative"="#00B1FF","Positive"="#8A40F5","General"="pink","PC1"=cols[1],"PC2"=cols[2],"PC3"=cols[3],"PC4"=cols[4],"PC5"=cols[5],"PRB"=cols2[1],"BPP"=cols2[2],"SADP"=cols2[3],"SZP"=cols2[4])) +
  scale_fill_manual(values=c("Cognitive"="#1FE20C","Negative"="#00B1FF","Positive"="#8A40F5","General"="pink","PC1"=cols[1],"PC2"=cols[2],"PC3"=cols[3],"PC4"=cols[4],"PC5"=cols[5],"PRB"=cols2[1],"BPP"=cols2[2],"SADP"=cols2[3],"SZP"=cols2[4])) +
  scale_shape_manual(values = c("1"=21,"2"=22,"3"=23, "4"=24, "5"=25, "6"=8)) +
  scale_y_continuous(limits=c(0, 1.02),breaks=seq(0,1,0.5),labels=format(seq(0,1,0.5), nsmall=1)) +
  scale_x_discrete(labels=c("PC1","PC2","PC3","PC4","PC5","Cog","Pos","Neg","Gen"),limits=c("PC1","PC2","PC3","PC4","PC5","Cognitive","Positive","Negative","General")) +
  ylab("r-value") + 
  xlab("Behavioral Map") +
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(15,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,10,0,0)))
dev.off()


###############################
###### 5-fold Validation ######
###############################
#### Main Text Summary Figures
PC1k1 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC1_K1trainN349_beta.txt")
PC1k2 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC1_K2trainN349_beta.txt")
PC1k3 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC1_K3trainN349_beta.txt")
PC1k4 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC1_K4trainN348_beta.txt")
PC1k5 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC1_K5trainN349_beta.txt")

PC2k1 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC2_K1trainN349_beta.txt")
PC2k2 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC2_K2trainN349_beta.txt")
PC2k3 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC2_K3trainN349_beta.txt")
PC2k4 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC2_K4trainN348_beta.txt")
PC2k5 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC2_K5trainN349_beta.txt")

PC3k1 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC3_K1trainN349_beta.txt")
PC3k2 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC3_K2trainN349_beta.txt")
PC3k3 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC3_K3trainN349_beta.txt")
PC3k4 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC3_K4trainN348_beta.txt")
PC3k5 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC3_K5trainN349_beta.txt")

PC4k1 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC4_K1trainN349_beta.txt")
PC4k2 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC4_K2trainN349_beta.txt")
PC4k3 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC4_K3trainN349_beta.txt")
PC4k4 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC4_K4trainN348_beta.txt")
PC4k5 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC4_K5trainN349_beta.txt")

PC5k1 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC5_K1trainN349_beta.txt")
PC5k2 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC5_K2trainN349_beta.txt")
PC5k3 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC5_K3trainN349_beta.txt")
PC5k4 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC5_K4trainN348_beta.txt")
PC5k5 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_PC5_K5trainN349_beta.txt")

Posk1 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Positive_K1trainN349_beta.txt")
Posk2 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Positive_K2trainN349_beta.txt")
Posk3 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Positive_K3trainN349_beta.txt")
Posk4 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Positive_K4trainN348_beta.txt")
Posk5 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Positive_K5trainN349_beta.txt")

Negk1 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Negative_K1trainN349_beta.txt")
Negk2 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Negative_K2trainN349_beta.txt")
Negk3 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Negative_K3trainN349_beta.txt")
Negk4 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Negative_K4trainN348_beta.txt")
Negk5 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Negative_K5trainN349_beta.txt")

Genk1 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_General_K1trainN349_beta.txt")
Genk2 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_General_K2trainN349_beta.txt")
Genk3 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_General_K3trainN349_beta.txt")
Genk4 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_General_K4trainN348_beta.txt")
Genk5 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_General_K5trainN349_beta.txt")

Cogk1 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Cognitive_K1trainN349_beta.txt")
Cogk2 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Cognitive_K2trainN349_beta.txt")
Cogk3 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Cognitive_K3trainN349_beta.txt")
Cogk4 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Cognitive_K4trainN348_beta.txt")
Cogk5 <- read.table("k-fold_cross-validation_maps/BSNIP_PSD_N436_GSR.udvarsme.CAB-NP-718_gbc_mFz_Cognitive_K5trainN349_beta.txt")

# -- Scatterplot of correlation between parcel coefficients
ggplot(data=PC3all, aes(x=V1,y=PC3k1$V1, col=PC3all$V1)) +
  geom_point(size=3, alpha=0.8) +
  geom_smooth(method = lm, se=FALSE, col="forest green",lwd=4) +
  scale_color_gradientn(colours=c("cyan","blue","black","red","yellow")) +
  theme1 +
  theme(axis.text.x = element_text(size = 40, color="black"), 
        axis.title.x = element_text(size = 40,margin=margin(15,0,0,0)),
        axis.text.y = element_text(size = 40, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
        axis.title.y = element_text(size = 40,margin=margin(0,15,0,0))) +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-5,5.2), breaks=seq(-5,5,5),labels=seq(-5,5,5)) +
  scale_y_continuous(expand = c(0.01, 0), limits=c(-5,5.2), breaks=seq(-5,5,5),labels=seq(-5,5,5)) +
  xlab("Full model betas") +
  ylab("Fold 1 betas") +
  coord_fixed() +
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(15,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,10,0,0)))
dev.off()

# -- Organize the dataframe for plotting purposes
corDF <- NULL
for (pc in c("Cog","Neg","Pos","Gen","PC1","PC2","PC3","PC4","PC5")){
  for (kval in 1:5){
  corDF <- rbind(corDF, cbind(pc, kval, as.numeric(cor.test(eval(parse(text=paste0(pc,"k",kval,"$V1"))),eval(parse(text=paste0(pc,"all$V1"))))$estimate[[1]])))
}}
corDF <- as.data.frame(corDF)
names(corDF) <- c("Map", "k", "r")
corDF$r<-dmm::unfactor(corDF$r)

# -- Plot the summary boxplot figure of 5-fold CV for all PCs and A Priori maps
cols <- rev(viridisLite::plasma(10))
cols <- cols[-1]
ggplot(corDF, aes(x=Map, y=r, fill=factor(Map))) +
  geom_boxplot(size=1) +
  theme1 + 
  scale_fill_manual(values=c("Cog"="#1FE20C","Neg"="#00B1FF","Pos"="#8A40F5","Gen"="pink","PC1"=cols[1],"PC2"=cols[2],"PC3"=cols[3],"PC4"=cols[4],"PC5"=cols[5],"PRB"=cols2[1],"BPP"=cols2[2],"SADP"=cols2[3],"SZP"=cols2[4])) +
  scale_color_manual(values=c("Cog"="#1FE20C","Neg"="#00B1FF","Pos"="#8A40F5","Gen"="pink","PC1"=cols[1],"PC2"=cols[2],"PC3"=cols[3],"PC4"=cols[4],"PC5"=cols[5],"PRB"=cols2[1],"BPP"=cols2[2],"SADP"=cols2[3],"SZP"=cols2[4])) +
  scale_y_continuous(limits=c(0, 1.02),breaks=seq(0,1,0.5),labels=format(seq(0,1,0.5), nsmall=1)) +
  scale_x_discrete(limits=c("PC1","PC2","PC3","PC4","PC5","Cog","Pos","Neg","Gen")) +
  ylab("r-value") +
  xlab("Behavioral Map") +
  theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(15,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,10,0,0)))
dev.off()







