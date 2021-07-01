# -- This script is a continuation of ../fig6_feature_selection_analyses/bsnip_featureselection.R
setwd("/Users/jielisaji/Dropbox/bbs_manuscript/fig7_pharmacological_neuroimaging_analyses")

# -- Get the demeaned delta GBC for BSNIP PSD subjects
bsnip_gbcdm <- t(read.table("../fig6_feature_selection_analyses/BSNIP_PSD_N436_deltaGBC.txt"))

# -- Ketamine for top PC3 parcels
pc3beta_thresh_clean <- pc3beta_thresh
pc3beta_thresh_clean[pc3beta_thresh_clean==0] <- NA
ketgbc_top39 <- read.table("GBC_Ketamine_Response_Top39Parcels.txt")
ketgbc_top39_clean <- ketgbc_top39
ketgbc_top39_clean[ketgbc_top39_clean==0] <- NA
cor(pc3beta_thresh_clean, ketgbc_top39_clean, use="pairwise.complete.obs")

# -- LSD for top PC3 parcels
lsdgbc_top39 <- read.table("GBC_LSD_Response_Top39Parcels.txt")
lsdgbc_top39_clean <- lsdgbc_top39
lsdgbc_top39_clean[lsdgbc_top39_clean==0] <- NA
cor(pc3beta_thresh_clean, lsdgbc_top39_clean, use="pairwise.complete.obs")

# -- Project independent data
pc3beta_threshmask = pc3beta_thresh != 0
repl_gbcdmrepl <- read.table("Replication_OCDSCZ_N39_GSR.udvarsme.CAB-NP-718_gbc_mFz_demean.txt")
repl_gbcdmrepl_top39parcels <- (matrix(rep(pc3beta_threshmask, 69),718,69) * repl_gbcdmrepl)

# -- Compute r across top PC3 parcels for all subjects in Replication sample
repl_gbcdmrepl_top39parcels_clean <- repl_gbcdmrepl_top39parcels[-which(repl_gbcdmrepl_top39parcels==0),]
pc3beta_thresh_clean <- pc3beta_thresh[-which(pc3beta_thresh==0)]
repl_gbcpc3featcor <- cor(repl_gbcdmrepl_top39parcels_clean, pc3beta_thresh_clean, method="spearman")

# -- Compute r across top PC5 parcels for all subjects in BSNIP sample
bsnip_gbcdmbsnip_top31parcels_clean <- bsnip_gbcdm[-which(pc5beta_thresh==0),]
pc5beta_thresh_clean <- pc5beta_thresh[-which(pc5beta_thresh==0)]
bsnip_gbcpc5featcor <- cor(bsnip_gbcdmbsnip_top31parcels_clean, pc5beta_thresh_clean, method="spearman")

# -- Similarity index bar for Subject X_PC3
X_PC3 <- data.frame("rank"="A", variable=factor(c("1","2","3")),"value"=c(-1-repl_gbcpc3featcor[42,1],repl_gbcpc3featcor[42,1], 1))
ggplot(X_PC3, aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept=0, lwd=2) +
  scale_x_discrete(expand=c(0.0001,0)) +
  scale_y_continuous(expand=c(0.001,0),limits=c(-1,1), breaks=c(-1,0,1),labels=c(" -1",0,"1 "), position="right") +
  scale_fill_manual(values=c("1"="grey","2"="#5A95D5","3"="grey")) +
  theme1 + 
  ylab("Spearman's rho") +
  theme(axis.line.x=element_blank(),
  		axis.text.x=element_blank(),
  		axis.title.x=element_blank(),
  		axis.ticks.x=element_blank(),
  		axis.title.y.right=element_text(size=35,angle=90,margin=margin(0,5,0,0)),
  		axis.text.y.right=element_text(size=35, hjust=0.5,margin=margin(0,0,0,5)))
dev.off()

# -- Similarity index bar for Subject Y_PC3
Y_PC3 <- data.frame("rank"="A", variable=factor(c("1","2","3")),"value"=c(-1, 1-repl_gbcpc3featcor[55,1],repl_gbcpc3featcor[55,1]))
ggplot(Y_PC3, aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept=0, lwd=2) +
  scale_x_discrete(expand=c(0.0001,0)) +
  scale_y_continuous(expand=c(0.001,0),limits=c(-1,1), breaks=c(-1,0,1),labels=c(" -1",0,"1 "), position="right") +
  scale_fill_manual(values=c("1"="grey","2"="grey","3"="#DCB13B")) +
  theme1 + 
  ylab("Spearman's rho") +
  theme(axis.line.x=element_blank(),
  		axis.text.x=element_blank(),
  		axis.title.x=element_blank(),
  		axis.ticks.x=element_blank(),
  		axis.title.y.right=element_text(size=35,angle=90,margin=margin(0,5,0,0)),
  		axis.text.y.right=element_text(size=35, hjust=0.5,margin=margin(0,0,0,5)))
dev.off()

# -- Similarity index bar for Subject Q_PC5
Q_PC5 <- data.frame("rank"="A", variable=factor(c("1","2","3")),"value"=c(-1-bsnip_gbcpc5featcor[37,1],bsnip_gbcpc5featcor[37,1], 1))
ggplot(Q_PC5, aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept=0, lwd=2) +
  scale_x_discrete(expand=c(0.0001,0)) +
  scale_y_continuous(expand=c(0.001,0),limits=c(-1,1), breaks=c(-1,0,1),labels=c(" -1",0,"1 "), position="right") +
  scale_fill_manual(values=c("1"="grey","2"="#5A95D5","3"="grey")) +
  theme1 + 
  ylab("Spearman's rho") +
  theme(axis.line.x=element_blank(),
  		axis.text.x=element_blank(),
  		axis.title.x=element_blank(),
  		axis.ticks.x=element_blank(),
  		axis.title.y.right=element_text(size=35,angle=90,margin=margin(0,5,0,0)),
  		axis.text.y.right=element_text(size=35, hjust=0.5,margin=margin(0,0,0,5)))
dev.off()

# -- Similarity index bar for Subject Z_PC5
Z_PC5 <- data.frame("rank"="A", variable=factor(c("1","2","3")),"value"=c(-1,1-bsnip_gbcpc5featcor[119,1],bsnip_gbcpc5featcor[119,1]))
ggplot(Z_PC5, aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept=0, lwd=2) +
  scale_x_discrete(expand=c(0.0001,0)) +
  scale_y_continuous(expand=c(0.001,0),limits=c(-1,1), breaks=c(-1,0,1),labels=c(" -1",0,"1 "), position="right") +
  scale_fill_manual(values=c("1"="grey","3"="#DCB13B","2"="grey")) +
  theme1 + 
  ylab("Spearman's rho") +
  theme(axis.line.x=element_blank(),
  		axis.text.x=element_blank(),
  		axis.title.x=element_blank(),
  		axis.ticks.x=element_blank(),
  		axis.title.y.right=element_text(size=35,angle=90,margin=margin(0,5,0,0)),
  		axis.text.y.right=element_text(size=35, hjust=0.5,margin=margin(0,0,0,5)))
dev.off()