rm(list = ls())
setwd("../fig5_cca_analyses")
library(ggplot2)
library(GGally)
library(CCA)
library(mixOmics)
library(CCP)
library(matrixStats)
library(ggridges)
library(fmsb)

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
mtsortdiag <- function(mt){
        mtf <- as.matrix(mt)
        mtb <- as.matrix(mt)

        ncs <- ncol(mt)-1

        for (nc in 1:ncs){
                if ((sum(abs(diag(mtf[nc:(nc+1),nc:(nc+1)])) < abs(diag(mtf[nc:(nc+1),(nc+1):nc]))) == 0) && (sum(abs(diag(mtf[nc:(nc+1),nc:(nc+1)])) < abs(diag(mtf[(nc+1):nc,nc:(nc+1)]))) == 0)){
                        #print("diagonal is maximal")
                } else {
                        if ( sum(abs(diag(mtf[nc:(nc+1),nc:(nc+1)]))) < sum(abs(diag(mtf[nc:(nc+1),(nc+1):nc]))) ){
                                #print("diagonal is not maximal")
                                mtf[,nc:(nc+1)] <- mtf[,(nc+1):nc]
                                names(mtf[,nc:(nc+1)]) <- names(mtf[,(nc+1):nc])
                        }
                }
        }

        for (nc in ncs:1){
                if ((sum(abs(diag(mtb[nc:(nc+1),nc:(nc+1)])) < abs(diag(mtb[nc:(nc+1),(nc+1):nc]))) == 0) && (sum(abs(diag(mtb[nc:(nc+1),nc:(nc+1)])) < abs(diag(mtb[(nc+1):nc,nc:(nc+1)]))) == 0)){
                        #print("diagonal is maximal")
                } else {
                        if ( sum(abs(diag(mtb[nc:(nc+1),nc:(nc+1)]))) < sum(abs(diag(mtb[nc:(nc+1),(nc+1):nc]))) ){
                                #print("diagonal is not maximal")
                                mtb[,nc:(nc+1)] <- mtb[,(nc+1):nc]
                                names(mtb[,nc:(nc+1)]) <- names(mtb[,(nc+1):nc])
                        }
                }
        }

        if ((sum(abs(diag(mtf)))>sum(abs(diag(mtb))))){
                return(mtf)
        } else {
                return(mtb)
        }
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

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
## Split-half CCA Replication #
###############################
# -- Read behavioral data
beh_raw <- read.table("../data_symptoms/BSNIP_PSD_N436_Behavior.dat", header=T)
beh <- beh_raw[,-1]
Group=beh_raw$Group

behpc_raw <- read.table("../fig1_symptom_pca_analyses/BSNIP_PSD436_AllPCScores.txt")
behpc <- behpc_raw[,2:6]

# -- Read neural data
gbcPSD_180pS <- t(read.table("BSNIP_PSD_N436_GBC_180P_SymmetrizedCortex.txt"))

nRep=1000
nB=5
H1shuffPC <- array(NA,dim=c(217,nB,nRep))
H2shuffPC <- array(NA,dim=c(219,nB,nRep))
H1indices <- matrix(NA, 217, nRep)
H2indices <- matrix(NA, 219, nRep)

shXcoef <- NULL
shYcoef <- NULL

PCloadings <- as.data.frame(read.table("../fig1_symptom_pca_analyses/BSNIP_N436_BACSPANSS_PCArotations.txt", sep=" ", header=TRUE))

for (shuff in 1:nRep){
	sczpind <- which(Group=="SCZP")
	sczpind.mixed <- sample(sczpind, size = 0.5*length(sczpind), replace = FALSE)
	
	sadpind <- which(Group=="SADP")
	sadpind.mixed <- sample(sadpind, size = 0.5*length(sadpind), replace = FALSE)
	
	bppind <- which(Group=="BPP")
	bppind.mixed <- sample(bppind, size = 0.5*length(bppind), replace = FALSE)
	
	h1.mixed <- sort(c(sczpind.mixed,sadpind.mixed,bppind.mixed))
	H1indices[,shuff] <- h1.mixed
	H2indices[,shuff] <- c(1:length(Group))[-h1.mixed]
	
	H1shuffPC[,,shuff] <- as.matrix(behpc[h1.mixed,])
	H2shuffPC[,,shuff] <- as.matrix(behpc[-h1.mixed,])
}

H1behw <- array(NA,dim=c(36,5,nRep))
H2behw <- array(NA,dim=c(36,5,nRep))
behwcor <- NULL

# -- Run for H1 and H2
for (shuff in 1:nRep){
        # -- Run CCA on H1
	H1ndat <- data.frame(gbcPSD_180pS[H1indices[,shuff],])
	H1bdat <- data.frame(H1shuffPC[,,shuff])
	names(H1bdat) <- c("PC1","PC2","PC3","PC4","PC5")
	ccaH1 <- cc(H1bdat, H1ndat)
	
        # -- Run CCA on H2
	H2ndat <- data.frame(gbcPSD_180pS[H2indices[,shuff],])
	H2bdat <- data.frame(H2shuffPC[,,shuff])
	names(H2bdat) <- c("PC1","PC2","PC3","PC4","PC5")
        ccaH2 <- cc(H2bdat, H2ndat)

	# -- Correlate H1 and H2 coefficients
        xcor <- data.frame(cor(ccaH1$xcoef, ccaH2$xcoef))
	ycor <- data.frame(cor(ccaH1$ycoef, ccaH2$ycoef))
	names(xcor) <- c("CV1","CV2","CV3","CV4","CV5")
	names(ycor) <- c("CV1","CV2","CV3","CV4","CV5")

	# -- Save the coefficients
        shXcoef <- c(shXcoef, diag(as.matrix(mtsortdiag(xcor))))
	shYcoef <- c(shYcoef, diag(as.matrix(mtsortdiag(ycor))))
	H1coef <- as.data.frame(ccaH1$xcoef)
	H2coef <- as.data.frame(ccaH2$xcoef)
	names(H1coef) <- c("CV1","CV2","CV3","CV4","CV5")
	names(H2coef) <- c("CV1","CV2","CV3","CV4","CV5")
	
	# -- Project the loadings for the original behavioral symptom measures on each CV
        H1behw[,,shuff] <- as.matrix(PCloadings[,1:5]) %*% as.matrix(H1coef[colnames(mtsortdiag(xcor))])
	H2behw[,,shuff] <- as.matrix(PCloadings[,1:5]) %*% as.matrix(H2coef[colnames(mtsortdiag(ycor))])
	behwcor <- c(behwcor, diag(cor(H1behw[,,shuff],H2behw[,,shuff])))
}

# -- Organize dataframes for plotting below
XcoefshRepDF <- data.frame("CV"=1:5, "R"=abs(shXcoef))
DFsummshXcoefRep <- summarySE(XcoefshRepDF, measurevar="R", groupvars=c("CV"))

YcoefshRepDF <- data.frame("CV"=1:5, "R"=abs(shYcoef))
DFsummshYcoefRep <- summarySE(YcoefshRepDF, measurevar="R", groupvars=c("CV"))

behwcorshRepDF <- data.frame("CV"=1:5, "R"=abs(behwcor))
DFsummshbehwRep <- summarySE(behwcorshRepDF, measurevar="R", groupvars=c("CV"))

# -- Plot the summary barchart of symptom PC loadings across 1k split-half runs for 5 CVs
pdf("CCA_SplitHalf_PCWeights.pdf",height=6,width=8)
ggplot(DFsummshXcoefRep, aes(x=CV, y=R)) +
  geom_hline(yintercept=1, col="red", lwd=1.5) +
  geom_bar(position="dodge", width=0.9, stat = "summary", fun.y = "mean",, fill=c("#485C77")) +
  geom_errorbar(aes(x=CV, ymin=R-se, ymax=R+se), width=.5, position=position_dodge(0.9), lwd=1.5) +
  scale_alpha_manual(values = c("1"="1","2"="0.9","3"="0.8", "4"="0.7", "5"="0.6")) +
  scale_y_continuous(expand=c(0.0001,0), limits=c(0.0, 1.1),breaks=seq(0,1,0.5),labels=format(seq(0,1,0.5), nsmall=1)) +
  scale_x_continuous(expand=c(0.000001,0),limits=c(0.45,5.5),breaks=seq(1,5,by=1),labels=1:5) +
  ylab("r-value") + 
  xlab("CV") +
  theme1 +
  theme(axis.text.x = element_text(size = 25, color="black"),
        axis.title.x = element_text(size = 25, margin=margin(15,0,0,0)),
        axis.title.y = element_text(size = 25,margin=margin(0,10,0,0)),
        axis.line = element_line(colour = "black", size = 1))
dev.off()

# -- Plot the summary barchart of neural loadings across 1k split-half runs for 5 CVs
pdf("CCA_SplitHalf_NeuralWeights.pdf",height=6,width=8)
ggplot(DFsummshYcoefRep, aes(x=CV, y=R)) +
  geom_hline(yintercept=1, col="red", lwd=1.5) +
  geom_bar(position="dodge", width=0.9, stat = "summary", fun.y = "mean",, fill=c("#485C77")) +
  geom_errorbar(aes(x=CV, ymin=R-se, ymax=R+se), width=.5, position=position_dodge(0.9), lwd=1.5) +
  scale_alpha_manual(values = c("1"="1","2"="0.9","3"="0.8", "4"="0.7", "5"="0.6")) +
  scale_y_continuous(expand=c(0.0001,0), limits=c(0.0, 1.1),breaks=seq(0,1,0.5),labels=c("","","")) +
  scale_x_continuous(expand=c(0.000001,0),limits=c(0.45,5.5),breaks=seq(1,5,by=1),labels=1:5) +
  ylab("r-value") + 
  xlab("CV") +
  theme1 +
  theme(axis.text.x = element_text(size = 25, color="black"),
        axis.title.x = element_text(size = 25, margin=margin(15,0,0,0)),
        axis.title.y = element_text(size = 25,margin=margin(0,10,0,0)),
        axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.line = element_line(colour = "black", size = 1))
dev.off()

pdf("CCA_SplitHalf_SymptomMeasureWeights.pdf",height=6,width=8)
ggplot(DFsummshbehwRep, aes(x=CV, y=R)) +
  geom_hline(yintercept=1, col="red", lwd=1.5) +
  geom_bar(position="dodge", width=0.9, stat = "summary", fun.y = "mean",, fill=c("#485C77")) +
  geom_errorbar(aes(x=CV, ymin=R-se, ymax=R+se), width=.5, position=position_dodge(0.9), lwd=1.5) +
  scale_alpha_manual(values = c("1"="1","2"="0.9","3"="0.8", "4"="0.7", "5"="0.6")) +
  scale_y_continuous(expand=c(0.0001,0), limits=c(0.0, 1.1),breaks=seq(0,1,0.5),labels=c("","","")) +
  scale_x_continuous(expand=c(0.000001,0),limits=c(0.45,5.5),breaks=seq(1,5,by=1),labels=1:5) +
  ylab("r-value") + 
  xlab("CV") +
  theme1 +
  theme(axis.text.x = element_text(size = 25, color="black"),
        axis.title.x = element_text(size = 25, margin=margin(15,0,0,0)),
        axis.title.y = element_text(size = 25,margin=margin(0,10,0,0)),
        axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.line = element_line(colour = "black", size = 1))
dev.off()

# -- Plot the scatterplot between N loadings for H1 and H2 in one exemplar run
pdf("CCA_SplitHalf_CV3_CorrelationScatterplot.pdf",height=6,width=6)
cv3neuralSHDF <- data.frame("H1"=ccaH1$ycoef[,3], "H2"=ccaH2$ycoef[,3])
ggplot(cv3neuralSHDF, aes(x=H1, y=H2, col=H2)) +
  geom_point(size=3, alpha=0.7)+
  scale_color_gradientn(colours=c("navy")) +
  scale_y_continuous(limits=c(-88,63),breaks=seq(-80,60,length=2),labels=seq(-80,60,length=2)) +
  scale_x_continuous(limits=c(-88,63),breaks=seq(-80,60,length=2),labels=seq(-80,60,length=2)) +
  xlab("H1 N Weights") +
  ylab("H2 N Weights") +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme1 + 
	theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(5,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,5,0,0)),
        axis.line = element_line(colour = "black", size = 1))
dev.off()

###############################
###### LOO CCA Validation #####
###############################
nCV=5
nParcels=180
nsubs=436

statX = array(NA,dim=c(nsubs,nCV))
statY = array(NA,dim=c(nsubs,nParcels))
statrho = array(NA,dim=c(nsubs,nCV))

scoresX <- array(NA,dim=c(nsubs-1,nCV,nsubs))
scoresY <- array(NA,dim=c(nsubs-1,nCV,nsubs))

coefX <- array(NA,dim=c(nCV,nCV,nsubs))
coefY <- array(NA,dim=c(nParcels,nCV,nsubs))

predX <- array(NA,dim=c(nsubs,nCV))
predY <- array(NA,dim=c(nsubs,nCV))

X <- behpc
Y <- gbcPSD_180pS

# -- Compute leave-one-subject-out CCA
for (i in 1:nsubs) {
	cci <- cc(X[-i, ], Y[-i, ]) # run CCA with all but subject i
	statX[i,] <- colMeans(cci$scores$corr.X.xscores^2)*(cci$cor^2) # -- Save the neural variance explained
	statY[i,] <- colMeans(cci$scores$corr.Y.yscores^2)*(cci$cor^2) # -- Save the behavioral variance explained
	statrho[i,] <- cci$cor # -- Save canonical correlation values
	# save the N and B scores
	scoresX[,,i] <- cci$scores$xscores
	scoresY[,,i] <- cci$scores$yscores
	# save the N and B coefs
	coefX[,,i] <- cci$xcoef
	coefY[,,i] <- cci$ycoef
	
	# predict left out subject
	predX[i,] <- as.matrix(X[i,]) %*% as.matrix(cci$xcoef)
	predY[i,] <- t(as.matrix(Y[i,])) %*% as.matrix(cci$ycoef)
}

# -- Plot exemplar CV3 leave-one-out CCA
cv3predDF <- data.frame("PredX"=predX[,3], "PredY"=predY[,3])

pdf("CCA_LeaveOneOut_PredScores_CV3.pdf", height=6, width=6)
ggplot(cv3predDF, aes(x=PredX, y=PredY, col=Group)) +
  geom_point(size=3, alpha=0.7)+
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_x_continuous(limits=c(-4.2,4.5),breaks=seq(-4,4,length=2),labels=seq(-4,4,length=2)) +
  scale_y_continuous(limits=c(-4.2,5.5),breaks=seq(-4,5,length=2),labels=seq(-4,5,length=2)) +
  xlab("Pred B Scores") +
  ylab("Pred N Scores") +
  coord_fixed() +
  theme(aspect.ratio=1) +
  theme1 + 
	theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(5,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,5,0,0)),
        axis.line = element_line(colour = "black", size = 1))
dev.off()

# -- Plot the summary of leave-one-out predicted CCA corr vs observed
ccaall <- cc(behpc, gbcPSD_180pS)
predDF <- data.frame("CV"=1:5,"r"=diag(cor(predX,predY)), "CanCor"=ccaall$cor)
pdf("CCA_LeaveOneOut_PredScores_Summary.pdf", height=6, width=8)
ggplot(predDF, aes(x=CV,y=r)) +
	#geom_hline(yintercept=1, col="red", lwd=1.5) +
	geom_segment(aes(x=c(0.5,1.5,2.5,3.5,4.5), y=CanCor, xend=c(1.5,2.5,3.5,4.5,5.5), yend=CanCor), size=1.5,col="red") +
	geom_col() + 
  	scale_alpha_manual(values = c("1"="1","2"="0.9","3"="0.8", "4"="0.7", "5"="0.6")) +
  	scale_y_continuous(expand=c(0.0001,0), limits=c(0.0, 1.1),breaks=seq(0,1,0.5),labels=format(seq(0,1,0.5), nsmall=1)) +
  	scale_x_continuous(expand=c(0.000001,0),limits=c(0.45,5.5),breaks=seq(1,5,by=1),labels=1:5) +
  	xlab("Canonical Variate") +
  	ylab("r-value") + 
  	#xlab("Value of k") + 
  	theme1 +
	theme(axis.text.x = element_text(size = 25, color="black"), 
        axis.title.x = element_text(size = 25, margin=margin(5,0,0,0)),
        axis.text.y = element_text(size = 25, color="black",angle=90,hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 25,margin=margin(0,5,0,0)),
        axis.line = element_line(colour = "black", size = 1))
dev.off()