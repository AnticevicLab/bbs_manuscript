rm(list = ls())
setwd("../fig2_symptom_pca_stability_analyses")
library(matrixStats)
library(ggplot2)
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
mtsortdiag <- function(mt){
        mtf <- as.matrix(mt)
        mtb <- as.matrix(mt)

        ncs <- ncol(mt)-1

        for (nc in 1:ncs){
                if ((sum(abs(diag(mtf[nc:(nc+1),nc:(nc+1)])) < abs(diag(mtf[nc:(nc+1),(nc+1):nc]))) == 0) && (sum(abs(diag(mtf[nc:(nc+1),nc:(nc+1)])) < abs(diag(mtf[(nc+1):nc,nc:(nc+1)]))) == 0)){
                        print("diagonal is maximal")
                } else {
                        if ( sum(abs(diag(mtf[nc:(nc+1),nc:(nc+1)]))) < sum(abs(diag(mtf[nc:(nc+1),(nc+1):nc]))) ){
                                print("diagonal is not maximal")
                                mtf[,nc:(nc+1)] <- mtf[,(nc+1):nc]
                                names(mtf[,nc:(nc+1)]) <- names(mtf[,(nc+1):nc])
                        }
                }
        }

        for (nc in ncs:1){
                if ((sum(abs(diag(mtb[nc:(nc+1),nc:(nc+1)])) < abs(diag(mtb[nc:(nc+1),(nc+1):nc]))) == 0) && (sum(abs(diag(mtb[nc:(nc+1),nc:(nc+1)])) < abs(diag(mtb[(nc+1):nc,nc:(nc+1)]))) == 0)){
                        print("diagonal is maximal")
                } else {
                        if ( sum(abs(diag(mtb[nc:(nc+1),nc:(nc+1)]))) < sum(abs(diag(mtb[nc:(nc+1),(nc+1):nc]))) ){
                                print("diagonal is not maximal")
                                mtb[,nc:(nc+1)] <- mtb[,(nc+1):nc]
                                names(mtb[,nc:(nc+1)]) <- names(mtb[,(nc+1):nc])
                        }
        }

        if ((sum(abs(diag(mtf)))>sum(abs(diag(mtb))))){
                return(mtf)
        } else {
                return(mtb)
        }
        }
}

test.pc.signif <- function(x,R=5000,s=36,...){
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
  # calculate the p-values
  pval<-apply(t(pve.perm)>pve,1,sum)/R
  return(list(pve=pve,pval=pval, pve.perm=pve.perm))
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

# -- Load data
beh_raw = read.table("../data_symptoms/BSNIP_PSD_N436_Behavior.dat", sep="\t", header=TRUE)
beh=beh_raw[1:436,2:37]
Group = as.factor(beh_raw[1:436,1])
# -- Set parameters
nRep <- 1000 # how many times to repeat the split-half PCA?
nB=ncol(beh)
R=1000
s=nB

# -- Initiate some variables to store the subject indices for each run, when assigning stratified split-half samples
H1shuffind <- array(NA,dim=c(217,nB,nRep))
H2shuffind <- array(NA,dim=c(219,nB,nRep))
H1indices <- matrix(NA, 217, nRep)
H2indices <- matrix(NA, 219, nRep)

# -- Assign stratified split-half samples
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
	
	H1shuffind[,,shuff] <- as.matrix(beh[h1.mixed,])
	H2shuffind[,,shuff] <- as.matrix(beh[-h1.mixed,])
}

# -- Initiate the variables in which we will store the outputs from each run of split-half PCA 
H1rot <- array(NA, dim=c(36,5,nRep))
H2rot <- array(NA, dim=c(36,5,nRep))
H1predScores <- array(NA, dim=c(217,5,nRep))
H2predScores <- array(NA, dim=c(219,5,nRep))
shdat <- NULL
H1Scores <- array(NA, dim=c(217,5,nRep))
H2Scores <- array(NA, dim=c(219,5,nRep))

# -- Perform the split-half CCA nRep times
for (shuff in 1:nRep){
        # 
	H1beh <- data.frame(H1shuffind[,,shuff])
	H2beh <- data.frame(H2shuffind[,,shuff])
	
	for (half in c("H1","H2")){
                # -- Test for significance etc.
		assign(paste0("PCA_PermTest_",half) , test.pc.signif(get(paste0(half,"beh")),R=R,s=s))
		assign(paste0("permutedPVE_",half) , eval(parse(text=paste0("PCA_PermTest_",half,"$pve.perm"))))
		assign(paste0("obsPVE_",half) , eval(parse(text=paste0("PCA_PermTest_",half,"$pve"))))
		assign(paste0("permutedpval_",half) ,  eval(parse(text=paste0("PCA_PermTest_",half,"$pval"))))
		
		criti <- floor(0.05 * R)
		if (criti < 1){
		 stop("Please specify a greater nperm to calculate p<0.05", call.=FALSE)
		}
		crit95 <- matrix(0,nrow=1,ncol=s)
		for(i in 1:s){
		 crit95[i] <- sort(permutedPVE_H1[1:R, i])[criti]
		}
		crit <- data.frame(PC=c(1:nB),crit=c(crit95))
		ncrit <- max(which(get(paste0("permutedpval_",half)) < 0.05))
		
                # -- Run PCA
		assign(paste0("pca",half), prcomp(get(paste0(half,"beh")), scale. = TRUE))
	}
        # -- Correlate the loadings of the first 5 PCs between H1 and H2 PCAs
	mt <- cor(pcaH1$rotation[,1:5],pcaH2$rotation[,1:5])

        # -- Sometimes, the order of the PCs is slightly switched between the two halves, because the PVE is quite similar. 
        # -- Hence we sort the order by maximizing the diagonal to get corresponding PC pairs
	mtsort <- mtsortdiag(mt)

	# -- Save the PCA scores so we can compare them later
        H1Scores[,,shuff] <- pcaH1$x[,1:5]
	H2Scores[,,shuff] <- pcaH2$x[,1:5]

        # -- Save the PCA rotations (i.e. loadings) so we can compare them later
	H1rot[,,shuff] <- pcaH1$rotation[,colnames(mtsort)]
	H2rot[,,shuff] <- pcaH2$rotation[,colnames(mtsort)]
	
        # -- Use the PCA loadings from H1 to predict H2 scores, and vice versa
	H2predScores[,,shuff] <- as.matrix(scale(H2beh)[,]) %*% as.matrix(pcaH1$rotation[,colnames(mtsort)])
	H1predScores[,,shuff] <- as.matrix(scale(H1beh)[,]) %*% as.matrix(pcaH2$rotation[,colnames(mtsort)])
	
	shdat <- c(shdat, diag(mtsort))
}
# -- Save the indices for each half for each run, so we know which subjects were in which half -- for mapping to neural features later
for (shuff in 1:nRep){
	for (half in c("H1","H2")){
		write.table(get(paste0(half,"indices"))[,shuff], file=paste0("BSNIP_StratifiedPRB_",half,"Indices_Run",shuff,".txt"), sep="", col.names=FALSE, row.names=FALSE)
}}

# -- Save the predicted scores for each run
for (shuff in 1:nRep){
	for (half in c("H1","H2")){
		for (pc in 1:5){
			write.table(get(paste0(half,"predScores"))[,pc,shuff], file=paste0("BSNIP_StratifiedPRB_",half,"PredScores","_dPC",pc,"Run",shuff,".csv"), col.names=FALSE, row.names=FALSE)
	}}}

# -- Save the observed scores for each half for each run
for (shuff in 1:nRep){
	for (half in c("H1","H2")){
		for (pc in 1:5){
			write.table(get(paste0(half,"Scores"))[,pc,shuff], file=paste0("BSNIP_StratifiedPRB_",half,"Scores","_dPC",pc,"Run",shuff,".csv"), col.names=FALSE, row.names=FALSE)
}}}

# -- Organize the dataframef for plotting
shRepDF <- data.frame("PC"=1:5, "R"=abs(shdat))
DFsummshRep <- summarySE(shRepDF, measurevar="R", groupvars=c("PC"))

# -- Plot the mean and SE r-value between H1 and H2 across all runs
pdf("BSNIP_SymptomPCA_SplitHalfComparison.pdf", width=8, height=6)
ggplot(DFsummshRep, aes(x=PC, y=R, alpha=factor(PC))) +
  geom_hline(yintercept=1, col="red", lwd=1.5) +
  stat_summary(position="dodge", geom="bar", width=1, fun = "mean", fill=c("#485C77")) +
  geom_errorbar(aes(x=PC, ymin=R-se, ymax=R+se), width=.5, position=position_dodge(0.9), lwd=1.5) +
  scale_alpha_manual(values = c("1"=1,"2"=0.9,"3"=0.8,"4"=0.7, "5"=0.6)) +
  scale_y_continuous(expand=c(0.0001,0), limits=c(0.0, 1.1),breaks=seq(0,1,0.5),labels=format(seq(0,1,0.5), nsmall=1)) +
  scale_x_continuous(expand=c(0.000001,0),limits=c(0.45,5.5),breaks=seq(1,5,by=1),labels=c("","","","","")) +
  ylab("r-value (Mean and SE)") + 
  xlab("PC") + 
  theme1 +
  theme(axis.text.x = element_text(size = 25, color="black"),
        axis.title.x = element_text(size = 30, margin=margin(15,0,0,0)),
        axis.title.y = element_text(size = 25,margin=margin(0,10,0,0)),
        axis.ticks= element_blank())
dev.off()