rm(list = ls())
setwd("../fig1_symptom_pca_analyses")
library(ggplot2)
library(scales)
library(matrixStats)
library(plot3D)
library(rgl)
library(ggridges)
# library(reshape2)
# library(fmsb)
# library(gplots)
# library(Hmisc)
# library(grid)

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

theme2<- theme(panel.border = element_rect(size = 2, colour = 'NA', fill='NA'), 
               axis.line = element_line(colour = "black", size = 2),
               axis.text.x = element_text(size = 40, color="black",margin=margin(10,20,0,0)),
               axis.title.x = element_blank(),
               axis.text.y = element_blank(),
               axis.title.y = element_text(size = 40,margin=margin(0,5,0,0)),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               plot.background = element_rect(fill = "transparent", colour = NA),
               panel.background = element_rect(fill = "transparent", colour = NA),
               legend.position="none",
               axis.ticks = element_line(colour = "black", size=2),
               axis.ticks.length = unit(.4, "cm"),
               axis.ticks.y = element_blank())               

# -- Functions
radarchartvar<-function (df, latcol="black", latlty=1, axistype = 0, seg = 4, pty = 16, pcol = 1:8, plty = 1:6, 
                         plwd = 1, pdensity = NULL, pangle = 45, pfcol = NA, cglty = 3, 
                         cglwd = 1, cglcol = "navy", axislabcol = "blue", title = "", 
                         maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlabels = NULL, 
                         vlcex = NULL, caxislabels = NULL, calcex = NULL, paxislabels = NULL, 
                         seglty = 1, seglwd = 3, segcol = "black", segmaxcol="black",segmincol="black",palcex = NULL, ...)
{
  if (!is.data.frame(df)) {
    cat("The data must be given as dataframe.\n")
    return()
  }
  if ((n <- length(df)) < 3) {
    cat("The number of variables must be 3 or more.\n")
    return()
  }
  if (maxmin == FALSE) {
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  plot(c(-1.2, 1.2), c(-1.2, 1.2), type = "n", frame.plot = FALSE, 
       axes = FALSE, xlab = "", ylab = "", main = title, asp = 1, 
       ...)
  theta <- seq(90, 450, length = n + 1) * pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  CGap <- ifelse(centerzero, 0, 1)
  for (i in 0:seg) {
    polygon(xx * (i + CGap)/(seg + CGap), yy * (i + CGap)/(seg + CGap), lty = seglty, lwd = seglwd, border = segcol)
    if (axistype == 1 | axistype == 3) 
      CAXISLABELS <- paste(i/seg * 100, "(%)")
    if (axistype == 4 | axistype == 5) 
      CAXISLABELS <- sprintf("%3.2f", i/seg)
    if (!is.null(caxislabels) & (i < length(caxislabels))) 
      CAXISLABELS <- caxislabels[i + 1]
    if (axistype == 1 | axistype == 3 | axistype == 4 | axistype == 
        5) {
      if (is.null(calcex)) 
        text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS, 
             col = axislabcol)
      else text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS, 
                col = axislabcol, cex = calcex)
    }
  }
  #polygon(xx * ((seg/2) + CGap)/(seg + CGap), yy * ((seg/2) + CGap)/(seg + CGap), lty = 1, lwd = 4.5, border = latcol)
  if (centerzero) {
    arrows(0, 0, xx * 1, yy * 1, lwd = cglwd, lty = cglty, 
           length = 0, col = cglcol)
  }
  else {
    arrows(xx/(seg + CGap), yy/(seg + CGap), xx * 1, yy * 
             1, lwd = cglwd, lty = cglty, length = 0, col = cglcol)
  }
  PAXISLABELS <- df[1, 1:n]
  if (!is.null(paxislabels)) 
    PAXISLABELS <- paxislabels
  if (axistype == 2 | axistype == 3 | axistype == 5) {
    if (is.null(palcex)) 
      text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol)
    else text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol, 
              cex = palcex)
  }
  VLABELS <- colnames(df)
  if (!is.null(vlabels)) 
    VLABELS <- vlabels
  if (is.null(vlcex)) 
    text(xx * 1.2, yy * 1.2, VLABELS)
  else text(xx * 1.2, yy * 1.2, VLABELS, cex = vlcex)
  series <- length(df[[1]])
  SX <- series - 2
  if (length(pty) < SX) {
    ptys <- rep(pty, SX)
  }
  else {
    ptys <- pty
  }
  if (length(pcol) < SX) {
    pcols <- rep(pcol, SX)
  }
  else {
    pcols <- pcol
  }
  if (length(plty) < SX) {
    pltys <- rep(plty, SX)
  }
  else {
    pltys <- plty
  }
  if (length(plwd) < SX) {
    plwds <- rep(plwd, SX)
  }
  else {
    plwds <- plwd
  }
  if (length(pdensity) < SX) {
    pdensities <- rep(pdensity, SX)
  }
  else {
    pdensities <- pdensity
  }
  if (length(pangle) < SX) {
    pangles <- rep(pangle, SX)
  }
  else {
    pangles <- pangle
  }
  if (length(pfcol) < SX) {
    pfcols <- rep(pfcol, SX)
  }
  else {
    pfcols <- pfcol
  }
  for (i in 3:series) {
    xxs <- xx
    yys <- yy
    scale <- CGap/(seg + CGap) + (df[i, ] - df[2, ])/(df[1, 
                                                         ] - df[2, ]) * seg/(seg + CGap)
    if (sum(!is.na(df[i, ])) < 3) {
      cat(sprintf("[DATA NOT ENOUGH] at %d\n%g\n", i, df[i, 
                                                         ]))
    }
    else {
      for (j in 1:n) {
        if (is.na(df[i, j])) {
          if (na.itp) {
            left <- ifelse(j > 1, j - 1, n)
            while (is.na(df[i, left])) {
              left <- ifelse(left > 1, left - 1, n)
            }
            right <- ifelse(j < n, j + 1, 1)
            while (is.na(df[i, right])) {
              right <- ifelse(right < n, right + 1, 1)
            }
            xxleft <- xx[left] * CGap/(seg + CGap) + 
              xx[left] * (df[i, left] - df[2, left])/(df[1, 
                                                         left] - df[2, left]) * seg/(seg + CGap)
            yyleft <- yy[left] * CGap/(seg + CGap) + 
              yy[left] * (df[i, left] - df[2, left])/(df[1, 
                                                         left] - df[2, left]) * seg/(seg + CGap)
            xxright <- xx[right] * CGap/(seg + CGap) + 
              xx[right] * (df[i, right] - df[2, right])/(df[1, 
                                                            right] - df[2, right]) * seg/(seg + CGap)
            yyright <- yy[right] * CGap/(seg + CGap) + 
              yy[right] * (df[i, right] - df[2, right])/(df[1, 
                                                            right] - df[2, right]) * seg/(seg + CGap)
            if (xxleft > xxright) {
              xxtmp <- xxleft
              yytmp <- yyleft
              xxleft <- xxright
              yyleft <- yyright
              xxright <- xxtmp
              yyright <- yytmp
            }
            xxs[j] <- xx[j] * (yyleft * xxright - yyright * 
                                 xxleft)/(yy[j] * (xxright - xxleft) - xx[j] * 
                                            (yyright - yyleft))
            yys[j] <- (yy[j]/xx[j]) * xxs[j]
          }
          else {
            xxs[j] <- 0
            yys[j] <- 0
          }
        }
        else {
          xxs[j] <- xx[j] * CGap/(seg + CGap) + xx[j] * 
            (df[i, j] - df[2, j])/(df[1, j] - df[2, j]) * 
            seg/(seg + CGap)
          yys[j] <- yy[j] * CGap/(seg + CGap) + yy[j] * 
            (df[i, j] - df[2, j])/(df[1, j] - df[2, j]) * 
            seg/(seg + CGap)
        }
      }
      if (is.null(pdensities)) {
        polygon(xxs, yys, lty = pltys[i - 2], lwd = plwds[i - 
                                                            2], border = pcols[i - 2], col = pfcols[i - 
                                                                                                      2])
      }
      else {
        polygon(xxs, yys, lty = pltys[i - 2], lwd = plwds[i - 
                                                            2], border = pcols[i - 2], density = pdensities[i - 
                                                                                                              2], angle = pangles[i - 2], col = pfcols[i - 
                                                                                                                                                         2])
      }
      points(xx * scale, yy * scale, pch = ptys[i - 2], 
             col = pcols[i - 2])
    }
  }
  polygon(xx * ((seg/2) + CGap)/(seg + CGap), yy * ((seg/2) + CGap)/(seg + CGap), lty = latlty, lwd = 4.5, border = latcol)
  polygon(xx * (0 + CGap)/(seg + CGap), yy * (0 + CGap)/(seg + CGap), lty = seglty, lwd = seglwd, border = segmincol)
  polygon(xx * (seg + CGap)/(seg + CGap), yy * (seg + CGap)/(seg + CGap), lty = seglty, lwd = seglwd, border = segmaxcol)
}
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

# -- Read data
beh_raw <-  read.table("../data_symptoms/BSNIP_PSD_N436_Behavior.dat", sep="\t", header=TRUE)
beh_PSD <- beh_raw[1:436, 2:37]

# -- Compute PCA and significance
pcaA = prcomp(beh_PSD, scale. = TRUE)
write.table(pcaA$rotation, file = "BSNIP_N436_BACSPANSS_PCArotations.txt")

prop_var <- data.frame(cbind((1:36), pcaA$sdev^2/sum(pcaA$sdev^2)))
names(prop_var) <- c("PC", "Variance")

R=5000
s=36
signif <- sign.pc(beh_PSD,R=R,s=s)
pve.perm <- signif$pve.perm
pve <- signif$pve
pval <- signif$pval # First 5 PCs are significant

# -- Compute 0.95 chance
crit95 <- matrix(0,nrow=1,ncol=s)
for(i in 1:s){
  crit95[i] <- sort(pve.perm[1:R, i])[250]
}
crit=data.frame(PC=c(1:36),crit=c(crit95))

# -- Plots
# -- Screeplot of variance explained by each PC
pdf("Symptom_PCA_Screeplot.pdf", width=10, height=6)
ggplot(data=prop_var, aes(x=PC, y=Variance)) +
  geom_line(data=crit, aes(x=PC,y=crit), lty=2, col="dark red",lwd=1.2) +
  geom_line(lwd=1.2)+
  geom_point(size=4,col="black",pch=21,bg="grey35") +
  geom_point(data=prop_var[1:5,1:2], aes(size=1/PC),col="black",pch=21,bg="green") +
  scale_size_area(name=prop_var[1:5,1:2], max_size=15) +
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
dev.off()

# -- Pie chart of variance attributed to significant vs non-significant PCs
pve_prop <- data.frame(PCs=c("1-5","6-36"),"Total Variance"=c(sum(pve[1:5]),sum(pve[6:36])))

pdf("Symptom_PCA_ProportionVarianceExplained_Piechart.pdf", width=10, height=6)
ggplot(data=pve_prop,aes(x=factor(1),y=Total.Variance,fill=factor(PCs))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta="y") +
  scale_fill_manual(values = c("1-5"="green","6-36"="grey35")) +
  geom_text(aes(y = Total.Variance/2 + c(cumsum(Total.Variance)[-length(Total.Variance)],0), 
                label = percent(Total.Variance)), size=8, col=c("black","white")) +
  theme1 +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line=element_blank(),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28),
        legend.box.spacing=unit(0.01,"inch"),
        legend.box.margin=margin(0,-15,0,0),
        legend.margin=margin(0,0,0,0),
        legend.position="left",
        panel.background = element_rect(fill = "transparent", colour = NA)) +
  guides(fill=guide_legend(keywidth=0.7,keyheight=0.7, default.unit="inch",title="PCs"))
dev.off()

scores_PSD = as.data.frame(cbind(beh_raw[1],pcaA$x))

# -- Project CON subjects into PSD symptom PC space
loadings <- pcaA$rotation
loadingsm <- as.matrix(loadings)
beh_CON <- read.table("../data_symptoms/BSNIP_CON_N202_behavior.dat", sep="\t", header=TRUE)
beh_CONm <-matrix(unlist(beh_CON[1:202,2:37]), ncol=36, byrow = FALSE)


PROBbehZMean <- colMeans(beh_raw[-1])
PROBbehZSD <- colSds(as.matrix(beh_raw[-1]))
CONbehZ_ZbyPROB <- sweep(sweep(beh_CONm, 2, PROBbehZMean, "-"), 2, PROBbehZSD,"/") # Z-score CON relative to PROB
CON_ZbyPROB_scores <- as.matrix(CONbehZ_ZbyPROB) %*% loadingsm
CON_ZbyPROB_scores_mean=colMeans(CON_ZbyPROB_scores)
CON_ZbyPROB_scores_sd=colSds(CON_ZbyPROB_scores)
CON_ZbyPROB_scores_Z <- scale(CON_ZbyPROB_scores)[,]

beh_PSD_mean <- colMeans(beh_PSD)
beh_PSD_sd <- colSds(as.matrix(beh_PSD))
beh_CON_ZbyPSD <- sweep(sweep(beh_CONm, 2, beh_PSD_mean, "-"), 2, beh_PSD_sd,"/") # Z-score CON relative to PSD

CON_ZbyPSD_scores <- as.matrix(beh_CON_ZbyPSD) %*% loadingsm # project CON into PSD PCA space
CON_ZbyPSD_scores_mean=colMeans(CON_ZbyPSD_scores) # mean of CON PC scores
CON_ZbyPSD_scores_sd=colSds(CON_ZbyPSD_scores) # sd of CON PC scores
CON_ZbyPSD_scores_Z <- scale(CON_ZbyPSD_scores)[,] # Z-score CON scores

# -- Make a dataframe with the groups
scoresAll=NULL
scoresAll$Group <- factor(scoresAll$Group, levels = c("CON","PSD","BPP","SADP","SCZP"))
PSDGscores_ZbyCONscores <- cbind(scores_PSD[1],sweep(sweep(scores_PSD[-1],2,CON_ZbyPSD_scores_mean),2,CON_ZbyPSD_scores_sd,"/")) # Standardize PSD by CON scores
scores_PSD_ZbyCONscores <-data.frame(Group=c("PSD"),PSDGscores_ZbyCONscores[-1])
CONscores_ZbyCONscores <- data.frame(Group=as.factor(c("CON")), CON_ZbyPSD_scores_Z)
scoresAll <- rbind(CONscores_ZbyCONscores,scores_PSD_ZbyCONscores,PSDGscores_ZbyCONscores)

# -- Write the scores to file
write.table(scores_PSD, file = "BSNIP_PSD436_AllPCScores.txt")
write.table(CON_ZbyPSD_scores, file = "BSNIP_CON202_AllPCScores.txt")

# -- Select the componenents, set up for plotting biplots
scoresPSD_PC123 <- pcaA$x[,1:3]
scoresCON_PC123 <- CON_ZbyPSD_scores[,1:3]
PC123 <-rbind(scoresPSD_PC123,scoresCON_PC123)
PC123 <- cbind(GCols=rbind(beh_raw[1],CONscores_ZbyCONscores[1]),PC123)
PC123$GCols[PC123$Group=='BPP']<-'#feb24c'
PC123$GCols[PC123$Group=='SADP']<-'#fc4e2a'
PC123$GCols[PC123$Group=='SCZP']<-'#800026'
PC123$GCols[PC123$Group=='CON']<-'black'

scoresAll$Group <- factor(scoresAll$Group, levels = c("SCZP","SADP","BPP","PSD","CON"))
scoresAll_means <- aggregate(scoresAll,by=list(scoresAll$Group), FUN="mean")

# -- Plot biplots
colours=c("#004E2E","#005F38","#006E45","#008953","#009E5D","#00B669",
          "#441978","#511D90","#6124AF","#6F26C3","#7C2CE2","#8E35F5","#A351FF",
          "#354252","#334862","#37506E","#406084","#41658F","#446C9C","#4572A6",
          "#5B013B","#6A0144","#7F0154","#950264","#AD016D","#C2017A","#CC0285","#D3038C","#DF0293","#ED039F",
          "#F702A6","#FF3DB9","#FF5ABE","#FF7FD4","#FFA1E4","#FEB5EB")

z <- c(rbind(colours, colours))
coords <- NULL

for (i in 1:nrow(pcaA$rotation)) {
  coords <- rbind(coords, rbind(c(0,0,0),16*pcaA$rotation[i,1:3]))
}
coordsCog <- (35*colMeans(pcaA$rotation[1:6,1:3]))
coordsPos <- (35*colMeans(pcaA$rotation[7:13,1:3]))
coordsNeg <- (35*colMeans(pcaA$rotation[14:20,1:3]))
coordsGen <- (35*colMeans(pcaA$rotation[21:36,1:3]))

par3d(cex=5, windowRect = c(20, 30, 800, 800))
rgl.viewpoint(  zoom = 1.2 )
plot3d(PC123[2:4], aspect=TRUE,
       xlab=expression("PC1"),ylab=expression("PC2"),zlab=expression("PC3"), 
       type="p", axes=F, size=15, col=PC123$GCols, alpha=0.7,
       xlim=c(-10,10), ylim=c(-10,10), zlim=c(-10,10))
lines3d(coords, col=z, lwd=4)
axis3d(edge="x--", lwd=3, axes.len=1, labels=c("-10","",0,"","10"),cex=2)
axis3d(edge="y--", lwd=3, axes.len=1, labels=c("","",0,"","10"),cex=2)
axis3d(edge="z", lwd=3, axes.len=1, labels=c("","",0,"","10"),cex=2)
grid3d("x")
grid3d("y")
grid3d("z")
view3d(userMatrix=rotationMatrix(2*pi * 1, 1, -1, -1))

PC123means <- aggregate(PC123, by=list(PC123$Group), FUN="mean")[,-2]
names(PC123means)[1] <- c("Group")
PC123means$GCols[PC123means$Group=='BPP']<-'#feb24c'
PC123means$GCols[PC123means$Group=='SADP']<-'#fc4e2a'
PC123means$GCols[PC123means$Group=='SCZP']<-'#800026'
PC123means$GCols[PC123means$Group=='CON']<-'black'

# -- Centroid view
par3d(cex=5, windowRect = c(20, 30, 800, 800))
rgl.viewpoint(  zoom = 1.2 )
plot3d(PC123means[2:4], aspect=TRUE,
       xlab="",ylab="",zlab="", 
       type="s", axes=F, size=5, col=PC123means$GCols, alpha=0.6,
       xlim=c(-10,10), ylim=c(-10,10), zlim=c(-10,10))
arrow3d(p0=c(0,0,0),p1=coordsCog, col="green",  lwd=50,alpha=1, type="rotation",barblen=0.08,width=0.6,n=3)
arrow3d(p0=c(0,0,0),p1=coordsPos, col="#8A40F5", lwd=50,alpha=1, type="rotation",barblen=0.08,width=0.6,n=3)
arrow3d(p0=c(0,0,0),p1=coordsNeg, col="#00B1FF",   lwd=50,alpha=1, type="rotation",barblen=0.08,width=0.6,n=3)
arrow3d(p0=c(0,0,0),p1=coordsGen, col="pink",   lwd=50,alpha=1, type="rotation",barblen=0.08,width=0.6,n=3)
axis3d(edge="x--", lwd=3, axes.len=1, labels=c("-10","",0,"","10"),cex=2)
axis3d(edge="y--", lwd=3, axes.len=1, labels=c("","",0,"","10"),cex=2)
axis3d(edge="z", lwd=3, axes.len=1, labels=c("","","",""),cex=2)
grid3d("x--")
grid3d("y")
grid3d("z++")
view3d(userMatrix=rotationMatrix(2*pi * 1, 1, -1, -1))

# -- Plot ridge plots
pdf("Symptom_PCA_RidgelinePlots_PC1.pdf", width=10, height=6)
ggplot(scoresAll, aes(x=PC1, y=Group,group=Group)) + 
  geom_density_ridges(aes(fill=Group), lwd=2, rel_min_height = 0.01) +
  scale_fill_manual(values = c("CON"="white","PSD"="black","BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-32,32), breaks=c(-30,-20,-10,0,10,20,30),labels=c(-30,"","",0,"","",30)) +
  ylab("PC1") +
  xlab("") +
  theme2 +
  geom_segment(aes(x = scoresAll_means$PC1[1], y=1.0,xend = scoresAll_means$PC1[1], yend = 2),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = scoresAll_means$PC1[2], y=2.0,xend = scoresAll_means$PC1[2], yend = 3),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = scoresAll_means$PC1[3], y=3.0,xend = scoresAll_means$PC1[3], yend = 4),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = scoresAll_means$PC1[4], y=4.0,xend = scoresAll_means$PC1[4], yend = 5),size=2, linetype = 1,color = "white") +
  geom_vline(xintercept=scoresAll_means$PC1[5], lty=2,col="grey",lwd=3)
dev.off()

pdf("Symptom_PCA_RidgelinePlots_PC2.pdf", width=10, height=6)
ggplot(scoresAll, aes(x=PC2, y=Group,group=Group)) + 
  geom_density_ridges(aes(fill=Group), lwd=2, rel_min_height = 0.01) +
  scale_fill_manual(values = c("CON"="white","PSD"="black","BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-5,5), breaks=c(-4,-2,0,2,4), labels=c(-4,"",0,"",4)) +
  ylab("PC2") +
  xlab("") +
  theme2 +
  geom_segment(aes(x =  scoresAll_means$PC2[1], y=1.00,xend = scoresAll_means$PC2[1], yend = 2.0),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  scoresAll_means$PC2[2], y=2.17,xend = scoresAll_means$PC2[2], yend = 3.0),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  scoresAll_means$PC2[3], y=3.25,xend = scoresAll_means$PC2[3], yend = 4.0),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  scoresAll_means$PC2[4], y=4.30,xend = scoresAll_means$PC2[4], yend = 5.0),size=2, linetype = 1,color = "white") +
  geom_vline(xintercept=scoresAll_means$PC2[5], lty=2,col="grey",lwd=3)
dev.off()

pdf("Symptom_PCA_RidgelinePlots_PC3.pdf", width=10, height=6)
ggplot(scoresAll, aes(x=PC3, y=Group,group=Group)) + 
  geom_density_ridges(aes(fill=Group), lwd=2, rel_min_height = 0.005) +
  scale_fill_manual(values = c("CON"="white","PSD"="black","BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-11,11), breaks=c(-10,-5,0,5,10), labels=c(-10,"",0,"",10)) +
  ylab("PC3") +
  xlab("") +
  theme2 +
  geom_segment(aes(x =  scoresAll_means$PC3[1], y=1.0,xend = scoresAll_means$PC3[1], yend = 2),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  scoresAll_means$PC3[2], y=2.0,xend = scoresAll_means$PC3[2], yend = 3),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  scoresAll_means$PC3[3], y=3.0,xend = scoresAll_means$PC3[3], yend = 4),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  scoresAll_means$PC3[4], y=4.0,xend = scoresAll_means$PC3[4], yend = 5),size=2, linetype = 1,color = "white") +
  geom_vline(xintercept=scoresAll_means$PC3[5], lty=2,col="grey",lwd=3)
dev.off()

pdf("Symptom_PCA_RidgelinePlots_PC4.pdf", width=10, height=6)
ggplot(scoresAll, aes(x=PC4, y=Group,group=Group)) + 
  geom_density_ridges(aes(fill=Group), lwd=2, rel_min_height = 0.01) +
  scale_fill_manual(values = c("CON"="white","PSD"="black","BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-11,11), breaks=c(-10,-5,0,5,10), labels=c(-10,"",0,"",10)) +
  ylab("PC4") +
  xlab("") +
  theme2 +
  geom_segment(aes(x =  scoresAll_means$PC4[1], y=1.0,xend = scoresAll_means$PC4[1], yend = 2),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  scoresAll_means$PC4[2], y=2.0,xend = scoresAll_means$PC4[2], yend = 3),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  scoresAll_means$PC4[3], y=3.0,xend = scoresAll_means$PC4[3], yend = 4),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  scoresAll_means$PC4[4], y=4.0,xend = scoresAll_means$PC4[4], yend = 5),size=2, linetype = 1,color = "white") +
  geom_vline(xintercept=scoresAll_means$PC4[5], lty=2,col="grey",lwd=3)
dev.off()

pdf("Symptom_PCA_RidgelinePlots_PC5.pdf", width=10, height=6)
ggplot(scoresAll, aes(x=PC5, y=Group,group=Group)) + 
  geom_density_ridges(aes(fill=Group), lwd=2, rel_min_height = 0.015) +
  scale_fill_manual(values = c("CON"="white","PSD"="black","BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-9,9), breaks=c(-8,-4,0,4,8), labels=c(-8,"",0,"",8)) +
  ylab("PC5") +
  xlab("") +
  theme2 +
  geom_segment(aes(x =  scoresAll_means$PC5[1], y=1.0,xend = scoresAll_means$PC5[1], yend = 2.0),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  scoresAll_means$PC5[2], y=2.0,xend = scoresAll_means$PC5[2], yend = 3.0),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  scoresAll_means$PC5[3], y=3.0,xend = scoresAll_means$PC5[3], yend = 4.0),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  scoresAll_means$PC5[4], y=4.2,xend = scoresAll_means$PC5[4], yend = 5.0),size=2, linetype = 1,color = "white") +
  geom_vline(xintercept=scoresAll_means$PC5[5], lty=2,col="grey",lwd=3)
dev.off()

# -- Plot radarplots of PC loadings
assign("PCloadingsRadarplot",data.frame(rbind(rep(0.5,36), rep(-0.5,36),t(pcaA$rotation[,5:1]))))
Sxcols=c("#004E2E","#005F38","#006E45","#008953","#009E5D","#00B669",
          "#441978","#511D90","#6124AF","#6F26C3","#7C2CE2","#8E35F5","#A351FF",
          "#354252","#334862","#37506E","#406084","#41658F","#446C9C","#4572A6",
          "#5B013B","#6A0144","#7F0154","#950264","#AD016D","#C2017A","#CC0285","#D3038C","#DF0293","#ED039F",
          "#F702A6","#FF3DB9","#FF5ABE","#FF7FD4","#FFA1E4","#FEB5EB")
pdf("Symptom_PCA_Loadings_Radarplot.pdf", width=6, height=6)
radarchartvar(PCloadingsRadarplot, latlty=1, latcol="darkred", axistype=0, seglty = 1, seglwd = 2, segmincol=c("grey"), segmaxcol=c("grey"), segcol=c("white"), cglty=1, cglwd=2, vlabels="", cglcol=c(Sxcols), pcol=c("grey30","grey40","grey50","grey60","grey70"), plwd=7, seg=6, plty=c(1,1,1,1,1))
dev.off()




