rm(list = ls())
setwd("../fig5_cca_analyses")

library(ggplot2)
library(GGally)
library(CCA)
library(CCP)
library(matrixStats)
library(ggridges)
library(fmsb)

# -- Functions
ccaVarperm <- function (X, Y, nboot = 999) 
  {
    if (nrow(Y) != nrow(X)) 
      stop(" Function p.perm: Input arrays must not have unequal number of rows.\n")
    N <- nrow(Y)
    ind <- 1:N
    p = dim(X)[2]
    q = dim(Y)[2]
    minpq = min(p, q)
    cc0 <- cc(X, Y) # run observed CCA
    rho0 <- cc0$cor # observed canonical correlations
    stat0rho <- as.matrix(rho0) # observed correlations
    stat0X <- colMeans(cc0$scores$corr.X.xscores^2)*rho0^2 # observed redundancy of X (var explained by Y)
    stat0Y <- colMeans(cc0$scores$corr.Y.yscores^2)*rho0^2 # observed redundancy of Y (var explained by X)
    statrho <- matrix(nrow = minpq, ncol = nboot) # initialise for correlations
    statX <- matrix(nrow = minpq, ncol = nboot) # initialise
    statY <- matrix(nrow = minpq, ncol = nboot) # initialise
    for (i in 1:nboot) { # run permutations with shuffled Y
      ind.mixed <- sample(ind, size = N, replace = FALSE)
      cci <- cc(X, Y[ind.mixed, ])
      statX[,i] <- colMeans(cci$scores$corr.X.xscores^2)*(cci$cor^2)
      statY[,i] <- colMeans(cci$scores$corr.Y.yscores^2)*(cci$cor^2)
      statrho[,i] <- cci$cor
    }
    pvalrho <- numeric(minpq)
    pvalX <- numeric(minpq)
    pvalY <- numeric(minpq)
    for (i in 1:minpq) {
      pvalrho[i] <- mean(statrho[i,]>= stat0rho[i]) # how many random shuffles exceeded the observed rho?
      pvalX[i] <- mean(statX[i,]>= stat0X[i])
      pvalY[i] <- mean(statY[i,]>= stat0Y[i])
    }
    mstatrho = rowMeans(statrho)
    mstatX = rowMeans(statX)
    mstatY = rowMeans(statY)
    ci_uprho<-rowQuantiles(statrho, probs = c(.95))
    ci_upX<-rowQuantiles(statX, probs = c(.95))
    ci_upY<-rowQuantiles(statY, probs = c(.95))
    ci_lowrho<-rowQuantiles(statrho, probs = c(.05))
    ci_lowX<-rowQuantiles(statX, probs = c(.05))
    ci_lowY<-rowQuantiles(statY, probs = c(.05))
    tabrho <- cbind(data.frame(stat0rho), ci_lowrho, mstatrho, ci_uprho, data.frame(pvalrho))
    tabX <- cbind(data.frame(stat0X), ci_lowX, mstatX, ci_upX, data.frame(pvalX))
    tabY <- cbind(data.frame(stat0Y), ci_lowY, mstatY, ci_upY, data.frame(pvalY))
    tabrhoXY <- list("rho"=tabrho, "X"=tabX, "Y"=tabY)
    return(tabrhoXY)
}

radarchartvar<-function (df, axistype = 0, seg = 4, pty = 16, pcol = 1:8, plty = 1:6, 
                         plwd = 1, pdensity = NULL, pangle = 45, pfcol = NA, cglty = 3, 
                         cglwd = 1, cglcol = "navy", axislabcol = "blue", title = "", 
                         maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlabels = NULL, 
                         vlcex = NULL, caxislabels = NULL, calcex = NULL, paxislabels = NULL, 
                         seglty = 1, seglwd = 3, segcol = "black", palcex = NULL, ...) 
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
    polygon(xx * ((seg/2) + CGap)/(seg + CGap), yy * ((seg/2) + CGap)/(seg + CGap), lty = 1, lwd = 4.5, border = "black")
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
  }

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
               axis.text.x = element_text(size = 25, color="black"), 
               axis.title.x = element_text(size = 30,margin=margin(15,0,0,0)),
               axis.text.y = element_text(size = 25, color="black",angle=0,hjust=0.5,margin=margin(0,10,0,0)), 
               axis.title.y = element_blank(),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               plot.background = element_rect(fill = "transparent", colour = NA),
               panel.background = element_rect(fill = "white", colour = NA),
               legend.position="none",
               axis.ticks = element_line(colour = "black", size=2),
               axis.ticks.length = unit(.4, "cm"))

theme3<- theme(panel.border = element_rect(size = 2, colour = 'NA', fill='NA'), 
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

# -- Read behavioral data
beh_raw <- read.table("../data_symptoms/BSNIP_PSD_N436_Behavior.dat", header=T)
beh <- beh_raw[,-1]
Group=beh_raw$Group

behpc_raw <- read.table("../fig1_symptom_pca_analyses/BSNIP_PSD436_AllPCScores.txt")
behpc <- behpc_raw[,2:6]

# -- Read neural data
gbcPSD_180pS <- t(read.table("BSNIP_PSD_N436_GBC_180P_SymmetrizedCortex.txt"))

##############################################################################################################################
# -- 180 neural vs. 36 symptom measure CCA 
ccNP180B36 <- cc(beh, gbcPSD_180pS)
ccNP180B36.permVartest<-ccaVarperm(beh, gbcPSD_180pS, nboot=10000)
p.adjust(ccNP180B36.permVartest$rho$pval, method="fdr")

signs<-rep(1,36)

# -- Organize for plotting
ccNP180B36behfacload <- as.data.frame(cbind("PC"=1:36, t(t(ccNP180B36$scores$corr.X.xscores)*signs)))
ccNP180B36behcoef <- as.data.frame(cbind("PC"=1:36, t(t(ccNP180B36$xcoef)*signs)))
ccNP180B36behcoeforig <- as.data.frame(cbind("PC"=1:36, ccNP180B36$xcoef))
names(ccNP180B36behcoef) <- c("PC", "CV1", "CV2", "CV3", "CV4", "CV5")
bwr <- colorRampPalette(c("blue","white","red"))(101)

ccNP180B36neurfacload <- as.data.frame(cbind("Parcel"=1:180, t(t(ccNP180B36$scores$corr.Y.yscores)*signs)))
ccNP180B36neurcoef <- as.data.frame(cbind("Parcel"=1:180, t(t(ccNP180B36$ycoef)*signs)))
names(ccNP180B36neurfacload) <- c("Parcel", "CV1", "CV2", "CV3", "CV4", "CV5")
names(ccNP180B36neurcoef) <- c("Parcel", "CV1", "CV2", "CV3", "CV4", "CV5")

ccNP180B36xscoresorig <- as.data.frame(ccNP180B36$scores$xscores)
ccNP180B36yscoresorig <- as.data.frame(ccNP180B36$scores$yscores)
ccNP180B36xscores <- as.data.frame(t(t(ccNP180B36$scores$xscores)*signs))
ccNP180B36yscores <- as.data.frame(t(t(ccNP180B36$scores$yscores)*signs))

# -- Plot the screeplot showing the canonical correlations for 180 vs. 36 CCA
pdf("CCA_ParcellatedN180B36_RhoScreeplot.pdf", width=10,height=6)
plotdatPN180B36 <- cbind("Mode"=1:dim(ccNP180B36.permVartest$rho)[1],ccNP180B36.permVartest$rho)
xcoordi <- c(rbind((1:dim(ccNP180B36.permVartest$rho)[1])-0.5,(1:dim(ccNP180B36.permVartest$rho)[1])-0.5), (dim(ccNP180B36.permVartest$rho)[1])+0.5)
xcoord <- head(c(xcoordi, rev(xcoordi))[-1], -1)
ycoordi <- c(rbind(ccNP180B36.permVartest$rho$ci_uprho, ccNP180B36.permVartest$rho$ci_uprho))
ycoord <- c(ycoordi, c(rbind(rev(ccNP180B36.permVartest$rho$ci_lowrho),rev(ccNP180B36.permVartest$rho$ci_lowrho))))
xycoords <- data.frame(cbind(xcoord,ycoord))
ggplot(plotdatPN180B36, aes(x=Mode, y=stat0rho)) +
  geom_polygon(data=xycoords, aes(x=xcoord, y=ycoord), alpha=0.3) +
  geom_line(aes(y=plotdatPN180B36$mstatrho), col="grey30", lwd=1,lty=2) +
  geom_line(lwd=1.2, col="red") +
  geom_point(bg="red",aes(size=stat0rho),col="black",pch=21) +
  scale_x_continuous(expand = c(0, 0),limits=c(0,37), breaks=seq(0,36,5),labels=seq(0,36,5)) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.88), breaks=seq(0, 0.8, by=0.1),labels=format(seq(0,0.8,0.1), nsmall = 1)) +
  xlab("Canonical Mode") +
  #ylab("Canonical Correlation") +
  theme2
dev.off()

# -- Plot the screeplot showing behavioral variance explained for 180 vs. 36 CCA
pdf("CCA_ParcellatedN180B36_BehaviorVarianceScreeplot.pdf", width=10,height=6)
plotdatPN180B36 <- cbind("Mode"=1:dim(ccNP180B36.permVartest$X)[1], ccNP180B36.permVartest$X)
xcoordi <- c(rbind((1:dim(ccNP180B36.permVartest$X)[1])-0.5,(1:dim(ccNP180B36.permVartest$X)[1])-0.5), (dim(ccNP180B36.permVartest$X)[1])+0.5)
xcoord <- head(c(xcoordi, rev(xcoordi))[-1], -1)
ycoordi <- c(rbind(ccNP180B36.permVartest$X$ci_upX, ccNP180B36.permVartest$X$ci_upX))
ycoord <- c(ycoordi, c(rbind(rev(ccNP180B36.permVartest$X$ci_lowX),rev(ccNP180B36.permVartest$X$ci_lowX))))
xycoords <- data.frame(cbind(xcoord,ycoord))
ggplot(plotdatPN180B36, aes(x=Mode, y=stat0X)) +
  geom_polygon(data=xycoords, aes(x=xcoord, y=ycoord), alpha=0.3) +
  geom_line(aes(y=plotdatPN180B36$mstatX), lwd=1, lty=2, col="grey30") +
  geom_line(lwd=1.2, col="forest green") +
  geom_point(bg="forest green",aes(size=stat0X),col="black",pch=21) +
  scale_x_continuous(expand = c(0, 0),limits=c(0,37), breaks=seq(0,36,5),labels=seq(0,36,5)) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.062), breaks=seq(0, 0.06, by=0.01),labels=format(seq(0,6,1), nsmall = 1)) +
  xlab("Canonical Mode") +
  theme2
dev.off()

# -- Plot the screeplot showing neural variance explained for 180 vs. 36 CCA
plotdatPN180B36 <- cbind("Mode"=1:dim(ccNP180B36.permVartest$X)[1], ccNP180B36.permVartest$X)
pdf("CCA_ParcellatedN180B36_NeuralVarianceScreeplot.pdf", width=10,height=6)
plotdatPN180B36 <- cbind("Mode"=1:dim(ccNP180B36.permVartest$Y)[1], ccNP180B36.permVartest$Y)
xcoordi <- c(rbind((1:dim(ccNP180B36.permVartest$Y)[1])-0.5,(1:dim(ccNP180B36.permVartest$Y)[1])-0.5), (dim(ccNP180B36.permVartest$Y)[1])+0.5)
xcoord <- head(c(xcoordi, rev(xcoordi))[-1], -1)
ycoordi <- c(rbind(ccNP180B36.permVartest$Y$ci_upY, ccNP180B36.permVartest$Y$ci_upY))
ycoord <- c(ycoordi, c(rbind(rev(ccNP180B36.permVartest$Y$ci_lowY),rev(ccNP180B36.permVartest$Y$ci_lowY))))
xycoords <- data.frame(cbind(xcoord,ycoord))
ggplot(plotdatPN180B36, aes(x=Mode, y=stat0Y)) +
  geom_polygon(data=xycoords, aes(x=xcoord, y=ycoord), alpha=0.3) +
  geom_line(aes(y=plotdatPN180B36$mstatY), lwd=1, lty=2, col="grey30") +
  geom_line(lwd=1.2, col="forest green") +
  geom_point(bg="forest green",aes(size=stat0Y),col="black",pch=21) +
  scale_x_continuous(expand = c(0, 0),limits=c(0,37), breaks=seq(0,36,5),labels=seq(0,36,5)) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.01), breaks=seq(0, 0.008, by=0.002),labels=format(seq(0,0.8,0.2), nsmall = 1)) +
  xlab("Canonical Mode") +
  theme2
dev.off()

# -- Pie chart of total behavioral var explained in 180 vs. 36 CCA
pdf("CCA_ParcellatedN180B36_BehaviorPiechart.pdf", width=6,height=6)
pve_prop9 <- data.frame("Type"=c("Explained","Not Explained"),"Total Variance"=c(sum(colMeans(ccNP180B36$scores$corr.X.yscores^2)), 1-sum(colMeans(ccNP180B36$scores$corr.X.yscores^2))))
ggplot(data=pve_prop9,aes(x=factor(1),y=Total.Variance,fill=factor(Type))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta="y", direction=1) +
  scale_fill_manual(values = c("Explained"="#2EBE31","Not Explained"="grey35")) +
  geom_text(aes(y = c(0.77,0.28),label = sprintf("%.2f%%", 100*pve_prop9$Total.Variance)), size=11, col=c("grey10","white")) +
  theme2 +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line=element_blank(),
        legend.position="none",
        panel.background = element_rect(fill = "transparent", colour = NA))
dev.off()

# -- Pie chart of total neural var explained in 180 vs. 36 CCA
pdf("CCA_ParcellatedN180B36_NeuralPiechart.pdf", width=6,height=6)
pve_prop9n <- data.frame("Type"=c("Explained","Not Explained"),"Total Variance"=c(sum(colMeans(ccNP180B36$scores$corr.Y.xscores^2)), 1-sum(colMeans(ccNP180B36$scores$corr.Y.xscores^2))))
ggplot(data=pve_prop9n,aes(x=factor(1),y=Total.Variance,fill=factor(Type))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta="y",direction=-1, start=1*pi/3) +
  scale_fill_manual(values = c("Explained"="#2EBE31","Not Explained"="grey35")) +
  geom_label(aes(y = c(0.97,0.45),label = sprintf("%.2f%%", 100*pve_prop9n$Total.Variance)), size=11, col=c("grey10","white"), fill=c("#2EBE31","grey35"), label.size=NA, alpha=c(1,1)) +
  theme2 +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line=element_blank(),
        legend.position="none",
        panel.background = element_rect(fill = "transparent", colour = NA))
dev.off()

##############################################################################################################################
# -- 180 neural vs. 5 PC CCA 
ccNP180B5 <- cc(behpc, gbcPSD_180pS) # Executing CCA for 180 symmetrized parcels
ccNP180B5.permVartest<-ccaVarperm(behpc, gbcPSD_180pS, nboot=10000)
p.adjust(ccNP180B5.permVartest$rho$pval, method="fdr")

#-- Flip signs for CV2,3!
signs<-c(1,-1,-1,1,1)

# -- Organize for plotting
ccNP180B5behfacload <- as.data.frame(cbind("PC"=1:5, t(t(ccNP180B5$scores$corr.X.xscores)*signs)))
ccNP180B5behcoef <- as.data.frame(cbind("PC"=1:5, t(t(ccNP180B5$xcoef)*signs)))
ccNP180B5behcoeforig <- as.data.frame(cbind("PC"=1:5, ccNP180B5$xcoef))
names(ccNP180B5behcoef) <- c("PC", "CV1", "CV2", "CV3", "CV4", "CV5")
names(ccNP180B5behfacload) <- c("PC", "CV1", "CV2", "CV3", "CV4", "CV5")
bwr <- colorRampPalette(c("blue","white","red"))(101)

ccNP180B5neurfacload <- as.data.frame(cbind("Parcel"=1:180, t(t(ccNP180B5$scores$corr.Y.yscores)*signs)))
ccNP180B5neurcoef <- as.data.frame(cbind("Parcel"=1:180, t(t(ccNP180B5$ycoef)*signs)))
names(ccNP180B5neurfacload) <- c("Parcel", "CV1", "CV2", "CV3", "CV4", "CV5")
names(ccNP180B5neurcoef) <- c("Parcel", "CV1", "CV2", "CV3", "CV4", "CV5")

ccNP180B5xscoresorig <- as.data.frame(ccNP180B5$scores$xscores)
ccNP180B5yscoresorig <- as.data.frame(ccNP180B5$scores$yscores)
ccNP180B5xscores <- as.data.frame(t(t(ccNP180B5$scores$xscores)*signs))
ccNP180B5yscores <- as.data.frame(t(t(ccNP180B5$scores$yscores)*signs))

# -- Plot loading of symptom PCs on each CV
for (cvno in 1:5){
  pdf(paste0("CCA_ParcellatedN180B5_CV",cvno,"_BehaviorFactorLoadings.pdf"), width=10,height=6)
  q<-ggplot(data=ccNP180B5behfacload, aes(x=PC, y=get(paste0("CV",cvno)),fill=get(paste0("CV",cvno)),order=PC)) +
    geom_col(col="black", lwd=2) +
    scale_fill_gradient2(low="blue",mid="white",high="red") +
    scale_y_continuous(expand = c(0, 0),limits=c(-1.1,1.1), breaks=seq(-0.5, 0.5, by=0.5),labels=format(seq(-0.5, 0.5, by=0.5), nsmall = 1)) +
    ylab("Loadings") +
    theme1 + theme(axis.text.x = element_text(size = 40, color="black"), 
                   axis.title.x = element_text(size = 40,margin=margin(15,0,0,0)),
                   axis.text.y = element_text(size = 40, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
                   axis.title.y = element_text(size = 40,margin=margin(0,15,0,0)))
  print(q)
  dev.off()
}

# -- Plot the screeplot showing the canonical correlations for 180 vs.5 CCA
pdf("CCA_ParcellatedN180B5_RhoScreeplot.pdf", width=10,height=6)
plotdatPN180B5 <- cbind("Mode"=1:dim(ccNP180B5.permVartest$rho)[1],ccNP180B5.permVartest$rho)
xcoordi <- c(rbind((1:dim(ccNP180B5.permVartest$rho)[1])-0.5,(1:dim(ccNP180B5.permVartest$rho)[1])-0.5), (dim(ccNP180B5.permVartest$rho)[1])+0.5)
xcoord <- head(c(xcoordi, rev(xcoordi))[-1], -1)
ycoordi <- c(rbind(ccNP180B5.permVartest$rho$ci_uprho, ccNP180B5.permVartest$rho$ci_uprho))
ycoord <- c(ycoordi, c(rbind(rev(ccNP180B5.permVartest$rho$ci_lowrho),rev(ccNP180B5.permVartest$rho$ci_lowrho))))
xycoords <- data.frame(cbind(xcoord,ycoord))
ggplot(plotdatPN180B5, aes(x=Mode, y=stat0rho)) +
  geom_polygon(data=xycoords, aes(x=xcoord, y=ycoord), alpha=0.3) +
  geom_line(aes(y=plotdatPN180B5$mstatrho), col="grey30", lwd=1,lty=2) +
  geom_line(lwd=1.2, col="red") +
  geom_point(bg="red",aes(size=stat0rho),col="black",pch=21) +
  scale_x_continuous(expand = c(0, 0),limits=c(0,6), breaks=c(1:5),labels=c(1:5)) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.88), breaks=seq(0, 0.8, by=0.1),labels=format(seq(0,0.8,0.1), nsmall = 1)) +
  xlab("Canonical Mode") +
  #ylab("Canonical Correlation") +
  theme2
dev.off()

# -- Plot the screeplot showing behavioral variance explained for 180 vs. 5 CCA, scaled by total symptom variance explained by the symptom PCs
proptotvar = 0.5093
pdf("CCA_ParcellatedN180B5_BehaviorVarianceScreeplot.pdf", width=10,height=6)
plotdatPN180B5 <- cbind("Mode"=1:dim(ccNP180B5.permVartest$X)[1], proptotvar*ccNP180B5.permVartest$X)
xcoordi <- c(rbind((1:dim(ccNP180B5.permVartest$X)[1])-0.5,(1:dim(ccNP180B5.permVartest$X)[1])-0.5), (dim(ccNP180B5.permVartest$X)[1])+0.5)
xcoord <- head(c(xcoordi, rev(xcoordi))[-1], -1)
ycoordi <- c(rbind(proptotvar*ccNP180B5.permVartest$X$ci_upX, proptotvar*ccNP180B5.permVartest$X$ci_upX))
ycoord <- c(ycoordi, c(rbind(rev(proptotvar*ccNP180B5.permVartest$X$ci_lowX),rev(proptotvar*ccNP180B5.permVartest$X$ci_lowX))))
xycoords <- data.frame(cbind(xcoord,ycoord))
ggplot(plotdatPN180B5, aes(x=Mode, y=stat0X)) +
  geom_polygon(data=xycoords, aes(x=xcoord, y=ycoord), alpha=0.3) +
  geom_line(aes(y=plotdatPN180B5$mstatX), lwd=1, lty=2, col="grey30") +
  geom_line(lwd=1.2, col="blue") +
  geom_point(bg="blue",aes(size=stat0X),col="black",pch=21) +
  scale_x_continuous(expand = c(0, 0),limits=c(0,6), breaks=c(1:5),labels=c(1:5)) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.062), breaks=seq(0, 0.06, by=0.01),labels=format(seq(0,6,1), nsmall = 1)) +
  xlab("Canonical Mode") +
  theme2
dev.off()

# -- Plot the screeplot showing neural variance explained for 180 vs. 5 CCA
pdf("CCA_ParcellatedN180B5_NeuralVarianceScreeplot.pdf", width=10,height=6)
plotdatPN180B5 <- cbind("Mode"=1:dim(ccNP180B5.permVartest$Y)[1], ccNP180B5.permVartest$Y)
xcoordi <- c(rbind((1:dim(ccNP180B5.permVartest$Y)[1])-0.5,(1:dim(ccNP180B5.permVartest$Y)[1])-0.5), (dim(ccNP180B5.permVartest$Y)[1])+0.5)
xcoord <- head(c(xcoordi, rev(xcoordi))[-1], -1)
ycoordi <- c(rbind(ccNP180B5.permVartest$Y$ci_upY, ccNP180B5.permVartest$Y$ci_upY))
ycoord <- c(ycoordi, c(rbind(rev(ccNP180B5.permVartest$Y$ci_lowY),rev(ccNP180B5.permVartest$Y$ci_lowY))))
xycoords <- data.frame(cbind(xcoord,ycoord))
ggplot(plotdatPN180B5, aes(x=Mode, y=stat0Y)) +
  geom_polygon(data=xycoords, aes(x=xcoord, y=ycoord), alpha=0.3) +
  geom_line(aes(y=plotdatPN180B5$mstatY), lwd=1, lty=2, col="grey30") +
  geom_line(lwd=1.2, col="blue") +
  geom_point(bg="blue",aes(size=stat0Y),col="black",pch=21) +
  scale_x_continuous(expand = c(0, 0),limits=c(0,6), breaks=c(1:5),labels=c(1:5)) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.011), breaks=seq(0, 0.01, by=0.002),labels=format(seq(0,1,0.2), nsmall = 1)) +
  xlab("Canonical Mode") +
  theme2
dev.off()

# -- Pie chart of total behavioral var explained in 180 vs. 5 CCA, scaled by total symptom variance explained by the symptom PCs
pdf("CCA_ParcellatedN180B5_BehaviorPiechart.pdf", width=10,height=6)
proptotvar = 0.5093
pve_prop10 <- data.frame("Type"=c("Explained","Not Explained"),"Total Variance"=c(proptotvar*sum(colMeans(ccNP180B5$scores$corr.X.yscores^2)), 1-proptotvar*sum(colMeans(ccNP180B5$scores$corr.X.yscores^2))))
ggplot(data=pve_prop10,aes(x=factor(1),y=Total.Variance,fill=factor(Type))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta="y", direction=-1, start=5*pi/3) +
  scale_fill_manual(values = c("Explained"="blue","Not Explained"="grey35")) +
  geom_text(aes(y = c(0.9,0.4),label = sprintf("%.2f%%", 100*pve_prop10$Total.Variance)),nudge_x=c(0.05,0), size=11, col=c("white","white"), alpha=c(1,1)) +
  theme2 +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line=element_blank(),
        legend.position="none",
        panel.background = element_rect(fill = "transparent", colour = NA))
dev.off()

# -- Pie chart of total neural var explained in 180 vs. 5 CCA
pdf("CCA_ParcellatedN180B5_NeuralPiechart.pdf", width=10,height=6)
pve_prop10n <- data.frame("Type"=c("Explained","Not Explained"),"Total Variance"=c(sum(colMeans(ccNP180B5$scores$corr.Y.xscores^2)), 1-sum(colMeans(ccNP180B5$scores$corr.Y.xscores^2))))
ggplot(data=pve_prop10n,aes(x=factor(1),y=Total.Variance,fill=factor(Type))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta="y",direction=-1, start=1*pi/3) +
  scale_fill_manual(values = c("Explained"="blue","Not Explained"="grey35")) +
  geom_label(aes(y = c(0.99,0.5),label = sprintf("%.2f%%", 100*pve_prop10n$Total.Variance)), size=11, col=c("white","white"), fill=c("blue","grey35"), label.size=NA, alpha=c(1,1)) +
  theme2 +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line=element_blank(),
        legend.position="none",
        panel.background = element_rect(fill = "transparent", colour = NA))
dev.off()

# -- Plot the scatterplot showing correlation between neural and behavioral scores for CV1
# -- Do this for 180 vs 36 CCA
scatdat9 <- data.frame(xscores=c(ccNP180B36$scores$xscores[,1]), yscores=c(ccNP180B36$scores$yscores[,1]), dx=factor(c("SADP","SCZP","BPP","SCZP","SADP","SCZP","SCZP","BPP","BPP","BPP","SCZP","BPP","SCZP","BPP","SCZP","SCZP","BPP","BPP","BPP","SADP","SCZP","SCZP","SCZP","BPP","BPP","SADP","BPP","SADP","BPP","SCZP","SADP","BPP","SADP","SCZP","BPP","BPP","SADP","SADP","BPP","SCZP","SCZP","SADP","BPP","SADP","SCZP","SADP","SADP","SADP","SADP","SCZP","SCZP","SADP","SADP","SCZP","BPP","SCZP","BPP","SCZP","SCZP","SADP","BPP","BPP","BPP","BPP","SADP","SCZP","BPP","BPP","SADP","BPP","BPP","BPP","BPP","SCZP","BPP","BPP","SADP","SCZP","SADP","SADP","SCZP","SADP","SCZP","SADP","BPP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","SADP","SCZP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","BPP","SADP","SADP","SCZP","SCZP","BPP","BPP","SCZP","SADP","SCZP","SCZP","SCZP","BPP","SCZP","SADP","SADP","BPP","BPP","SCZP","SCZP","BPP","SCZP","BPP","SADP","SCZP","SADP","BPP","SADP","SADP","BPP","BPP","SADP","SADP","BPP","BPP","SADP","BPP","BPP","SCZP","SCZP","SCZP","BPP","BPP","BPP","BPP","SADP","SADP","SADP","SCZP","SCZP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SCZP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SADP","SCZP","BPP","BPP","SCZP","SADP","SCZP","SADP","SCZP","BPP","SCZP","BPP","BPP","SADP","BPP","BPP","BPP","BPP","BPP","BPP","BPP","BPP","SADP","SCZP","SCZP","BPP","SADP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SADP","SCZP","SCZP","BPP","BPP","BPP","BPP","BPP","BPP","SADP","SCZP","SADP","SADP","SCZP","SCZP","BPP","SCZP","BPP","BPP","SADP","BPP","SADP","SCZP","SADP","SCZP","SADP","SCZP","BPP","SCZP","SCZP","BPP","SADP","SADP","SADP","BPP","SADP","SCZP","SADP","SADP","BPP","SADP","SCZP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SCZP","SADP","BPP","SCZP","BPP","BPP","SCZP","BPP","BPP","SCZP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SADP","SCZP","BPP","SCZP","SCZP","SADP","SADP","BPP","BPP","SADP","SADP","BPP","BPP","SCZP","SCZP","SADP","SADP","BPP","BPP","BPP","BPP","SCZP","SCZP","SCZP","SCZP","SCZP","BPP","SADP","SCZP","SCZP","SCZP","BPP","SADP","SCZP","SADP","SCZP","SCZP","BPP","BPP","SCZP","SCZP","SADP","BPP","SADP","BPP","SADP","SADP","SADP","SADP","SCZP","BPP","SADP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SADP","SADP","SADP","SCZP","BPP","SCZP","SCZP","SCZP","SADP","BPP","BPP","SCZP","SCZP","BPP","SADP","SADP","SCZP","SADP","SADP","SADP","SCZP","SADP","SCZP","SADP","SADP","SADP","SADP","SCZP","BPP","SCZP","BPP","SADP","SCZP","SCZP","BPP","SADP","SCZP","BPP","SCZP","BPP","SCZP","SADP","BPP","SCZP","SCZP","BPP","SCZP","BPP","BPP","BPP","SADP","SCZP","BPP","SCZP","BPP","SCZP","BPP","SCZP","SADP","SCZP","SADP","BPP","SADP","SCZP","SCZP","SCZP","SADP","BPP","SADP","BPP","SCZP","BPP","BPP","BPP","BPP","SADP","SADP","SADP","SCZP","SADP","SCZP","BPP","SADP","SCZP","BPP","SCZP","BPP","BPP","BPP","SCZP")))
pdf("CCA_ParcellatedN180B36_CV1_Scores_Scatterplot.pdf", width=6,height=6)
ggplot(scatdat9,aes(x=xscores,y=yscores, col=dx)) +
  geom_point(alpha=0.7, size=3) +
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_x_continuous(expand = c(0, 0),limits=c(-4.6,4.6), breaks=seq(-4,4,8),labels=seq(-4,4,8)) +
  scale_y_continuous(expand = c(0, 0),limits=c(-4.6,4.6), breaks=seq(-4,4,8),labels=seq(-4,4,8)) +
  xlab("Behavioral CV1 score") +
  ylab("Neural CV1 score") +
  theme1
dev.off()

# -- Do this for 180 vs 5 CCA
scatdat10 <- data.frame(xscores=c(ccNP180B5$scores$xscores[,1]), yscores=c(ccNP180B5$scores$yscores[,1]), dx=factor(c("SADP","SCZP","BPP","SCZP","SADP","SCZP","SCZP","BPP","BPP","BPP","SCZP","BPP","SCZP","BPP","SCZP","SCZP","BPP","BPP","BPP","SADP","SCZP","SCZP","SCZP","BPP","BPP","SADP","BPP","SADP","BPP","SCZP","SADP","BPP","SADP","SCZP","BPP","BPP","SADP","SADP","BPP","SCZP","SCZP","SADP","BPP","SADP","SCZP","SADP","SADP","SADP","SADP","SCZP","SCZP","SADP","SADP","SCZP","BPP","SCZP","BPP","SCZP","SCZP","SADP","BPP","BPP","BPP","BPP","SADP","SCZP","BPP","BPP","SADP","BPP","BPP","BPP","BPP","SCZP","BPP","BPP","SADP","SCZP","SADP","SADP","SCZP","SADP","SCZP","SADP","BPP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","SADP","SCZP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","BPP","SADP","SADP","SCZP","SCZP","BPP","BPP","SCZP","SADP","SCZP","SCZP","SCZP","BPP","SCZP","SADP","SADP","BPP","BPP","SCZP","SCZP","BPP","SCZP","BPP","SADP","SCZP","SADP","BPP","SADP","SADP","BPP","BPP","SADP","SADP","BPP","BPP","SADP","BPP","BPP","SCZP","SCZP","SCZP","BPP","BPP","BPP","BPP","SADP","SADP","SADP","SCZP","SCZP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SCZP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SADP","SCZP","BPP","BPP","SCZP","SADP","SCZP","SADP","SCZP","BPP","SCZP","BPP","BPP","SADP","BPP","BPP","BPP","BPP","BPP","BPP","BPP","BPP","SADP","SCZP","SCZP","BPP","SADP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SADP","SCZP","SCZP","BPP","BPP","BPP","BPP","BPP","BPP","SADP","SCZP","SADP","SADP","SCZP","SCZP","BPP","SCZP","BPP","BPP","SADP","BPP","SADP","SCZP","SADP","SCZP","SADP","SCZP","BPP","SCZP","SCZP","BPP","SADP","SADP","SADP","BPP","SADP","SCZP","SADP","SADP","BPP","SADP","SCZP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SCZP","SADP","BPP","SCZP","BPP","BPP","SCZP","BPP","BPP","SCZP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SADP","SCZP","BPP","SCZP","SCZP","SADP","SADP","BPP","BPP","SADP","SADP","BPP","BPP","SCZP","SCZP","SADP","SADP","BPP","BPP","BPP","BPP","SCZP","SCZP","SCZP","SCZP","SCZP","BPP","SADP","SCZP","SCZP","SCZP","BPP","SADP","SCZP","SADP","SCZP","SCZP","BPP","BPP","SCZP","SCZP","SADP","BPP","SADP","BPP","SADP","SADP","SADP","SADP","SCZP","BPP","SADP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SADP","SADP","SADP","SCZP","BPP","SCZP","SCZP","SCZP","SADP","BPP","BPP","SCZP","SCZP","BPP","SADP","SADP","SCZP","SADP","SADP","SADP","SCZP","SADP","SCZP","SADP","SADP","SADP","SADP","SCZP","BPP","SCZP","BPP","SADP","SCZP","SCZP","BPP","SADP","SCZP","BPP","SCZP","BPP","SCZP","SADP","BPP","SCZP","SCZP","BPP","SCZP","BPP","BPP","BPP","SADP","SCZP","BPP","SCZP","BPP","SCZP","BPP","SCZP","SADP","SCZP","SADP","BPP","SADP","SCZP","SCZP","SCZP","SADP","BPP","SADP","BPP","SCZP","BPP","BPP","BPP","BPP","SADP","SADP","SADP","SCZP","SADP","SCZP","BPP","SADP","SCZP","BPP","SCZP","BPP","BPP","BPP","SCZP")))
pdf("CCA_ParcellatedN180B5_CV1_Scores_Scatterplot.pdf", width=6,height=6)
ggplot(scatdat10,aes(x=xscores,y=yscores, col=dx)) +
  geom_point(alpha=0.7, size=3) +
  scale_colour_manual(values = c("BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_x_continuous(expand = c(0, 0),limits=c(-4.6,3.6), breaks=seq(-4,3,7),labels=seq(-4,3,7)) +
  scale_y_continuous(expand = c(0, 0),limits=c(-4.6,3.6), breaks=seq(-4,3,7),labels=seq(-4,3,7)) +
  xlab("Behavioral CV1 score") +
  ylab("Neural CV1 score") +
  theme1
dev.off()

# -- Plot the ridgeline plots showing distributions of scores in each diagnostic category
### ! Flipped signs to follow convention in figure!
gbc180beh5_ccascoresPSDG <- cbind(Group, as.data.frame(t(t(ccNP180B5$scores$yscores)*signs)))

# -- Read control data
behpc_CON_raw <- read.table("../fig1_symptom_pca_analyses/BSNIP_CON202_AllPCScores.txt")
behpc_CON <- behpc_CON_raw[,1:5]
gbcCON_180pS <- t(read.table("BSNIP_CON_N202_GBC_180P_SymmetrizedCortex.txt"))

# -- Project controls!
gbc180beh5_ccascoresCON <- as.data.frame(as.matrix(gbcCON_180pS)%*%t(t(ccNP180B5$ycoef)*signs))

# -- Normalize all groups to controls
gbc180beh5_ccascoresCON_Z <- scale(gbc180beh5_ccascoresCON)[,]
gbc180beh5_ccascoresCON_sd <- colSds(as.matrix(gbc180beh5_ccascoresCON))
gbc180beh5_ccascoresCON_m <- colMeans(as.matrix(gbc180beh5_ccascoresCON))

gbc180beh5_ccascoresPSDG_Z <- cbind(gbc180beh5_ccascoresPSDG[1],sweep(sweep(gbc180beh5_ccascoresPSDG[-1],2,gbc180beh5_ccascoresCON_m),2,gbc180beh5_ccascoresCON_sd,"/"))
gbc180beh5_ccascoresALL_Z <- rbind(gbc180beh5_ccascoresPSDG_Z, cbind(Group="PSD", gbc180beh5_ccascoresPSDG_Z[-1]), cbind(Group="CON", as.data.frame(gbc180beh5_ccascoresCON_Z)))
gbc180beh5_ccascoresALL_Z$Group <- factor(gbc180beh5_ccascoresALL_Z$Group, levels = c("SCZP","SADP","BPP","PSD","CON"))
gbc180beh5_ccascoresALL_Zmeans <- aggregate(gbc180beh5_ccascoresALL_Z,by=list(gbc180beh5_ccascoresALL_Z$Group), FUN="mean")

# -- Plot CV score distributions for each of the behavioral CVs in the 180 vs 5 CCA
pdf("CCA_ParcellatedN180B5_CV1_NormalizedScoresDistributions.pdf", width=8,height=6)
ggplot(gbc180beh5_ccascoresALL_Z, aes(x=V1, y=Group,group=Group)) + 
  geom_density_ridges(aes(fill=Group), lwd=2, rel_min_height = 0.0000005, scale=1.5) +
  scale_fill_manual(values = c("CON"="white","PSD"="black","BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-5,5), breaks=seq(-5,5,1),labels=c(-5,"","","","",0,"","","","",5)) +
  ylab("CV1") +
  xlab("") +
  theme3 +
  geom_segment(aes(x = gbc180beh5_ccascoresALL_Zmeans$V1[1], y=1.0,xend = gbc180beh5_ccascoresALL_Zmeans$V1[1], yend = 2),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = gbc180beh5_ccascoresALL_Zmeans$V1[2], y=2.0,xend = gbc180beh5_ccascoresALL_Zmeans$V1[2], yend = 3),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = gbc180beh5_ccascoresALL_Zmeans$V1[3], y=3.0,xend = gbc180beh5_ccascoresALL_Zmeans$V1[3], yend = 4),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = gbc180beh5_ccascoresALL_Zmeans$V1[4], y=4.0,xend = gbc180beh5_ccascoresALL_Zmeans$V1[4], yend = 5),size=2, linetype = 1,color = "white") +
  geom_vline(xintercept=gbc180beh5_ccascoresALL_Zmeans$V1[5], lty=2,col="grey",lwd=3)
dev.off()

pdf("CCA_ParcellatedN180B5_CV2_NormalizedScoresDistributions.pdf", width=8,height=6)
ggplot(gbc180beh5_ccascoresALL_Z, aes(x=V2, y=Group,group=Group)) + 
  geom_density_ridges(aes(fill=Group), lwd=2, rel_min_height = 0.001, scale=1) +
  scale_fill_manual(values = c("CON"="white","PSD"="black","BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-5,5), breaks=seq(-5,5,1),labels=c(-5,"","","","",0,"","","","",5)) +
  ylab("CV2") +
  xlab("") +
  theme3 +
  geom_segment(aes(x = gbc180beh5_ccascoresALL_Zmeans$V2[1], y=1.0,xend = gbc180beh5_ccascoresALL_Zmeans$V2[1], yend = 2),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = gbc180beh5_ccascoresALL_Zmeans$V2[2], y=2.0,xend = gbc180beh5_ccascoresALL_Zmeans$V2[2], yend = 3),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = gbc180beh5_ccascoresALL_Zmeans$V2[3], y=3.0,xend = gbc180beh5_ccascoresALL_Zmeans$V2[3], yend = 4),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = gbc180beh5_ccascoresALL_Zmeans$V2[4], y=4.0,xend = gbc180beh5_ccascoresALL_Zmeans$V2[4], yend = 5),size=2, linetype = 1,color = "white") +
  geom_vline(xintercept=gbc180beh5_ccascoresALL_Zmeans$V2[5], lty=2,col="grey",lwd=3)
dev.off()

pdf("CCA_ParcellatedN180B5_CV3_NormalizedScoresDistributions.pdf", width=8,height=6)
ggplot(gbc180beh5_ccascoresALL_Z, aes(x=V3, y=Group,group=Group)) + 
  geom_density_ridges(aes(fill=Group), lwd=2, rel_min_height = 0.001, scale=1.5) +
  scale_fill_manual(values = c("CON"="white","PSD"="black","BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-5,5), breaks=seq(-5,5,1),labels=c(-5,"","","","",0,"","","","",5)) +
  ylab("CV3") +
  xlab("") +
  theme3 +
  geom_segment(aes(x = gbc180beh5_ccascoresALL_Zmeans$V3[1], y=1.0,xend = gbc180beh5_ccascoresALL_Zmeans$V3[1], yend = 2),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = gbc180beh5_ccascoresALL_Zmeans$V3[2], y=2.0,xend = gbc180beh5_ccascoresALL_Zmeans$V3[2], yend = 3),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = gbc180beh5_ccascoresALL_Zmeans$V3[3], y=3.0,xend = gbc180beh5_ccascoresALL_Zmeans$V3[3], yend = 4),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = gbc180beh5_ccascoresALL_Zmeans$V3[4], y=4.0,xend = gbc180beh5_ccascoresALL_Zmeans$V3[4], yend = 5),size=2, linetype = 1,color = "white") +
  geom_vline(xintercept=gbc180beh5_ccascoresALL_Zmeans$V3[5], lty=2,col="grey",lwd=3)
dev.off()

pdf("CCA_ParcellatedN180B5_CV4_NormalizedScoresDistributions.pdf", width=8,height=6)
ggplot(gbc180beh5_ccascoresALL_Z, aes(x=V4, y=Group,group=Group)) + 
  geom_density_ridges(aes(fill=Group), lwd=2, rel_min_height = 0.001, scale=1.5) +
  scale_fill_manual(values = c("CON"="white","PSD"="black","BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-5,5), breaks=seq(-5,5,1),labels=c(-5,"","","","",0,"","","","",5)) +
  ylab("CV4") +
  xlab("") +
  theme3 +
  geom_segment(aes(x =  gbc180beh5_ccascoresALL_Zmeans$V4[1], y=1.0,xend = gbc180beh5_ccascoresALL_Zmeans$V4[1], yend = 2),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  gbc180beh5_ccascoresALL_Zmeans$V4[2], y=2.0,xend = gbc180beh5_ccascoresALL_Zmeans$V4[2], yend = 3),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  gbc180beh5_ccascoresALL_Zmeans$V4[3], y=3.0,xend = gbc180beh5_ccascoresALL_Zmeans$V4[3], yend = 4),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  gbc180beh5_ccascoresALL_Zmeans$V4[4], y=4.0,xend = gbc180beh5_ccascoresALL_Zmeans$V4[4], yend = 5),size=2, linetype = 1,color = "white") +
  geom_vline(xintercept=gbc180beh5_ccascoresALL_Zmeans$V4[5], lty=2,col="grey",lwd=3)
dev.off()

pdf("CCA_ParcellatedN180B5_CV5_NormalizedScoresDistributions.pdf", width=8,height=6)
ggplot(gbc180beh5_ccascoresALL_Z, aes(x=V5, y=Group,group=Group)) + 
  geom_density_ridges(aes(fill=Group), lwd=2, rel_min_height = 0.00001, scale=1.5) +
  scale_fill_manual(values = c("CON"="white","PSD"="black","BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0), limits=c(-5,5), breaks=seq(-5,5,1),labels=c(-5,"","","","",0,"","","","",5)) +
  ylab("CV5") +
  xlab("") +
  theme3 +
  geom_segment(aes(x =  gbc180beh5_ccascoresALL_Zmeans$V5[1], y=1.0,xend = gbc180beh5_ccascoresALL_Zmeans$V5[1], yend = 2),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  gbc180beh5_ccascoresALL_Zmeans$V5[2], y=2.0,xend = gbc180beh5_ccascoresALL_Zmeans$V5[2], yend = 3),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  gbc180beh5_ccascoresALL_Zmeans$V5[3], y=3.0,xend = gbc180beh5_ccascoresALL_Zmeans$V5[3], yend = 4),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x =  gbc180beh5_ccascoresALL_Zmeans$V5[4], y=4.0,xend = gbc180beh5_ccascoresALL_Zmeans$V5[4], yend = 5),size=2, linetype = 1,color = "white") +
  geom_vline(xintercept=gbc180beh5_ccascoresALL_Zmeans$V5[5], lty=2,col="grey",lwd=3)
dev.off()

# -- Compute coef of the original 36 symptom measures on CVs in the 180 vs. 5 CCA
PCArot <- read.table("/Users/jielisaji/Dropbox/N-BRIDGE/BSNIP_Analyses/Data/PCArotations.txt")
CVbehload <- as.data.frame(as.matrix(PCArot[,1:5]) %*% as.matrix(ccNP180B5behcoef[2:6]))

# -- Radar plots
for (cvno in 1:5){
  assign(paste0("CV",cvno,"rad"),data.frame(rbind(rep(0.3,36), rep(-0.3,36),CVbehload[,cvno])))
}
assign("CV2rad",data.frame(rbind(rep(0.2,36), rep(-0.2,36),CVbehload[,2])))
assign("CV4rad",data.frame(rbind(rep(0.2,36), rep(-0.2,36),CVbehload[,4])))

Sxcols=c("#004E2E","#005F38","#006E45","#008953","#009E5D","#00B669",
         "#441978","#511D90","#6124AF","#6F26C3","#7C2CE2","#8E35F5","#A351FF",
         "#354252","#334862","#37506E","#406084","#41658F","#446C9C","#4572A6",
         "#5B013B","#6A0144","#7F0154","#950264","#AD016D","#C2017A","#CC0285","#D3038C","#DF0293","#ED039F",
         "#F702A6","#FF3DB9","#FF5ABE","#FF7FD4","#FFA1E4","#FEB5EB")

for (cvno in 1:5){
  pdf(paste0("CCA_ParcellatedN180B5_CV",cvno,"_RadarplotBehaviorLoadings.pdf"), width=6,height=6)
  radarchartvar(get(paste0("CV",cvno,"rad")), axistype=0, seglty = 1, seglwd = 2, segcol=c("grey"), cglty=1, cglwd=2, vlabels="", cglcol=c(Sxcols), pcol="grey30", plwd=8, seg=2)
  dev.off()
}