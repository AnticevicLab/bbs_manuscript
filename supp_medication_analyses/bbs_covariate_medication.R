rm(list = ls())
setwd("../supp_medication_analyses")
library(ggplot2)
library(rgl)
library(plot3D)

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

# -- PCA
beh_raw = read.table("../data_symptoms/BSNIP_PSD_N436_Behavior.dat", sep="\t", header=TRUE)
beh=beh_raw[1:436,2:37]
pcall = prcomp(beh, scale. = TRUE)

# -- Load CPZ equivalents
cpz <-c(NA,NA,0,506.9624712,416.1805903,0,625.8587106,326.8894131,188.9605615,0,NA,0,1164.562201,625.8587106,NA,337.233036,526.3539215,0,0,236.6680929,0,303.9468537,0,115.9356336,NA,611.019627,0,NA,543.0207344,524.9902048,131.2750346,NA,149.5598222,128.1309935,29.38811684,29.38811684,337.233036,337.233036,882.1682567,337.233036,730.4868471,747.6728529,326.8894131,0,NA,0,337.39336,154.3870315,337.233036,NA,1164.562201,350,625.8587106,NA,0,NA,0,0,NA,72.91565446,106.6911072,68.83067798,236.6680929,NA,508.6441123,NA,113.0825866,NA,0,NA,188.9605615,236.6680929,236.6680929,207.553076,NA,508.6441123,761.1287423,154.3870315,0,46.54193233,68.83067798,506.9624712,106.6911072,0,0,NA,337.233036,1088.018155,337.39336,NA,NA,236.6680929,240.6382223,NA,836.6714439,289.2892972,1168.578501,591.3726652,149.5598222,0,625.8587106,1464.306751,200,NA,377.1754848,170.0325108,13.02717258,NA,289.0518709,337.39336,NA,508.6441123,66.29691945,613.2193132,1245.157376,188.9605615,1117.469076,591.3726652,678.4877642,0,204.5261474,0,188.9605615,13.02717258,NA,658.2039344,337.233036,NA,0,2251.336674,209.2707475,NA,NA,440.2070472,1369.168417,154.3870315,108.665738,0,337.39336,131.2750346,NA,152.8654464,NA,873.9695353,NA,NA,326.8894131,NA,240.7103817,477.3784747,545.4085634,NA,337.233036,508.6441123,450.2915427,326.8894131,625.8587106,524.9902048,72.91565446,NA,907.7640161,1068.010939,154.3870315,240.7103817,239.4339733,NA,NA,NA,461.4016042,326.8894131,239.4339733,0,574.0614529,508.6441123,188.9605615,NA,0,543.0207344,245.4187056,506.9624712,68.83067798,0,106.6911072,0,337.39336,656.4222496,NA,NA,236.6680929,0,337.233036,NA,207.553076,239.4339733,188.9605615,244.9733843,239.4339733,NA,207.553076,523.4073542,1268.091214,NA,836.6714439,0,240.7103817,NA,508.6441123,105.035718,1887.2239,1148.474398,NA,100.0727801,907.7640161,625.8587106,560.5256714,244.9733843,989.0136125,1369.168417,907.7640161,NA,0,NA,892.8769434,524.9902048,239.4339733,188.9605615,22.70104123,264.876345,66.29691945,NA,0,NA,154.3870315,106.6911072,149.5598222,107.658405,NA,NA,1295.203981,72.91565446,240.7103817,326.8894131,836.6714439,NA,47.29896519,NA,NA,0,0,524.9902048,NA,NA,NA,154.3870315,NA,0,41.03105604,NA,NA,245.4187056,NA,823.8577162,68.83067798,326.8894131,1193.648331,NA,625.8587106,543.0207344,188.9605615,0,102.1893252,236.6680929,NA,564.6882581,118.6383515,0,207.553076,66.29691945,25,5.774688675,1131.38713,506.9624712,395.8487412,0,154.3870315,239.4339733,506.9624712,1112.290164,70.19857882,0,761.1287423,NA,508.6441123,NA,NA,129.5795702,264.876345,989.0136125,337.233036,0,1761.509691,NA,337.39336,NA,NA,416.1805903,0,154.3870315,331.0414174,508.6441123,NA,907.7640161,524.9902048,22.82626845,66.29691945,595.1887836,207.553076,NA,47.29896519,490.4887255,771.4510368,457.7391817,NA,895.69084,0,149.5598222,239.4339733,326.8894131,337.233036,149.5598222,NA,0,506.9624712,907.7640161,0,NA,416.1805903,337.39336,66.29691945,337.233036,239.4339733,188.9605615,0,72.91565446,16.26449941,NA,240.7103817,2692.394751,NA,239.4339733,207.553076,416.1805903,NA,NA,188.9605615,29.38811684,NA,326.8894131,0,761.1287423,NA,NA,491.6200675,154.3870315,922.7304971,NA,NA,NA,175.2067175,0,149.5598222,1369.168417,907.7640161,506.9624712,395.8487412,337.233036,326.8894131,873.3471273,NA,91.60744411,477.1754848,508.6441123,66.29691945,39.61734402,NA,240.7103817,188.9605615,591.3726652,337.233036,NA,NA,1865.786297,625.8587106,66.29691945,66.29691945,0,802.2672605,0,146.5467723,NA,0,154.3870315,337.233036,188.9605615,524.9902048,337.233036,NA,0,0,NA,508.6441123,0,244.9733843,NA,NA,NA,2228.084368,0,326.8894131,297.5943918,NA,0,543.0207344,456.6325955,78.47524527,907.7640161,438.4102284,0,NA,0,907.7640161,NA,0,0,0,395.8487412)
cpzdf <- data.frame("CPZ"=cpz, "PC1"=pcall$x[,1],"PC2"=pcall$x[,2],"PC3"=pcall$x[,3],"PC4"=pcall$x[,4],"PC5"=pcall$x[,5])

mcols[medstatus=='FALSE']<-'grey30'
mcols[medstatus=='TRUE']<-'#DCD200'
#mcols <- na.omit(mcols)

# -- Plot correlation scatterplots
for (pcno in 1:5){
  pdf(paste0("Symptom_PCA_MedicationEffect_Scatterplot_PC", pcno, ".pdf"), height=6,width=6)
	q <- ggplot(cpzdf, aes(x=CPZ, y=get(paste0("PC",pcno)), aes(color=mcols))) +
	geom_point(color=mcols) +
	scale_y_continuous(limits=c(-10,6)) +
	stat_smooth(method = "lm", size = 3, se=FALSE, col="red",fullrange = FALSE,na.rm=TRUE) +
	xlab("") +
	ylab(paste0("PC",pcno, " Score")) +
	theme1
	print(q)
  dev.off()
}

# -- Test significance
t.test(na.omit(cpzdf[cpzdf$CPZ!=0,2]), na.omit(cpzdf[cpzdf$CPZ==0,2]))
t.test(na.omit(cpzdf[cpzdf$CPZ!=0,3]), na.omit(cpzdf[cpzdf$CPZ==0,3]))
t.test(na.omit(cpzdf[cpzdf$CPZ!=0,4]), na.omit(cpzdf[cpzdf$CPZ==0,4]))
t.test(na.omit(cpzdf[cpzdf$CPZ!=0,5]), na.omit(cpzdf[cpzdf$CPZ==0,5]))
t.test(na.omit(cpzdf[cpzdf$CPZ!=0,6]), na.omit(cpzdf[cpzdf$CPZ==0,6]))

cpzm <- data.frame("Group"=c(rep("Med",279), rep("NoMed",59)),"PC1"=c(na.omit(cpzdf[cpzdf$CPZ!=0,2]), na.omit(cpzdf[cpzdf$CPZ==0,2])),"PC2"=c(na.omit(cpzdf[cpzdf$CPZ!=0,3]), na.omit(cpzdf[cpzdf$CPZ==0,3])),"PC3"=c(na.omit(cpzdf[cpzdf$CPZ!=0,4]), na.omit(cpzdf[cpzdf$CPZ==0,4])),"PC4"=c(na.omit(cpzdf[cpzdf$CPZ!=0,5]), na.omit(cpzdf[cpzdf$CPZ==0,5])),"PC5"=c(na.omit(cpzdf[cpzdf$CPZ!=0,6]), na.omit(cpzdf[cpzdf$CPZ==0,6])))

for (pcno in 1:5){
  pdf(paste0("Symptom_PCA_MedicationEffect_Barplot_PC", pcno, ".pdf"), height=6,width=6)
	assign(paste0("cpzmpc",pcno,"summ"), summarySE(cpzm, measurevar=paste0("PC",pcno), groupvars="Group"))
	q <- ggplot(data=get(paste0("cpzmpc",pcno,"summ")), aes(x=Group,y=get(paste0("PC",pcno)), fill=Group)) +
  	geom_col() +
    geom_errorbar(aes(x=Group, ymin=get(paste0("PC",pcno))-se, ymax=get(paste0("PC",pcno))+se), col="black", width=0.2, position=position_dodge(0.9), lwd=1.5) +  
  	scale_fill_manual(values = c("NoMed"="#DCD200","Med"="grey30")) +
  	scale_x_discrete(limits=c("NoMed","Med"),labels = c("","")) +
  	ylab(paste0("PC",pcno, " Score")) +
  	xlab("") +
 	theme1 + theme(axis.ticks.x = element_blank())
	print(q)
  dev.off()
}

PC123 <- pcall$x[,1:3]
medstatus <- (cpz == 0)
PC123med <- data.frame(PC123[1:436,], "Med"=medstatus, "MCols"=medstatus)
PC123med$MCols[PC123med$Med=='FALSE']<-'grey30'
PC123med$MCols[PC123med$Med=='TRUE']<-'#DCD200'
PC123med <- na.omit(PC123med)

colours=c("#004E2E","#005F38","#006E45","#008953","#009E5D","#00B669",
          "#441978","#511D90","#6124AF","#6F26C3","#7C2CE2","#8E35F5","#A351FF",
          "#354252","#334862","#37506E","#406084","#41658F","#446C9C","#4572A6",
          "#5B013B","#6A0144","#7F0154","#950264","#AD016D","#C2017A","#CC0285","#D3038C","#DF0293","#ED039F",
          "#F702A6","#FF3DB9","#FF5ABE","#FF7FD4","#FFA1E4","#FEB5EB")

z <- c(rbind(colours, colours))
coords <- NULL

for (i in 1:nrow(pcall$rotation)) {
  coords <- rbind(coords, rbind(c(0,0,0),16*pcall$rotation[i,1:3]))
}
coordsCog <- (35*colMeans(pcall$rotation[1:6,1:3]))
coordsPos <- (35*colMeans(pcall$rotation[7:13,1:3]))
coordsNeg <- (35*colMeans(pcall$rotation[14:20,1:3]))
coordsGen <- (35*colMeans(pcall$rotation[21:36,1:3]))

# -- Plot 3D bipot
par3d(cex=5, windowRect = c(20, 30, 800, 800))
rgl.viewpoint(  zoom = 1.2 )
plot3d(PC123med[1:3], aspect=TRUE,
       xlab=expression("PC1"),ylab=expression("PC2"),zlab=expression("PC3"), 
       type="p", axes=F, size=15, col=PC123med$MCols, alpha=0.7,
       xlim=c(-10,10), ylim=c(-10,10), zlim=c(-10,10))
arrow3d(p0=c(0,0,0),p1=coordsCog, col="green",  lwd=50,alpha=1, type="rotation",barblen=0.08,width=0.6,n=4)
arrow3d(p0=c(0,0,0),p1=coordsPos, col="#8A40F5",lwd=50,alpha=1, type="rotation",barblen=0.08,width=0.6,n=4)
arrow3d(p0=c(0,0,0),p1=coordsNeg, col="#00B1FF",lwd=50,alpha=1, type="rotation",barblen=0.08,width=0.6,n=4)
arrow3d(p0=c(0,0,0),p1=coordsGen, col="pink",   lwd=50,alpha=1, type="rotation",barblen=0.08,width=0.6,n=3)
axis3d(edge="x--", lwd=3, axes.len=1, labels=c("-10","",0,"","10"),cex=2)
axis3d(edge="y--", lwd=3, axes.len=1, labels=c("","",0,"","10"),cex=2)
axis3d(edge="z", lwd=3, axes.len=1, labels=c("","",0,"","10"),cex=2)
grid3d("x")
grid3d("y")
grid3d("z")
#view3d(userMatrix=rotationMatrix(2*pi * 1, 1, -1, -1))
#movie3d(spin3d(axis = c(0,1,0), rpm = 2), duration=30,  type = "gif", dir='.')

