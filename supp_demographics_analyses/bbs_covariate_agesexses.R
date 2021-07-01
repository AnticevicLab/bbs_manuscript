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

#-- Get sex, sex, ses (from available B-SNIP data)
beh_raw = read.table("Data/BSNIP_PROB_N436_Behavior.dat", sep="\t", header=TRUE)
beh=beh_raw[1:436,3:38]
panss <-beh[1:436,7:36]
bacs <- beh[1:436,1:6]
pcall = prcomp(beh, scale. = TRUE)
#pcall <- read.table("Data/PCAscoresPROB_origorder.txt", sep="\t", header=TRUE)

# -- Load Sex
sex <-c("F","M","F","M","M","F","M","M","F","M","F","F","M","F","M","M","M","F","F","M","F","F","F","F","F","F","F","F","F","M","M","M","F","F","F","F","M","F","M","F","M","F","M","F","F","M","F","F","F","M","M","F","M","M","M","M","M","M","M","F","F","F","F","M","F","M","M","F","F","F","F","F","F","F","M","F","M","F","F","F","M","F","M","M","F","M","M","M","M","M","M","M","M","M","M","F","M","M","F","F","M","M","M","F","F","F","F","F","F","M","F","M","M","M","F","M","M","M","F","M","M","M","F","F","F","F","F","M","M","M","M","M","F","M","F","F","F","F","M","M","F","M","F","F","M","F","F","F","M","M","M","M","M","M","M","M","M","M","F","F","M","M","F","M","F","F","F","F","F","M","M","M","M","F","F","F","F","F","F","M","F","F","M","F","F","M","M","F","M","F","M","F","F","M","M","M","F","M","F","M","M","M","M","M","M","M","F","F","F","M","F","F","M","F","M","F","M","F","F","F","F","F","F","M","F","M","M","F","F","M","F","F","F","M","F","M","M","F","M","F","F","M","F","M","F","M","M","F","M","M","F","F","M","M","F","F","M","M","M","F","M","M","F","F","M","F","M","F","M","M","F","M","F","F","F","F","M","F","F","F","F","M","M","F","F","M","M","M","M","F","M","M","F","F","F","M","M","M","F","F","M","M","F","M","F","M","M","F","F","M","M","M","F","F","F","M","F","M","F","F","M","M","M","M","F","F","M","M","F","M","F","F","M","F","M","M","F","F","F","M","F","M","F","M","F","M","M","M","M","F","F","F","M","M","F","F","F","M","M","F","F","F","M","F","M","M","F","F","M","F","F","F","M","F","F","M","F","F","M","M","M","M","F","F","M","F","M","M","F","M","M","F","M","F","F","F","F","F","F","F","F","M","F","F","M","M","M","F","F","M","M","F","F","F","F","M","M","M","F","M","M","M","F","M","M","M","M","M","F","F","F","M","F","F","F","F")
sexdf <- data.frame("Sex"=sex, "PC1"=pcall$x[,1],"PC2"=pcall$x[,2],"PC3"=pcall$x[,3],"PC4"=pcall$x[,4],"PC5"=pcall$x[,5])

sexdf_m <- aggregate(sexdf, by=list(sexdf$Sex), FUN="mean")[,-2]

for (pcno in 1:5){
  q <- ggplot(sexdf, aes(x=get(paste0("PC",pcno)),y=Sex, group=Sex)) + 
  geom_density_ridges(aes(fill=Sex), lwd=2, rel_min_height = 0.0005) +
  scale_fill_manual(values = c("F"="grey70", "M"="grey70")) +
  geom_segment(aes(x = eval(parse(text=paste0("sexdf_m$PC", pcno)))[2], y=1.0,xend = eval(parse(text=paste0("sexdf_m$PC", pcno)))[2], yend = 2.0),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = eval(parse(text=paste0("sexdf_m$PC", pcno)))[1], y=2.0,xend = eval(parse(text=paste0("sexdf_m$PC", pcno)))[1], yend = 3.0),size=2, linetype = 1,color = "white") +
  scale_y_discrete(limits=c("M","F"), expand=expansion(mult=c(2,2), add=c(-1.5,0))) +
	xlab(paste0("PC",pcno, " Score")) +
	ylab("Sex") +
	theme1 + 
      theme(axis.text.x = element_text(size = 35, color="black", margin=margin(5,0,0,0), vjust=1), 
        axis.title.x = element_text(size = 35, color="black"),
        axis.text.y = element_text(size = 35, color="black", angle=90, hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 35,margin=margin(0,5,0,0)),
        axis.line= element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size=2),
        axis.ticks.length = unit(.4, "cm"))       	
  print(q)
}

# -- Bar plots for mean effect of sex
sexPC1df_summ <- summarySE(sexdf, measurevar="PC1", groupvars=c("Sex"))
sexPC2df_summ <- summarySE(sexdf, measurevar="PC2", groupvars=c("Sex"))
sexPC3df_summ <- summarySE(sexdf, measurevar="PC3", groupvars=c("Sex"))
sexPC4df_summ <- summarySE(sexdf, measurevar="PC4", groupvars=c("Sex"))
sexPC5df_summ <- summarySE(sexdf, measurevar="PC5", groupvars=c("Sex"))

for (pcno in 1:5){
q <-  ggplot(data=get(paste0("sexPC",pcno,"df_summ")), aes(x=Sex,y=get(paste0("PC",pcno)), fill=Sex)) +
      geom_col() +
      scale_fill_manual(values = c("F"="grey70","M"="grey70")) +
      theme1 +
      geom_errorbar(aes(x=Sex, ymin=get(paste0("PC",pcno))-se, ymax=get(paste0("PC",pcno))+se), col="black", width=0.2, position=position_dodge(0.9), lwd=1.5) +  
      scale_y_continuous(expand = c(0.001, 0), limits=c(-0.48,.48), breaks=seq(-0.4,0.4,0.2),labels=format(c(-0.4,"",0.0,"",0.4), nsmall=1)) +
      xlab("Sex") +
      ylab(paste0("PC",pcno," Score")) +
      theme(axis.text.x = element_text(size = 35, color="black", margin=margin(5,0,0,0), vjust=1), 
        axis.title.x = element_text(size = 35, color="black"),
        axis.text.y = element_text(size = 35, color="black", angle=90, hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 35,margin=margin(0,5,0,0)),
        axis.line= element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size=2),
        axis.ticks.length = unit(.4, "cm"))       
print(q)
}

# -- Test significance
t.test(sexdf$PC1[sexdf$Sex=="F"], sexdf$PC1[sexdf$Sex=="M"])
t.test(sexdf$PC2[sexdf$Sex=="F"], sexdf$PC2[sexdf$Sex=="M"])
t.test(sexdf$PC3[sexdf$Sex=="F"], sexdf$PC3[sexdf$Sex=="M"])
t.test(sexdf$PC4[sexdf$Sex=="F"], sexdf$PC4[sexdf$Sex=="M"])
t.test(sexdf$PC5[sexdf$Sex=="F"], sexdf$PC5[sexdf$Sex=="M"])

# -- Load Age
age <-c(49,23,60,19,25,36,21,36,21,44,51,32,21,47,33,22,46,23,57,21,35,49,48,40,37,24,41,42,36,47,17,46,62,27,29,56,47,45,45,56,41,33,19,55,45,45,37,35,57,18,38,42,31,22,56,23,27,20,30,30,27,41,22,34,43,23,60,53,23,55,27,44,50,44,28,30,61,40,38,55,45,24,33,26,23,25,28,15,19,28,35,26,24,25,26,27,27,25,51,16,19,45,53,48,44,27,54,42,29,30,48,33,25,26,42,24,22,47,40,33,18,51,25,40,44,50,40,52,53,26,45,36,32,42,40,35,19,45,26,30,48,25,19,43,55,25,15,61,42,33,19,47,19,33,38,49,40,30,48,28,21,52,23,36,26,51,25,52,46,46,27,28,26,25,42,50,25,35,25,56,62,22,18,29,38,18,21,43,19,49,38,57,57,52,18,35,27,24,28,31,45,28,37,39,31,37,43,48,17,34,42,16,21,50,50,39,25,48,50,25,64,41,56,30,57,51,27,25,59,45,30,26,34,22,38,27,18,57,36,49,52,27,46,23,42,24,19,26,29,18,26,26,18,29,29,25,32,22,28,57,26,33,50,21,17,19,33,32,25,61,44,19,54,39,52,27,24,31,44,23,42,26,25,38,47,35,43,38,57,45,26,19,22,29,38,37,21,21,22,51,26,33,55,25,24,20,22,30,49,33,24,51,20,17,50,32,22,51,47,23,26,51,37,43,18,51,32,22,34,21,22,53,40,39,40,54,51,62,32,20,61,27,65,23,26,30,44,52,21,39,46,35,25,56,25,34,35,25,18,50,37,33,27,46,55,18,59,28,23,38,38,37,29,31,44,27,30,28,46,27,36,23,27,44,22,25,27,34,22,27,44,40,31,47,52,30,23,22,22,48,21,44,24,23,29,26,53,59,36,28,27,50,34,18,23,33,24,43,46,22,17,56,29,37,37,48,53,30,56,25,42,43,23,28,38,56)
agedf <- data.frame("Age"=age, "PC1"=pcall$x[,1],"PC2"=pcall$x[,2],"PC3"=pcall$x[,3],"PC4"=pcall$x[,4],"PC5"=pcall$x[,5])
Group=c("SADP","SCZP","BPP","SCZP","SADP","SCZP","SCZP","BPP","BPP","BPP","SCZP","BPP","SCZP","BPP","SCZP","SCZP","BPP","BPP","BPP","SADP","SCZP","SCZP","SCZP","BPP","BPP","SADP","BPP","SADP","BPP","SCZP","SADP","BPP","SADP","SCZP","BPP","BPP","SADP","SADP","BPP","SCZP","SCZP","SADP","BPP","SADP","SCZP","SADP","SADP","SADP","SADP","SCZP","SCZP","SADP","SADP","SCZP","BPP","SCZP","BPP","SCZP","SCZP","SADP","BPP","BPP","BPP","BPP","SADP","SCZP","BPP","BPP","SADP","BPP","BPP","BPP","BPP","SCZP","BPP","BPP","SADP","SCZP","SADP","SADP","SCZP","SADP","SCZP","SADP","BPP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","SADP","SCZP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","BPP","SADP","SADP","SCZP","SCZP","BPP","BPP","SCZP","SADP","SCZP","SCZP","SCZP","BPP","SCZP","SADP","SADP","BPP","BPP","SCZP","SCZP","BPP","SCZP","BPP","SADP","SCZP","SADP","BPP","SADP","SADP","BPP","BPP","SADP","SADP","BPP","BPP","SADP","BPP","BPP","SCZP","SCZP","SCZP","BPP","BPP","BPP","BPP","SADP","SADP","SADP","SCZP","SCZP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SCZP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SADP","SCZP","BPP","BPP","SCZP","SADP","SCZP","SADP","SCZP","BPP","SCZP","BPP","BPP","SADP","BPP","BPP","BPP","BPP","BPP","BPP","BPP","BPP","SADP","SCZP","SCZP","BPP","SADP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SADP","SCZP","SCZP","BPP","BPP","BPP","BPP","BPP","BPP","SADP","SCZP","SADP","SADP","SCZP","SCZP","BPP","SCZP","BPP","BPP","SADP","BPP","SADP","SCZP","SADP","SCZP","SADP","SCZP","BPP","SCZP","SCZP","BPP","SADP","SADP","SADP","BPP","SADP","SCZP","SADP","SADP","BPP","SADP","SCZP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SCZP","SADP","BPP","SCZP","BPP","BPP","SCZP","BPP","BPP","SCZP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SADP","SCZP","BPP","SCZP","SCZP","SADP","SADP","BPP","BPP","SADP","SADP","BPP","BPP","SCZP","SCZP","SADP","SADP","BPP","BPP","BPP","BPP","SCZP","SCZP","SCZP","SCZP","SCZP","BPP","SADP","SCZP","SCZP","SCZP","BPP","SADP","SCZP","SADP","SCZP","SCZP","BPP","BPP","SCZP","SCZP","SADP","BPP","SADP","BPP","SADP","SADP","SADP","SADP","SCZP","BPP","SADP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SADP","SADP","SADP","SCZP","BPP","SCZP","SCZP","SCZP","SADP","BPP","BPP","SCZP","SCZP","BPP","SADP","SADP","SCZP","SADP","SADP","SADP","SCZP","SADP","SCZP","SADP","SADP","SADP","SADP","SCZP","BPP","SCZP","BPP","SADP","SCZP","SCZP","BPP","SADP","SCZP","BPP","SCZP","BPP","SCZP","SADP","BPP","SCZP","SCZP","BPP","SCZP","BPP","BPP","BPP","SADP","SCZP","BPP","SCZP","BPP","SCZP","BPP","SCZP","SADP","SCZP","SADP","BPP","SADP","SCZP","SCZP","SCZP","SADP","BPP","SADP","BPP","SCZP","BPP","BPP","BPP","BPP","SADP","SADP","SADP","SCZP","SADP","SCZP","BPP","SADP","SCZP","BPP","SCZP","BPP","BPP","BPP","SCZP")

for (pcno in 1:5){
	q <- ggplot(agedf, aes(col=Group, x=Age, y=get(paste0("PC",pcno)))) +
	geom_point(aes(col=Group),size=3) +
  scale_colour_manual(values = c("CON"="white","PROB"="black","BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
	scale_y_continuous(limits=c(-10,6)) +
	stat_smooth(method = "lm", size = 3, se=FALSE, col="red",fullrange = FALSE,na.rm=TRUE) +
	xlab("Age") +
	ylab(paste0("PC",pcno, " Score")) +
  theme1+
  theme(axis.text.x = element_text(size = 35, color="black", margin=margin(5,0,0,0), vjust=1), 
        axis.title.x = element_text(size = 35, color="black"),
        axis.text.y = element_text(size = 35, color="black", angle=90, hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 35,margin=margin(0,5,0,0)),
        axis.line= element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size=2),
        axis.ticks.length = unit(.4, "cm"))
	print(q)
}

# -- Test SES significance
cor.test(sesdf$SES, sesdf$PC1, method=c("pearson"))
cor.test(sesdf$SES, sesdf$PC2, method=c("pearson"))
cor.test(sesdf$SES, sesdf$PC3, method=c("pearson"))
cor.test(sesdf$SES, sesdf$PC4, method=c("pearson"))
cor.test(sesdf$SES, sesdf$PC5, method=c("pearson"))

# -- Load SES (Hollingshead)
ses <- c(NA,69,51,54,65,61,58,18,61,11,58,33,61,37,65,65,29,43,NA,61,66,73,58,33,65,57,40,69,40,65,73,51,11,44,22,47,18,65,65,58,54,11,37,26,47,44,NA,44,26,69,NA,65,36,51,11,61,33,54,62,44,36,61,22,54,55,65,15,47,58,41,50,18,61,22,47,22,65,44,61,61,55,69,29,57,61,54,54,73,69,29,33,36,54,44,52,29,33,58,40,NA,69,50,65,69,18,40,61,62,40,58,58,22,57,44,44,37,65,61,48,NA,61,51,61,18,65,37,57,51,51,54,69,26,44,54,44,61,61,40,58,37,15,69,61,33,18,58,73,26,44,58,65,65,51,29,26,65,22,61,36,22,61,69,65,54,40,22,NA,54,37,40,47,NA,NA,65,44,33,57,73,NA,33,22,61,69,25,62,65,58,40,58,40,33,57,47,11,65,58,36,44,40,26,55,62,NA,36,22,22,57,51,69,65,40,69,61,40,22,69,65,61,40,54,33,54,36,65,40,44,29,NA,33,69,NA,40,44,47,44,65,69,33,44,47,33,29,44,54,18,NA,NA,25,58,44,54,44,69,18,61,58,29,73,NA,47,62,69,NA,58,69,65,29,40,65,18,69,65,40,44,22,15,58,57,61,36,33,65,65,61,18,65,18,54,25,54,40,69,65,47,33,22,51,61,69,33,65,22,62,62,57,65,33,33,73,58,61,NA,61,69,37,29,36,22,25,44,65,29,61,65,65,26,11,61,44,65,58,22,51,26,51,53,37,29,58,61,65,40,58,36,47,NA,65,44,65,73,22,54,54,61,44,29,40,65,58,58,25,65,61,44,40,69,65,54,61,41,33,44,40,51,40,29,22,33,NA,65,NA,61,61,40,58,54,47,29,44,29,NA,33,58,29,33,43,73,65,61,40,51,40,29,61,54,59,69,NA,33,61,43,65,40,40,NA,29,58,40,22,NA,69,NA,58,NA,44,40,59,36,54,66,59,44,29,54,NA,18)
sesclass <- c(NA,5,4,4,5,4,4,2,4,1,4,3,4,3,5,5,2,3,NA,4,5,5,4,3,5,4,3,5,3,5,5,4,1,3,2,3,2,5,5,4,4,1,3,2,3,3,NA,3,2,5,NA,5,3,4,1,4,3,4,4,3,3,4,2,4,4,5,1,3,4,3,4,2,4,2,3,2,5,3,4,4,4,5,2,4,4,4,4,5,5,2,3,3,4,3,4,2,3,4,3,NA,5,4,5,5,2,3,4,4,3,4,4,2,4,3,3,3,5,4,4,NA,4,4,4,2,5,3,4,4,4,4,5,2,3,4,3,4,4,3,4,3,1,5,4,3,2,4,5,2,3,4,5,5,4,2,2,5,2,4,3,2,4,5,5,4,3,2,NA,4,3,3,3,NA,NA,5,3,3,4,5,NA,3,2,4,5,2,4,5,4,3,4,3,3,4,3,1,5,4,3,3,3,2,4,4,NA,3,2,2,4,4,5,5,3,5,4,3,2,5,5,4,3,4,3,4,3,5,3,3,2,NA,3,5,NA,3,3,3,3,5,5,3,3,3,3,2,3,4,2,NA,NA,2,4,3,4,3,5,2,4,4,2,5,NA,3,4,5,NA,4,5,5,2,3,5,2,5,5,3,3,2,1,4,4,4,3,3,5,5,4,2,5,2,4,2,4,3,5,5,3,3,2,4,4,5,3,5,2,4,4,4,5,3,3,5,4,4,NA,4,5,3,2,3,2,2,3,5,2,4,5,5,2,1,4,3,5,4,2,4,2,4,4,3,2,4,4,5,3,4,3,3,NA,5,3,5,5,2,4,4,4,3,2,3,5,4,4,2,5,4,3,3,5,5,4,4,3,3,3,3,4,3,2,2,3,NA,5,NA,4,4,3,4,4,3,2,3,2,NA,3,4,2,3,3,5,5,4,3,4,3,2,4,4,4,5,NA,3,4,3,5,3,3,NA,2,4,3,2,NA,5,NA,4,NA,3,3,4,3,4,5,4,3,2,4,NA,2)

sesdf  <- data.frame("SES"=ses,  "PC1"=pcall$x[,1],"PC2"=pcall$x[,2],"PC3"=pcall$x[,3],"PC4"=pcall$x[,4],"PC5"=pcall$x[,5])
sesclassdf  <- data.frame("SES"=sesclass, "PC1"=pcall$x[,1],"PC2"=pcall$x[,2],"PC3"=pcall$x[,3],"PC4"=pcall$x[,4],"PC5"=pcall$x[,5])
sesclassdf$SES <- as.factor(sesclassdf$SES)
Group=c("SADP","SCZP","BPP","SCZP","SADP","SCZP","SCZP","BPP","BPP","BPP","SCZP","BPP","SCZP","BPP","SCZP","SCZP","BPP","BPP","BPP","SADP","SCZP","SCZP","SCZP","BPP","BPP","SADP","BPP","SADP","BPP","SCZP","SADP","BPP","SADP","SCZP","BPP","BPP","SADP","SADP","BPP","SCZP","SCZP","SADP","BPP","SADP","SCZP","SADP","SADP","SADP","SADP","SCZP","SCZP","SADP","SADP","SCZP","BPP","SCZP","BPP","SCZP","SCZP","SADP","BPP","BPP","BPP","BPP","SADP","SCZP","BPP","BPP","SADP","BPP","BPP","BPP","BPP","SCZP","BPP","BPP","SADP","SCZP","SADP","SADP","SCZP","SADP","SCZP","SADP","BPP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","SADP","SCZP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","BPP","SADP","SADP","SCZP","SCZP","BPP","BPP","SCZP","SADP","SCZP","SCZP","SCZP","BPP","SCZP","SADP","SADP","BPP","BPP","SCZP","SCZP","BPP","SCZP","BPP","SADP","SCZP","SADP","BPP","SADP","SADP","BPP","BPP","SADP","SADP","BPP","BPP","SADP","BPP","BPP","SCZP","SCZP","SCZP","BPP","BPP","BPP","BPP","SADP","SADP","SADP","SCZP","SCZP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SCZP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SADP","SCZP","BPP","BPP","SCZP","SADP","SCZP","SADP","SCZP","BPP","SCZP","BPP","BPP","SADP","BPP","BPP","BPP","BPP","BPP","BPP","BPP","BPP","SADP","SCZP","SCZP","BPP","SADP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SADP","SCZP","SCZP","BPP","BPP","BPP","BPP","BPP","BPP","SADP","SCZP","SADP","SADP","SCZP","SCZP","BPP","SCZP","BPP","BPP","SADP","BPP","SADP","SCZP","SADP","SCZP","SADP","SCZP","BPP","SCZP","SCZP","BPP","SADP","SADP","SADP","BPP","SADP","SCZP","SADP","SADP","BPP","SADP","SCZP","SCZP","SCZP","BPP","SCZP","SCZP","SCZP","SCZP","SADP","BPP","SCZP","SCZP","SADP","BPP","SCZP","BPP","BPP","SCZP","BPP","BPP","SCZP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SCZP","SADP","SCZP","SCZP","SCZP","SADP","SCZP","BPP","SCZP","SCZP","SADP","SADP","BPP","BPP","SADP","SADP","BPP","BPP","SCZP","SCZP","SADP","SADP","BPP","BPP","BPP","BPP","SCZP","SCZP","SCZP","SCZP","SCZP","BPP","SADP","SCZP","SCZP","SCZP","BPP","SADP","SCZP","SADP","SCZP","SCZP","BPP","BPP","SCZP","SCZP","SADP","BPP","SADP","BPP","SADP","SADP","SADP","SADP","SCZP","BPP","SADP","BPP","SCZP","SADP","BPP","BPP","BPP","SCZP","SADP","SADP","SADP","SCZP","BPP","SCZP","SCZP","SCZP","SADP","BPP","BPP","SCZP","SCZP","BPP","SADP","SADP","SCZP","SADP","SADP","SADP","SCZP","SADP","SCZP","SADP","SADP","SADP","SADP","SCZP","BPP","SCZP","BPP","SADP","SCZP","SCZP","BPP","SADP","SCZP","BPP","SCZP","BPP","SCZP","SADP","BPP","SCZP","SCZP","BPP","SCZP","BPP","BPP","BPP","SADP","SCZP","BPP","SCZP","BPP","SCZP","BPP","SCZP","SADP","SCZP","SADP","BPP","SADP","SCZP","SCZP","SCZP","SADP","BPP","SADP","BPP","SCZP","BPP","BPP","BPP","BPP","SADP","SADP","SADP","SCZP","SADP","SCZP","BPP","SADP","SCZP","BPP","SCZP","BPP","BPP","BPP","SCZP")

for (pcno in 1:5){
	q <- ggplot(sesdf, aes(col=Group, x=SES, y=get(paste0("PC",pcno)))) +
	geom_point(aes(col=Group),size=3) +
  scale_colour_manual(values = c("CON"="white","PROB"="black","BPP"="#feb24c","SADP"="#fc4e2a","SCZP"="#800026")) +
	scale_x_continuous(limits=c(10,82)) +
  scale_y_continuous(limits=c(-10,6)) +
	stat_smooth(method = "lm", size = 3, se=FALSE, col="red",fullrange = FALSE,na.rm=TRUE) +
	xlab("SES") +
	ylab(paste0("PC",pcno, " Score")) +
  theme1+
  theme(axis.text.x = element_text(size = 35, color="black", margin=margin(5,0,0,0), vjust=1), 
        axis.title.x = element_text(size = 35, color="black"),
        axis.text.y = element_text(size = 35, color="black", angle=90, hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 35,margin=margin(0,5,0,0)),
        axis.line= element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size=2),
        axis.ticks.length = unit(.4, "cm"))
	print(q)
}

# -- Test SES significance
cor.test(sesdf$SES, sesdf$PC1, method=c("pearson"))
cor.test(sesdf$SES, sesdf$PC2, method=c("pearson"))
cor.test(sesdf$SES, sesdf$PC3, method=c("pearson"))
cor.test(sesdf$SES, sesdf$PC4, method=c("pearson"))
cor.test(sesdf$SES, sesdf$PC5, method=c("pearson"))

sesclassdf_m <- aggregate(sesclassdf, by=list(sesclassdf$SES), FUN="mean")[,-2]

sesPC1df_summ <- summarySE(sesclassdf, measurevar="PC1", groupvars=c("SES"))
sesPC2df_summ <- summarySE(sesclassdf, measurevar="PC2", groupvars=c("SES"))
sesPC3df_summ <- summarySE(sesclassdf, measurevar="PC3", groupvars=c("SES"))
sesPC4df_summ <- summarySE(sesclassdf, measurevar="PC4", groupvars=c("SES"))
sesPC5df_summ <- summarySE(sesclassdf, measurevar="PC5", groupvars=c("SES"))

for (pcno in 1:5){
q <-  ggplot(data=get(paste0("sesPC",pcno,"df_summ")), aes(x=SES,y=get(paste0("PC",pcno)))) +
      geom_col(col="grey50",fill="grey50") +
      scale_x_discrete(limits=c("1","2","3","4","5"),labels=c("1","2","3","4","5")) +
      theme1 +
      geom_errorbar(aes(x=SES, ymin=get(paste0("PC",pcno))-se, ymax=get(paste0("PC",pcno))+se), col="black", width=0.2, position=position_dodge(0.9), lwd=1.5) +  
      scale_y_continuous(expand = c(0.001, 0), limits=c(-2,2.4), breaks=seq(-2,2,0.5),labels=format(c(-2.0,"","","",0.0,"","","",2.0), nsmall=1)) +
      xlab("SES Class") +
      ylab(paste0("PC",pcno," Score")) +
      theme(axis.text.x = element_text(size = 35, color="black", margin=margin(5,0,0,0), vjust=1), 
        axis.title.x = element_text(size = 35, color="black"),
        axis.text.y = element_text(size = 35, color="black", angle=90, hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 35,margin=margin(0,5,0,0)),
        axis.line= element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size=2),
        axis.ticks.length = unit(.4, "cm"))    
  print(q)
}

# -- Test significance

sesnonaclassdf <- na.omit(sesclassdf)
summary(aov(PC1 ~ SES, data = sesnonaclassdf))
summary(aov(PC2 ~ SES, data = sesnonaclassdf))
summary(aov(PC3 ~ SES, data = sesnonaclassdf))
summary(aov(PC4 ~ SES, data = sesnonaclassdf))
summary(aov(PC5 ~ SES, data = sesnonaclassdf))

for (pcno in 1:5){
	q <- ggplot(sesnonaclassdf, aes(x=get(paste0("PC",pcno)),y=SES, group=SES)) +
  geom_density_ridges(fill=c("grey50"), lwd=2, rel_min_height = 0.0001) +
  geom_segment(aes(x = eval(parse(text=paste0("sesclassdf_m$PC", pcno)))[1], y=1.0,xend = eval(parse(text=paste0("sesclassdf_m$PC", pcno)))[1], yend = 2.0),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = eval(parse(text=paste0("sesclassdf_m$PC", pcno)))[2], y=2.0,xend = eval(parse(text=paste0("sesclassdf_m$PC", pcno)))[2], yend = 3.0),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = eval(parse(text=paste0("sesclassdf_m$PC", pcno)))[3], y=3.0,xend = eval(parse(text=paste0("sesclassdf_m$PC", pcno)))[3], yend = 4.0),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = eval(parse(text=paste0("sesclassdf_m$PC", pcno)))[4], y=4.0,xend = eval(parse(text=paste0("sesclassdf_m$PC", pcno)))[4], yend = 5.0),size=2, linetype = 1,color = "white") +
  geom_segment(aes(x = eval(parse(text=paste0("sesclassdf_m$PC", pcno)))[5], y=5.0,xend = eval(parse(text=paste0("sesclassdf_m$PC", pcno)))[5], yend = 6.0),size=2, linetype = 1,color = "white") +
	scale_y_discrete(limits=c(1,2,3,4,5), labels=c(1,2,3,4,5)) + 
  xlab(paste0("PC",pcno, " Score")) +
	ylab("SES Class") +
	theme1 + 
      theme(axis.text.x = element_text(size = 35, color="black", margin=margin(5,0,0,0), vjust=1), 
        axis.title.x = element_text(size = 35, color="black"),
        axis.text.y = element_text(size = 35, color="black", angle=90, hjust=0.5,margin=margin(0,5,0,0)), 
        axis.title.y = element_text(size = 35,margin=margin(0,5,0,0)),
        axis.line= element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size=2),
        axis.ticks.length = unit(.4, "cm"))    
	print(q)
}