# Scripts to produce the figures of
# Joly et al.
# Flexible methods for estimating genetic distance from nucleotide data 

require(ggplot2)
require(gridExtra)
require(plyr)

# Go to data directory

# summarySE ---------------------------------------------------------------

# function required for the script.

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
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
                     median = median   (xx[[col]], na.rm=na.rm),
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


# Figure 1 -----------------------------------------------------------------

results1 <- read.table("Fig1a_results.txt", header=TRUE, dec=".")
results1 <- as.data.frame(results1)
results2 <- read.table("Fig1b_results.txt", header=TRUE, dec=".")
results2 <- as.data.frame(results2)

# Get expected sequence divergence
results1$seqdiv <- results1[,3]*2+(0.001)
results2$seqdiv <- results2[,3]*2+(0.01)

#apply cutoff for plotting
results1[results1[,1]>0.06,1] <- 0.06
results2[results2[,1]>0.08,1] <- 0.08

colPalette <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f")

#plotting
p <- ggplot(results1, aes(factor(divergence), distance))
plot1 <- p + geom_boxplot(aes(fill = factor(divergence)),outlier.size=0.5,outlier.shape="") +
  geom_hline(aes(yintercept=(seqdiv),colour=factor(seqdiv)),size=0.5,linetype=2) + 
  facet_wrap(~ method,ncol=2) + ylim(0,0.06) + guides(fill=FALSE) + theme_bw() +
  xlab(expression(paste("Population divergence (",tau,")"))) + ylab("Genetic distance") +
  ggtitle(expression(paste(theta," = 0.001"))) +
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        plot.title = element_text(size=15,vjust=1.5)) +
  scale_fill_manual(values=colPalette) + scale_colour_manual(values=colPalette)

p <- ggplot(results2, aes(factor(divergence), distance))
plot2 <- p + geom_boxplot(aes(fill = factor(divergence)),outlier.size=0.5,outlier.shape="") +
  geom_hline(aes(yintercept=(seqdiv),colour=factor(seqdiv)),size=0.5,linetype=2) + 
  facet_wrap(~ method,ncol=2) + ylim(0,0.08) + guides(fill=FALSE) + theme_bw() +
  xlab(expression(paste("Population divergence (",tau,")"))) + ylab("Genetic distance") + 
  ggtitle(expression(paste(theta," = 0.01"))) +
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        plot.title = element_text(size=15,vjust=1.5)) + 
  scale_fill_manual(values=colPalette) + scale_colour_manual(values=colPalette)

pdf("Fig1.pdf",width=12,height=9)
grid.arrange(plot1, plot2, ncol=2)
dev.off()
#postscript("Fig1.ps",width=12,height=9)
#grid.arrange(plot1, plot2, ncol=2)
#dev.off()


# Figure 2 ---------------------------------------------------------------

results1 <- read.table("Fig1a_results.txt", header=TRUE, dec=".")
results1 <- as.data.frame(results1)
results2 <- read.table("Fig1b_results.txt", header=TRUE, dec=".")
results2 <- as.data.frame(results2)

# Get expected sequence divergence
results1$seqdiv <- results1[,3]*2+(0.001)
results2$seqdiv <- results2[,3]*2+(0.01)

colPalette7 <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f",
                 "#e5c494")

res1.med <- summarySE(results1, measurevar="distance", groupvars=c("method","divergence"))
res1.med <- data.frame(res1.med,theta=rep("paste(theta,\" = 0.001\")",nrow(res1.med)))
polygon1 <- data.frame(x=c(0,0.002,0.004,0.008,0.012,0.02,0.02,0.012,0.008,0.004,0.002,0),
                       y=c(0,0.004,0.008,0.016,0.024,0.04,0.041,0.025,0.017,0.009,0.005,0.001))
p <- ggplot(res1.med, aes(x=divergence,y=median))
plot1c <- p + geom_polygon(data=polygon1,aes(x=x,y=y), fill="gray85", linetype=0, inherit.aes = FALSE) +
  geom_line(aes(colour=factor(method),linetype=factor(method)),size=0.6) +  
  xlab(expression(paste("Population divergence (",tau,")"))) +
  ylab("Median genetic distance") + 
  facet_grid(~ theta,labeller=label_parsed) + theme_bw() + 
  guides(colour = guide_legend(title = "Distance method", override.aes = list(size=0.8)), 
         linetype  = guide_legend(title = "Distance method")) +
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        legend.key=element_rect(linetype=0),legend.key.width=unit(2,units="cm")) +
  scale_colour_manual(values=colPalette7)

res2.med <- summarySE(results2, measurevar="distance", groupvars=c("method","divergence"))
res2.med <- data.frame(res2.med,theta=rep("paste(theta,\" = 0.01\")",nrow(res2.med)))
polygon2 <- data.frame(x=c(0,0.002,0.004,0.008,0.012,0.02,0.02,0.012,0.008,0.004,0.002,0),
                       y=c(0,0.004,0.008,0.016,0.024,0.04,0.05,0.034,0.026,0.018,0.014,0.01))
p <- ggplot(res2.med, aes(x=divergence,y=median))
plot1d <- p + geom_polygon(data=polygon2,aes(x=x,y=y), fill="gray85", linetype=0, inherit.aes = FALSE) +
  geom_line(aes(colour=factor(method),linetype=factor(method)),size=0.6) +
  xlab(expression(paste("Population divergence (",tau,")"))) + 
  ylab("Median genetic distance") +
  facet_grid(~ theta,labeller=label_parsed) + theme_bw() + 
  guides(colour = guide_legend(title = "Distance method", override.aes = list(size=0.8)), 
         linetype  = guide_legend(title = "Distance method")) +
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        legend.key=element_rect(linetype=0),legend.key.width=unit(2,units="cm")) +
  scale_colour_manual(values=colPalette7)

pdf("Fig2.pdf",width=8,height=9)
grid.arrange(plot1c, plot1d, ncol=1)
dev.off()

#postscript("Fig2.ps",width=8,height=9)
#grid.arrange(plot1c, plot1d, ncol=1)
#dev.off()


# Figure S1 -----------------------------------------------------------------

# Load data with rate variation
resultsS1 <- read.table("FigS1a_results.txt", header=TRUE, dec=".")
resultsS1 <- as.data.frame(resultsS1)
resultsS2 <- read.table("FigS1b_results.txt", header=TRUE, dec=".")
resultsS2 <- as.data.frame(resultsS2)

# Get expected sequence divergence
resultsS1$seqdiv <- resultsS1[,3]*2+(0.001)
resultsS2$seqdiv <- resultsS2[,3]*2+(0.01)

# apply cutoff for plotting
resultsS1[resultsS1[,1]>0.06,1] <- 0.06
resultsS2[resultsS2[,1]>0.08,1] <- 0.08

colPalette <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f")

#plotting
p <- ggplot(resultsS1, aes(factor(divergence), distance))
plotS1 <- p + geom_boxplot(aes(fill = factor(divergence)),outlier.size=0.5,outlier.shape="") +
  geom_hline(aes(yintercept=(seqdiv),colour=factor(seqdiv)),size=0.5,linetype=2) + 
  facet_wrap(~ method,ncol=2) + ylim(0,0.06) + guides(fill=FALSE) + theme_bw() +
  xlab(expression(paste("Population divergence (",tau,")"))) + ylab("Genetic distance") +
  ggtitle(expression(paste(theta," = 0.001"))) +
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        plot.title = element_text(size=15,vjust=1.5)) +
  scale_fill_manual(values=colPalette) + scale_colour_manual(values=colPalette)

p <- ggplot(resultsS2, aes(factor(divergence), distance))
plotS2 <- p + geom_boxplot(aes(fill = factor(divergence)),outlier.size=0.5,outlier.shape="") +
  geom_hline(aes(yintercept=(seqdiv),colour=factor(seqdiv)),size=0.5,linetype=2) + 
  facet_wrap(~ method,ncol=2) + ylim(0,0.08) + guides(fill=FALSE) + theme_bw() +
  xlab(expression(paste("Population divergence (",tau,")"))) + ylab("Genetic distance") + 
  ggtitle(expression(paste(theta," = 0.01"))) +
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        plot.title = element_text(size=15,vjust=1.5)) + 
  scale_fill_manual(values=colPalette) + scale_colour_manual(values=colPalette)

pdf("FigS1.pdf",width=12,height=9)
grid.arrange(plotS1, plotS2, ncol=2)
dev.off()

#postscript("FigS1.ps",width=12,height=9)
#grid.arrange(plotS1, plotS2, ncol=2)
#dev.off()


# Figure 3 ---------------------------------------------------------------


# Reload data without rate variation
results1 <- read.table("Fig1a_results.txt", header=TRUE, dec=".")
results1 <- as.data.frame(results1)
results2 <- read.table("Fig1b_results.txt", header=TRUE, dec=".")
results2 <- as.data.frame(results2)

# Get expected sequence divergence
results1$seqdiv <- results1[,3]*2+(0.001)
results2$seqdiv <- results2[,3]*2+(0.01)

# reload data with rate variation
resultsS1 <- read.table("FigS1a_results.txt", header=TRUE, dec=".")
resultsS1 <- as.data.frame(resultsS1)
resultsS2 <- read.table("FigS1b_results.txt", header=TRUE, dec=".")
resultsS2 <- as.data.frame(resultsS2)

# Get expected sequence divergence
resultsS1$seqdiv <- resultsS1[,3]*2+(0.001)
resultsS2$seqdiv <- resultsS2[,3]*2+(0.01)


# Summarize the data
x1 <- ddply(results1, .(method, divergence), summarize,
            sd = round(sd(distance),5))
x2 <- ddply(results2, .(method, divergence), summarize,
            sd = round(sd(distance),5))
x1<-cbind("no",0.001,x1)
colnames(x1)[1:2]<-c("Rate_variation","theta")
x2<-cbind("no",0.01,x2)
colnames(x2)[1:2]<-c("Rate_variation","theta")

xS1 <- ddply(resultsS1, .(method, divergence), summarize,
            sd = round(sd(distance),5))
xS2 <- ddply(resultsS2, .(method, divergence), summarize,
            sd = round(sd(distance),5))
xS1<-cbind("yes",0.001,xS1)
colnames(xS1)[1:2]<-c("Rate_variation","theta")
xS2<-cbind("yes",0.01,xS2)
colnames(xS2)[1:2]<-c("Rate_variation","theta")
x <- rbind(x1,x2,xS1,xS2)
p <- ggplot(x, aes(divergence, sd))
plot3 <- p + geom_line(aes(colour=factor(theta),linetype=Rate_variation),size=0.6) +
  facet_wrap(~method,nrow=2) + theme_bw() + scale_colour_manual(values=colPalette) +
  xlab(expression(paste("Population divergence (",tau,")"))) + 
  ylab("Standard deviation of sequence divergence") + 
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        legend.key=element_rect(linetype=0)) + 
  guides(colour = guide_legend(title = "Theta", override.aes = list(size=1.2)), 
         linetype  = guide_legend(title = "Rate variation"))

pdf("Fig3.pdf",width=12,height=6)
plot3
dev.off()


# Figure S2 -----------------------------------------------

results1 <- read.table("recomb.short.a_results.txt", header=TRUE, dec=".")
results1 <- as.data.frame(results1)
results2 <- read.table("recomb.short.b_results.txt", header=TRUE, dec=".")
results2 <- as.data.frame(results2)

# Get expected sequence divergence
results1$seqdiv <- results1[,3]*2+(0.001)
results2$seqdiv <- results2[,3]*2+(0.01)

#apply cutoff for plotting
results1[results1[,1]>0.06,1] <- 0.06
results2[results2[,1]>0.08,1] <- 0.08

colPalette <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f")

#plotting
p <- ggplot(results1, aes(factor(divergence), distance))
plot1 <- p + geom_boxplot(aes(fill = factor(divergence)),outlier.size=0.5,outlier.shape="") +
  geom_hline(aes(yintercept=(seqdiv),colour=factor(seqdiv)),size=0.5,linetype=2) + 
  facet_wrap(~ method,ncol=2) + ylim(0,0.06) + guides(fill=FALSE) + theme_bw() +
  xlab(expression(paste("Population divergence (",tau,")"))) + ylab("Genetic distance") +
  ggtitle(expression(paste(theta," = 0.001"))) +
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        plot.title = element_text(size=15,vjust=1.5)) +
  scale_fill_manual(values=colPalette) + scale_colour_manual(values=colPalette)

p <- ggplot(results2, aes(factor(divergence), distance))
plot2 <- p + geom_boxplot(aes(fill = factor(divergence)),outlier.size=0.5,outlier.shape="") +
  geom_hline(aes(yintercept=(seqdiv),colour=factor(seqdiv)),size=0.5,linetype=2) + 
  facet_wrap(~ method,ncol=2) + ylim(0,0.08) + guides(fill=FALSE) + theme_bw() +
  xlab(expression(paste("Population divergence (",tau,")"))) + ylab("Genetic distance") + 
  ggtitle(expression(paste(theta," = 0.01"))) +
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        plot.title = element_text(size=15,vjust=1.5)) + 
  scale_fill_manual(values=colPalette) + scale_colour_manual(values=colPalette)

pdf("FigS2.pdf",width=12,height=9)
grid.arrange(plot1, plot2, ncol=2)
dev.off()


# Figure 4 ------------------------------------------------

# Recombination

results1 <- read.table("recomb.none_results.txt", header=TRUE, dec=".")
results1 <- as.data.frame(results1)
results2 <- read.table("recomb.high_results.txt", header=TRUE, dec=".")
results2 <- as.data.frame(results2)

# Get expected sequence divergence
results1$seqdiv <- results1[,3]*2+(0.001)
results2$seqdiv <- results2[,3]*2+(0.001)

#qqplot with ggplot

variable=2
all.results <- rbind(
  data.frame(
    method="GENPOFAD",
    no_recomb=results1[intersect(which(results1[,2]=="GENPOFAD"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="GENPOFAD"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="GENPOFAD"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="GENPOFAD"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),
  data.frame(
    method="MATCHSTATES",
    no_recomb=results1[intersect(which(results1[,2]=="MATCHSTATES"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="MATCHSTATES"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="MATCHSTATES"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="MATCHSTATES"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="MIN",
    no_recomb=results1[intersect(which(results1[,2]=="MIN"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="MIN"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="MIN"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="MIN"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="MRCA",
    no_recomb=results1[intersect(which(results1[,2]=="MRCA"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="MRCA"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="MRCA"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="MRCA"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="NEI",
    no_recomb=results1[intersect(which(results1[,2]=="NEI"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="NEI"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="NEI"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="NEI"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="PBC",
    no_recomb=results1[intersect(which(results1[,2]=="PBC"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="PBC"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="PBC"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="PBC"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="PP",
    no_recomb=results1[intersect(which(results1[,2]=="PP"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="PP"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="PP"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="PP"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  )
)
p <- ggplot(all.results,aes(no_recomb,recomb))
plot_recomb1 <- p + 
  geom_vline(xintercept=0.005,colour="gray80",size=1.2) +
  geom_point(colour="#8da0cb") + 
  facet_wrap(~method,ncol=1) + theme_bw() + 
  xlab(expression(paste(italic(r)," = ", 0))) + 
  ylab(expression(paste(italic(r)," = ", 2^-8))) + 
  ggtitle(expression(paste(tau," = 0.002"))) + theme_bw() +
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8),
        plot.title = element_text(size=15,vjust=1.5)) + 
  geom_abline(intercept=0, slope=1)

variable=4
all.results <- rbind(
  data.frame(
    method="GENPOFAD",
    no_recomb=results1[intersect(which(results1[,2]=="GENPOFAD"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="GENPOFAD"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="GENPOFAD"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="GENPOFAD"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),
  data.frame(
    method="MATCHSTATES",
    no_recomb=results1[intersect(which(results1[,2]=="MATCHSTATES"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="MATCHSTATES"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="MATCHSTATES"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="MATCHSTATES"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="MIN",
    no_recomb=results1[intersect(which(results1[,2]=="MIN"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="MIN"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="MIN"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="MIN"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="MRCA",
    no_recomb=results1[intersect(which(results1[,2]=="MRCA"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="MRCA"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="MRCA"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="MRCA"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="NEI",
    no_recomb=results1[intersect(which(results1[,2]=="NEI"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="NEI"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="NEI"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="NEI"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="PBC",
    no_recomb=results1[intersect(which(results1[,2]=="PBC"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="PBC"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="PBC"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="PBC"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="PP",
    no_recomb=results1[intersect(which(results1[,2]=="PP"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="PP"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="PP"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="PP"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  )
)
p <- ggplot(all.results,aes(no_recomb,recomb))
plot_recomb2 <- p + 
  geom_vline(xintercept=0.017,colour="gray80",size=1.2) +
  geom_point(colour="#8da0cb") + 
  facet_wrap(~method,ncol=1)  + theme_bw() +
  xlab(expression(paste(italic(r)," = ", 0))) + 
  ylab(expression(paste(italic(r)," = ", 2^-8))) + 
  ggtitle(expression(paste(tau," = 0.008"))) + theme_bw() +
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8),
        plot.title = element_text(size=15,vjust=1.5)) + 
  geom_abline(intercept=0, slope=1)

variable=6
all.results <- rbind(
  data.frame(
    method="GENPOFAD",
    no_recomb=results1[intersect(which(results1[,2]=="GENPOFAD"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="GENPOFAD"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="GENPOFAD"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="GENPOFAD"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),
  data.frame(
    method="MATCHSTATES",
    no_recomb=results1[intersect(which(results1[,2]=="MATCHSTATES"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="MATCHSTATES"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="MATCHSTATES"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="MATCHSTATES"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="MIN",
    no_recomb=results1[intersect(which(results1[,2]=="MIN"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="MIN"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="MIN"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="MIN"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="MRCA",
    no_recomb=results1[intersect(which(results1[,2]=="MRCA"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="MRCA"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="MRCA"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="MRCA"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="NEI",
    no_recomb=results1[intersect(which(results1[,2]=="NEI"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="NEI"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="NEI"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="NEI"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="PBC",
    no_recomb=results1[intersect(which(results1[,2]=="PBC"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="PBC"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="PBC"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="PBC"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  ),  
  data.frame(
    method="PP",
    no_recomb=results1[intersect(which(results1[,2]=="PP"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results1[intersect(which(results1[,2]=="PP"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])],
    recomb=results2[intersect(which(results2[,2]=="PP"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1][order(results2[intersect(which(results2[,2]=="PP"),seq(from=variable,by=6,length.out=nrow(results1)/6)),1])]
  )
)
p <- ggplot(all.results,aes(no_recomb,recomb))
plot_recomb3 <- p + 
  geom_vline(xintercept=0.041,colour="gray80",size=1.2) +
  geom_point(colour="#8da0cb") + 
  facet_wrap(~method,ncol=1) + 
  xlab(expression(paste(italic(r)," = ", 0))) + 
  ylab(expression(paste(italic(r)," = ", 2^-8))) + 
  ggtitle(expression(paste(tau," = 0.02"))) + theme_bw() +
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8),
        plot.title = element_text(size=15,vjust=1.5)) + 
  geom_abline(intercept=0, slope=1)

pdf("Fig4.pdf",width=8,height=12)
grid.arrange(plot_recomb1, plot_recomb2, plot_recomb3, ncol=3)
dev.off()


# Figure 5 -----------------------------------------------------------------

Fig5a_old <- read.table("Fig5a.old_results.txt", header=TRUE, dec=".")
Fig5a_old <- as.data.frame(cbind(Fig5a_old,age=rep("old",nrow(Fig5a_old))))
Fig5a_mid <- read.table("Fig5a.mid_results.txt", header=TRUE, dec=".")
Fig5a_mid <- as.data.frame(cbind(Fig5a_mid,age=rep("intermediate",nrow(Fig5a_mid)))) 
Fig5a_young <- read.table("Fig5a.young_results.txt", header=TRUE, dec=".")
Fig5a_young <- as.data.frame(cbind(Fig5a_young,age=rep("young",nrow(Fig5a_young)))) 

Fig5a <- rbind(Fig5a_young,Fig5a_mid,Fig5a_old)
#remove na values (due to divide by 0)
Fig5a[is.na(Fig5a[,1]),1] <- 0.5

Fig5b_old <- read.table("Fig5b.old_results.txt", header=TRUE, dec=".")
Fig5b_old <- as.data.frame(cbind(Fig5b_old,age=rep("old",nrow(Fig5b_old))))
Fig5b_mid <- read.table("Fig5b.mid_results.txt", header=TRUE, dec=".")
Fig5b_mid <- as.data.frame(cbind(Fig5b_mid,age=rep("intermediate",nrow(Fig5b_mid)))) 
Fig5b_young <- read.table("Fig5b.young_results.txt", header=TRUE, dec=".")
Fig5b_young <- as.data.frame(cbind(Fig5b_young,age=rep("young",nrow(Fig5b_young)))) 

Fig5b <- rbind(Fig5b_young,Fig5b_mid,Fig5b_old)
#remove na values (due to divide by 0)
Fig5b[is.na(Fig5b[,1]),1] <- 0.5

p <- ggplot(Fig5a, aes(factor(method), index))
Plot3a <- p  + geom_hline(yintercept=0.5, linetype=1, size=2, colour="grey85") +
  geom_boxplot(aes(fill = factor(method)),outlier.size=1,outlier.shape="") + 
  ylim(0,1) + facet_wrap(~ age,ncol=1) + coord_flip() + guides(fill=FALSE) +
  xlab("") + ylab("Hybrid index") + theme_bw() + scale_fill_manual(values=colPalette7) +
  scale_colour_manual(values=colPalette7)

p <- ggplot(Fig5b, aes(factor(method), index))
Plot3b <- p  + geom_hline(yintercept=1/3, linetype=1, size=2, colour="grey85") +
  geom_boxplot(aes(fill = factor(method)),outlier.size=1,outlier.shape="") + 
  ylim(0,1) + facet_wrap(~ age,ncol=1) + coord_flip() + guides(fill=FALSE) +
  xlab("") + ylab("Hybrid index") + theme_bw() + scale_fill_manual(values=colPalette7) +
  scale_colour_manual(values=colPalette7)

grid.arrange(Plot3a, Plot3b, ncol=2)

pdf("Fig5.pdf",width=10,height=8)
grid.arrange(Plot3a, Plot3b, ncol=2)
dev.off()



# Figure 6 ---------------------------------------------------------------

# Monte Carlo sampling

results1 <- read.table("Fig1a_results.txt", header=TRUE, dec=".")
results1 <- as.data.frame(results1)
results1$seqdiv <- results1[,3]*2+(0.001)
Pres1 <- results1[results1$seqdiv==0.025,]
nmethods <- length(levels(Pres1[,2]))
nreplicates <- c(1,2,5,10,20,40)
nreps <- 500

sd.res <- numeric(0)
method.res<-character(0)
nreplicates.res <- numeric(0)

for (i in 1:nmethods) {
  for (j in 1:length(nreplicates)) {
    for (k in 1:nreps) {
      sd.res <- c(sd.res,mean(sample(Pres1[Pres1[,2]==levels(Pres1[,2])[i],1],nreplicates[j])))
      method.res <- c(method.res,levels(Pres1[,2])[i])
      nreplicates.res <- c(nreplicates.res,nreplicates[j])
    }
  }
}

Pres1.sd <- data.frame(mean=sd.res,method=method.res,nreplicates=nreplicates.res)
dfc <- summarySE(Pres1.sd, measurevar="mean", groupvars=c("method","nreplicates"))
dfc1 <- as.data.frame(cbind(dfc,exp=rep("Distance",nrow(dfc)))) 
colPalette7 <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f",
                 "#e5c494")
p <- ggplot(dfc1, aes(nreplicates, sd))
Plot6a <- p + geom_line(aes(colour=factor(method),linetype=factor(method)),size=0.6) +  xlab("Number of markers") +
  ylab("Standard deviation") + labs(colour="Method") + facet_wrap(~exp) + theme_bw()  +
  guides(colour = guide_legend(title = "Distance method", override.aes = list(size=0.8)), 
         linetype  = guide_legend(title = "Distance method")) +
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        legend.key=element_rect(linetype=0),legend.key.width=unit(1,units="cm")) +
  scale_colour_manual(values=colPalette7)


#open hybrid index
Pres2<-read.table("Fig5a.mid_results.txt", header=TRUE, dec=".")
Pres2 <- as.data.frame(Pres2)
Pres2[is.na(Pres2)] <- 0.5
nmethods <- length(levels(Pres2[,2]))
nreplicates <- c(1,2,5,10,20,40)
nreps <- 500

sd.res <- numeric(0)
method.res<-character(0)
nreplicates.res <- numeric(0)

for (i in 1:nmethods) {
  for (j in 1:length(nreplicates)) {
    for (k in 1:nreps) {
      #results.sd <- rbind(results.sd,c(sd(sample(res[res[,2]==levels(res[,2])[i],1],nreplicates[j])),levels(res[,2])[i],nreplicates[j]))
      sd.res <- c(sd.res,mean(sample(Pres2[Pres2[,2]==levels(Pres2[,2])[i],1],nreplicates[j])))
      method.res <- c(method.res,levels(Pres2[,2])[i])
      nreplicates.res <- c(nreplicates.res,nreplicates[j])
    }
  }
}

Pres2.sd <- data.frame(mean=sd.res,method=method.res,nreplicates=nreplicates.res)
dfc <- summarySE(Pres2.sd, measurevar="mean", groupvars=c("method","nreplicates"))
dfc2 <- as.data.frame(cbind(dfc,exp=rep("Hybrid index",nrow(dfc)))) 
p <- ggplot(dfc2, aes(nreplicates, sd))
Plot6b <- p + geom_line(aes(colour=factor(method),linetype=factor(method)),size=0.6) +
  xlab("Number of markers") +
  ylab("Standard deviation") + labs(colour="Method") + facet_wrap(~exp) + theme_bw() + 
  guides(colour = guide_legend(title = "Distance method", override.aes = list(size=0.8)), 
         linetype  = guide_legend(title = "Distance method")) +
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.2),
        legend.key=element_rect(linetype=0),legend.key.width=unit(1,units="cm")) +
  scale_colour_manual(values=colPalette7)

pdf("Fig6.pdf",width=5,height=7)
grid.arrange(Plot6a, Plot6b, ncol=1)
dev.off()

