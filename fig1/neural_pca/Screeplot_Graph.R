
install.packages('gplots')

library(ggplot2)
library(Hmisc)
library(gplots)
library(plyr)
library(scales)

# set working directory
setwd("~/Dropbox/DP5_N-Bridge/Figure_Calculations/GSR/N40_Neural_PCA_mVWMWB1d")


# ---------------------------------------------------------------------------------------->
### Neural pGBC PCA Screeplots + Var plots
# ---------------------------------------------------------------------------------------->

# read in dataframes: 1. explained variance 2. pcrit

DIFFpgbc = data.frame(cbind(read.table("ket-pla_pca_explainedvar.dat", header=FALSE)),read.table("ket-pla_pca_explainedvar_5kshuffle_pcrit95.dat", header=FALSE))

# insert column of PC 1-39 at begining of dataframe


DIFFpgbc["PC"] <- 1:nrow(DIFFpgbc) 
DIFFpgbc <- DIFFpgbc[c('PC', 'V1', 'V1.1')]

# rename columns
names(DIFFpgbc) <- c("PC", "% Variance","Pcrit")

# define theme
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
               #legend.title=element_blank(),
               #legend.text=element_text(size=20),
               #legend.key = element_blank(),
               #legend.background = element_rect(fill = "white", colour = "white"),
               #legend.box.background = element_rect(fill = "white", colour = "white"),
               legend.position="none",
               axis.ticks = element_line(colour = "black", size=2),
               axis.ticks.length = unit(.4, "cm"))


# -------------------->
# DELTA
# -------------------->


# difference scree plot
# ***N.B. need to choose those % variance values that are bigger than pcrit***

pdf("Difference_Neural_PCA_screeplot_15_mainfigure.pdf", width=10,height=7)
ggplot(data=DIFFpgbc[1:39,], aes(x=PC, y=`% Variance`)) +
  geom_line(aes(x=PC,y=Pcrit), lty=2, col="dark red",lwd=1.2) +
  geom_line(lwd=1.2)+
  geom_point(data=DIFFpgbc[1:5,],col="black",pch=21,bg="#666666",aes(size=`% Variance`)) + # choose those % variance values that are bigger than pcrit
  scale_size_area(name=DIFFpgbc$`% Variance`, max_size=18) +
  geom_point(size=4,data=DIFFpgbc[6:39,],col="black",pch=21,bg="white") +
  ylab('% Variance Explained') +
  xlab('Principal Component') +
  scale_x_continuous(limits=c(0,40), labels=c(1,39),breaks=seq(1,39,38)) + # breaks=seq(1,39,38) is tick mark showing every 38 spots starting at 1 ending at 39, 
  scale_y_continuous(limits=c(0,20), labels=c(0,5,10,15,20)) +
  #scale_y_continuous(expand = c(0.001, 0),limits=c(0,22),labels=c(0,5,10,15,20), breaks=c(0,5,10,15,20)) +
  #geom_text(data = label.df, label = "***", cex=9, nudge_y=0.01) +
  theme1 + 
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        axis.title.y = element_text(size = 35),
        axis.text.y = element_text(size = 30),
        axis.title.x = element_text(size = 35,margin=margin(-20,-30,0,0)),
        axis.text.x = element_text(size = 30, margin=margin(10,0,0,0)))
dev.off()




