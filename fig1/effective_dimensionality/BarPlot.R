
install.packages('gplots')

library(ggplot2)
library(Hmisc)
library(gplots)
library(plyr)
library(scales)

# set working directory
setwd("~/Dropbox/DP5_N-Bridge/Figure_Calculations/GSR/effective_dimensionality")

# read dataframes
df = data.frame(read.csv("participationratio_19_7_21.csv"))

# bar chart 
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

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(df, varname="PR", 
                    groupnames=c("condition"))

# 1. select placebo

pla <- df2[c(10,11,12,13,14),]

pla$condition <- factor(pla$condition,levels = c("pla_ket", "pla_kstudy_high", "pla_kstudy_low", "pla_lsd", "pla_psi"))

pdf("high_low_placebo.pdf", width=5,height=5)
ggplot(data=pla, aes(x=condition, y=PR, fill = condition)) +
  geom_bar(stat="identity", width = 0.5) +
  theme1 +
  geom_errorbar(aes(ymin=PR-sd, ymax=PR+sd), width=0, size = 1.5,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("pla_ket" = "#666666", "pla_kstudy_high" = "darkgrey", "pla_kstudy_low" = "lightgrey", "pla_lsd" = "#A99BBA", "pla_psi" = "#85719C"))
dev.off()

#t-test high low
pla_ttest <- df[df$condition == 'pla_kstudy_high' | df$condition == 'pla_kstudy_low',]
pla_ttest$condition <- factor(pla_ttest$condition,levels = c("pla_kstudy_high", "pla_kstudy_low"))

res <- t.test(PR ~ condition, data = pla_ttest, var.equal = TRUE)
res

data_summary(pla_ttest, varname="PR", groupnames=c("condition"))

# 1. select substance

sub <- df2[c(4,5,6,9,15),]

sub$condition <- factor(sub$condition,levels = c("ket", "ket_sub_high", "ket_sub_low", "lsd", "psi"))

pdf("high_low_subcebo.pdf", width=5,height=5)
ggplot(data=sub, aes(x=condition, y=PR, fill = condition)) +
  geom_bar(stat="identity", width = 0.5) +
  theme1 +
  geom_errorbar(aes(ymin=PR-sd, ymax=PR+sd), width=0, size = 1.5,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("ket" = "#666666", "ket_sub_high" = "darkgrey", "ket_sub_low" = "lightgrey", "lsd" = "#A99BBA", "psi" = "#85719C"))
dev.off()

#t-test high low
sub_ttest <- df[df$condition == 'ket_sub_high' | df$condition == 'ket_sub_low',]

res <- t.test(PR ~ condition, data = sub_ttest, var.equal = TRUE)
res

data_summary(sub_ttest, varname="PR", groupnames=c("condition"))




df3 <- df2[!(df2$condition == "pla"),] 

df3$condition <- factor(df3$condition,levels = c("lsd", "psi", "ket"))

pdf("pr_ket_psi_lsd_v2.pdf", width=3.5,height=5)
ggplot(data=df3, aes(x=condition, y=PR, fill = condition)) +
  geom_bar(stat="identity", width = 0.5) +
  theme1 +
  geom_errorbar(aes(ymin=PR-sd, ymax=PR+sd), width=0, size = 1.5,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("ket" = "#666666", "lsd" = "#A99BBA", "psi" = "#85719C"))
dev.off()

df4 <- df2[!(df2$condition == "pla"),] 

df4 <- df4[c(1,2,3,5,6),]

df4$condition <- factor(df4$condition,levels = c("delta_lsd", "delta_psi", "delta_ket", "ket23", "ket23_high"))

pdf("pr_ket23_ket_psi_lsd.pdf", width=5,height=5)
ggplot(data=df4, aes(x=condition, y=PR, fill = condition)) +
  geom_bar(stat="identity", width = 0.5) +
  theme1 +
  geom_errorbar(aes(ymin=PR-sd, ymax=PR+sd), width=0, size = 1.5,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("delta_ket" = "#666666", "ket23" = "darkgrey", "ket23_high" = "lightgrey", "delta_lsd" = "#A99BBA", "delta_psi" = "#85719C"))
dev.off()

#t-test high low
delta_ttest <- df[df$condition == 'ket23' | df$condition == 'ket23_high',]

res <- t.test(PR ~ condition, data = delta_ttest, var.equal = TRUE)
res

data_summary(delta_ttest, varname="PR", groupnames=c("condition"))




# read dataframes
df = data.frame(read.csv("participationratio_19_7_21.csv"))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(df, varname="PR", 
                    groupnames=c("condition"))

df3 <-df2[!(df2$condition=="ket" | df2$condition=="lsd" | df2$condition=="psi"),]

df3$condition <- factor(df3$condition,levels = c("pla_ket", "delta_ket", "pla_psi", "delta_psi", "pla_lsd", "delta_lsd"))

pdf("pr_PLACEBO_vs_delta.pdf", width=5,height=7)
ggplot(data=df3, aes(x=condition, y=PR, fill = condition)) +
  geom_bar(stat="identity", width = 0.5) +
  theme1 +
  geom_errorbar(aes(ymin=PR-sd, ymax=PR+sd), width=0, size = 1.5,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("pla_ket" = "#666666", "delta_ket" = "#666666", "pla_lsd" = "#A99BBA", "delta_lsd" = "#A99BBA", "pla_psi" = "#85719C", "delta_psi" = "#85719C"))
dev.off()

