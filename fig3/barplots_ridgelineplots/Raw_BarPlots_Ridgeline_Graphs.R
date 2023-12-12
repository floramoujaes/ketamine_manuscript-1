# Install Libraries:
install.packages("lme4")
install.packages("ez")
install.packages("nlme")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("aspace")
install.packages("PerformanceAnalytics")
install.packages("Hmisc")
install.packages("psych")
install.packages("lmerTest")
install.packages("effsize")
install.packages("reshape")
install.packages("reshape2")

library(effsize)
library(ez)
library(lme4)
library(nlme)
library(ggplot2)
library(dplyr) 
library(aspace) 
library(PerformanceAnalytics) 
library(car)
library(grid)
library(gridExtra)
library(reshape)
library(lmerTest)
library(ggridges)
library(ggplot2)
library(reshape)
library(reshape2)

# set working directory
setwd("~/Dropbox/DP5_N-Bridge/Figure_Calculations/GSR/N-Bridge/behavdelta_neurodelta_GSR_mVWMWB1d/behaviour_graphs")


# 1. plot bar chart of placebo, ketamine & delta raw scores
# 2. Ridgeline Z score graphs
# 3. Ridgeline PC1/2 graphs

# input files for 1 & 2 in /gpfs/loomis/pi/n3/Studies/Anticevic.DP5/analysis/n-bridge/Input_Files/behavior
# input files for 3 in /gpfs/loomis/pi/n3/Studies/Anticevic.DP5/analysis/n-bridge/behavdelta_neurodelta_GSR_mVWMWB1d/analysis/figures/NBRIDGE_behavdelta_neurodelta_GSR_mVWMWB1d_prep_CG_BehaviorPCARidgeline_PC1.tsv

# panss dataset
panss <- read.csv("R_DP5_PREPOSTKET.csv", header=TRUE)
# select the 40 participants (first 40 rows)
panss <- panss[1:40,]
panss[panss == 923] <- NA

# panss delta 
panss$p_delta <- (panss$panss_p_post - panss$panss_p_pre)
panss$n_delta <- (panss$panss_n_psot - panss$panss_n_pre)
panss$g_delta <- (panss$panss_g_post - panss$panss_g_pre)


# cognition dataset
cog <- read.csv("DP5_HCS_sWM_Angle_Means_N40.csv", header=TRUE)
# drop first column
cog <- cog[,2:4]
colnames(cog)[1] <- "Subject"
colnames(cog)[2] <- "Group"
colnames(cog)[3] <- "Angular_Distance"
# sign flip values
cog$Angular_Distance <- cog$Angular_Distance * -1


#-------------------------------------------------->
# 1. plot bar chart of placebo, ketamine & delta raw scores
#-------------------------------------------------->

# COGNITION

# plot ket & pla
# define column order
cog$Group <- factor(cog$Group,levels = c("HCS_Placebo","HCS_Ketamine"))

pdf("Bar_chart_raw_plaket_cognitive_flip15.pdf", width=1.5,height=3.5)
ggplot(cog, aes(cog, x=Group, y=Angular_Distance, fill=Group)) + 
  stat_summary(fun="mean", geom="bar", colour='white', size=1) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 1, color = "black") +
  theme_classic() + 
  scale_fill_manual(values=c("darkgrey","#666666")) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(-15, 0))
dev.off()

x <- cog$Angular_Distance[cog$Group == "HCS_Ketamine"]
y <- cog$Angular_Distance[cog$Group == "HCS_Placebo"]
t.test(x, y, paired = TRUE, alternative = "two.sided")
# t = -3.5475, df = 39, p-value = 0.001031

cohen.d(x, y, paired=TRUE, na.rm=FALSE, 
        hedges.correction = FALSE, conf.level = 0.95)
# d estimate: -0.6145254 (medium)

# plot delta 
# need to unmelt data to create the delta
cog_d <- dcast(data = cog,formula = Subject~Group,fun.aggregate = sum,value.var = "Angular_Distance")
cog_d$c_delta <- (cog_d$HCS_Ketamine - cog_d$HCS_Placebo)
cog_melt <- melt(cog_d, id=c("Subject"))
cog_melt_d <- cog_melt[cog_melt$variable == 'c_delta',]

pdf("Bar_chart_raw_delta_cognitive_flip15.pdf", width=1.5,height=3.7)
ggplot(cog_melt_d, aes(cog_melt_d, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", aes(width=0.7), position=position_dodge(width=0.95), colour='black', size=2) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  scale_fill_manual(values=c("#76B4C1")) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(-15, 0))
dev.off()


# PANSS POSITIVE

# panss pos ket pla
# panss pos subset data
colnames(panss)
panss_sub <- panss[,c("ID", "panss_p_pre", "panss_p_post")]
panss_sub_melt <- melt(panss_sub, id=c("ID"))
# define column order
panss_sub_melt$variable <- factor(panss_sub_melt$variable,levels = c("panss_p_pre","panss_p_post"))

pdf("Bar_chart_raw_ketpla_positive.pdf", width=1.5,height=3.5)
ggplot(panss_sub_melt, aes(panss_sub_melt, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", colour='white', size=1) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 1, color = "black") +
  theme_classic() + 
  scale_fill_manual(values=c("darkgrey","#666666")) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 15))
dev.off()

x <- panss_sub_melt$value[panss_sub_melt$variable == "panss_p_post"]
y <- panss_sub_melt$value[panss_sub_melt$variable == "panss_p_pre"]
t.test(x, y, paired = TRUE, alternative = "two.sided")
# t = 11.208, df = 35, p-value = 3.955e-13

cohen.d(x, y, paired=TRUE, na.rm=TRUE, 
        hedges.correction = FALSE, conf.level = 0.95)
# d estimate: 2.446756 (large)

# panss pos delta
panss_sub <- panss[,c("ID", "p_delta")]
panss_sub_melt <- melt(panss_sub, id=c("ID"))

pdf("Bar_chart_raw_delta_positive.pdf", width=1.5,height=3.7)
ggplot(panss_sub_melt, aes(panss_sub_melt, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", aes(width=0.7), position=position_dodge(width=0.95), colour='black', size=2) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  scale_fill_manual(values=c("#76B4C1")) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 15))
dev.off()

# PANSS NEGATIVE

# panss pos ket pla
# panss pos subset data
colnames(panss)
panss_sub <- panss[,c("ID", "panss_n_pre", "panss_n_psot")]
panss_sub_melt <- melt(panss_sub, id=c("ID"))
# define column order
panss_sub_melt$variable <- factor(panss_sub_melt$variable,levels = c("panss_n_pre","panss_n_psot"))

pdf("Bar_chart_raw_ketpla_negative.pdf", width=1.5,height=3.5)
ggplot(panss_sub_melt, aes(panss_sub_melt, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", colour='white', size=1) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 1, color = "black") +
  theme_classic() + 
  scale_fill_manual(values=c("darkgrey","#666666")) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 15))
dev.off()

x <- panss_sub_melt$value[panss_sub_melt$variable == "panss_n_psot"]
y <- panss_sub_melt$value[panss_sub_melt$variable == "panss_n_pre"]
t.test(x, y, paired = TRUE, alternative = "two.sided")
# t = 7.2319, df = 34, p-value = 2.276e-08

cohen.d(x, y, paired=TRUE, na.rm=TRUE, 
        hedges.correction = FALSE, conf.level = 0.95)
# d estimate: 1.658179 (large)


# panss pos delta
panss_sub <- panss[,c("ID", "n_delta")]
panss_sub_melt <- melt(panss_sub, id=c("ID"))

pdf("Bar_chart_raw_delta_negative.pdf", width=1.5,height=3.7)
ggplot(panss_sub_melt, aes(panss_sub_melt, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", aes(width=0.7), position=position_dodge(width=0.95), colour='black', size=2) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  scale_fill_manual(values=c("#76B4C1")) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 15))
dev.off()


# PANSS GENERAL

# panss pos ket pla
# panss pos subset data
colnames(panss)
panss_sub <- panss[,c("ID", "panss_g_pre", "panss_g_post")]
panss_sub_melt <- melt(panss_sub, id=c("ID"))
# define column order
panss_sub_melt$variable <- factor(panss_sub_melt$variable,levels = c("panss_g_pre","panss_g_post"))

pdf("Bar_chart_raw_ketpla_general.pdf", width=1.5,height=3.5)
ggplot(panss_sub_melt, aes(panss_sub_melt, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", colour='white', size=1) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 1, color = "black") +
  theme_classic() + 
  scale_fill_manual(values=c("darkgrey","#666666")) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 30))
dev.off()

x <- panss_sub_melt$value[panss_sub_melt$variable == "panss_g_post"]
y <- panss_sub_melt$value[panss_sub_melt$variable == "panss_g_pre"]
t.test(x, y, paired = TRUE, alternative = "two.sided")
# t = 10.133, df = 34, p-value = 8.307e-12

cohen.d(x, y, paired=TRUE, na.rm=TRUE, 
        hedges.correction = FALSE, conf.level = 0.95)
# d estimate: 2.585368 (large)

# panss pos delta
panss_sub <- panss[,c("ID", "g_delta")]
panss_sub_melt <- melt(panss_sub, id=c("ID"))

pdf("Bar_chart_raw_delta_general.pdf", width=1.5,height=3.7)
ggplot(panss_sub_melt, aes(panss_sub_melt, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", aes(width=0.7), position=position_dodge(width=0.95), colour='black', size=2) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  scale_fill_manual(values=c("#76B4C1")) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 30))
dev.off()



#------------------------------------------------------->
# 2. Ridgeline Z score graphs: treat ket, pla as one sample
#------------------------------------------------------->

# cognition dataset
cog <- read.csv("DP5_HCS_sWM_Angle_Means_N40.csv", header=TRUE)
# drop first column
cog <- cog[,2:4]
colnames(cog)[1] <- "Subject"
colnames(cog)[2] <- "Group"
colnames(cog)[3] <- "Angular_Distance"
# sign flip values
cog$Angular_Distance <- cog$Angular_Distance * -1

# COGNITION Z Score
# delta
cog_d <- dcast(data = cog,formula = Subject~Group,fun.aggregate = sum,value.var = "Angular_Distance")

# melt
cog_m <- melt(cog_d, id=c("Subject"))

# z score ket & pla values
cog_m$z_score <- scale(cog_m$value)

# factor
cog_m$variable <- factor(cog_m$variable,levels = c("HCS_Placebo","HCS_Ketamine"))

pdf("Ridgeline_Z_ketpla_cognitive_flip50_v2.pdf", width=5.6,height=3.5)
ggplot(cog_m, aes(z_score, variable, group = variable)) + geom_density_ridges(aes(fill=variable), lwd=1.5, rel_min_height = 0.01) +
  # add colours
  scale_fill_manual(values = c("#666666","darkgray")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) +
  coord_cartesian(xlim = c(-6, 4))
dev.off()

x <- cog_m$z_score[cog_m$variable == "HCS_Placebo"]
y <- cog_m$z_score[cog_m$variable == "HCS_Ketamine"]
Var.test(x, y,paired=TRUE)
# t = -5.3622, df = 38, p-value = 4.263e-06

cohen.d(x, y, paired=TRUE, na.rm=TRUE, 
        hedges.correction = FALSE, conf.level = 0.95)
# d estimate: 0.6145254 (medium)

# delta
ket <- cog_m[cog_m$variable == 'HCS_Ketamine',]
colnames(ket)[3] <- "ket_value"
colnames(ket)[4] <- "ket_zscore"

pla <- cog_m[cog_m$variable == 'HCS_Placebo',]
colnames(pla)[3] <- "pla_value"
colnames(pla)[4] <- "pla_zscore"

# combine the dataframes
delta <- cbind(ket, pla)

# subtract ket z - pla z
delta$diff_zscore <- delta$ket_zscore - delta$pla_zscore

# subset
myvars <- c("Subject", "pla_zscore", "ket_zscore", "diff_zscore")
delta_s <- delta[myvars]
write.csv(delta_s, "delta_s.csv")


# melt by subject
df2 <- read.csv("delta_s.csv", header=TRUE)
df3 <- subset(df2, select = -c(X))
delta_m2 <- melt(df3, id=c("Subject"))

# select delta
delta_m2$variable <- factor(delta_m2$variable,levels = c("diff_zscore","ket_zscore","pla_zscore"))

pdf("Ridgeline_Z_all_cognitive_flip50_V2.pdf", width=5.6,height=3.5)
ggplot(delta_m2, aes(value, variable, group = variable)) + geom_density_ridges(aes(fill=variable), lwd=2, rel_min_height = 0.01, scale = 2) +
  # add colours
  scale_fill_manual(values = c("#76B4C1","#008299","#978d85")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) +
  coord_cartesian(xlim = c(-6, 4))
dev.off()

Z_ketpla <- delta_m2[delta_m2$variable == 'ket_zscore' | delta_m2$variable == 'pla_zscore',]
pdf("Ridgeline_Z_ketpla_cognitive_flip50_V2.pdf", width=5.6,height=3.5)
ggplot(Z_ketpla, aes(value, variable, group = variable)) + geom_density_ridges(aes(fill=variable), lwd=2, rel_min_height = 0.01, scale = 2) +
  # add colours
  scale_fill_manual(values = c("gray", "white")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) +
  coord_cartesian(xlim = c(-6, 4))
dev.off()


Z_ketpla <- delta_m2[delta_m2$variable == 'diff_zscore',]
pdf("Ridgeline_Z_delta_cognitive_flip50_V2.pdf", width=5.6,height=2.2)
ggplot(Z_ketpla, aes(value, variable, group = variable)) + geom_density_ridges(aes(fill=variable), lwd=2, rel_min_height = 0.01, scale = 2) +
  # add colours
  scale_fill_manual(values = c("#76B4C1")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) +
  coord_cartesian(xlim = c(-6, 4))
dev.off()



# PANSS POS Z Score

# ** clear workspace ***

# panss dataset
panss <- read.csv("R_DP5_PREPOSTKET.csv", header=TRUE)
# select the 40 participants (first 40 rows)
panss <- panss[1:40,]
panss[panss == 923] <- NA

# z score
panss_zp <- panss[,c("ID", "panss_p_pre", "panss_p_post")]

# melt
panss_zp_m <- melt(panss_zp, id=c("ID"))

# z score ket & pla values
panss_zp_m$z_score <- scale(panss_zp_m$value)

# delta
ket <- panss_zp_m[panss_zp_m$variable == 'panss_p_post',]
colnames(ket)[3] <- "ket_value"
colnames(ket)[4] <- "ket_zscore"

pla <- panss_zp_m[panss_zp_m$variable == 'panss_p_pre',]
colnames(pla)[3] <- "pla_value"
colnames(pla)[4] <- "pla_zscore"

# combine the dataframes
delta <- cbind(ket, pla)

# subtract ket z - pla z
delta$diff_zscore <- delta$ket_zscore - delta$pla_zscore

# subset
myvars <- c("ID", "pla_zscore", "ket_zscore", "diff_zscore")
delta_s <- delta[myvars]
write.csv(delta_s, "delta_s.csv")


# melt by subject
df2 <- read.csv("delta_s.csv", header=TRUE)
df3 <- subset(df2, select = -c(X))
delta_m2 <- melt(df3, id=c("ID"))

# select delta
delta_m2$variable <- factor(delta_m2$variable,levels = c("diff_zscore","ket_zscore","pla_zscore"))

pdf("Ridgeline_Z_all_positive_V2.pdf", width=5.6,height=3.5)
ggplot(delta_m2, aes(value, variable, group = variable)) + geom_density_ridges(aes(fill=variable), lwd=2, rel_min_height = 0.01, scale = 2) +
  # add colours
  scale_fill_manual(values = c("#76B4C1","#008299","#978d85")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) +
  coord_cartesian(xlim = c(-6, 4))
dev.off()

Z_ketpla <- delta_m2[delta_m2$variable == 'ket_zscore' | delta_m2$variable == 'pla_zscore',]
pdf("Ridgeline_Z_ketpla_positive_flip50_V2.pdf", width=5.6,height=3.5)
ggplot(Z_ketpla, aes(value, variable, group = variable)) + geom_density_ridges(aes(fill=variable), lwd=1.5, rel_min_height = 0.01, scale = 2) +
  # add colours
  scale_fill_manual(values = c("#666666","darkgray")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) +
  coord_cartesian(xlim = c(-6, 4))
dev.off()

library('PairedData') # pitman morgan test
x <- delta_m2$value[delta_m2$variable == "ket_zscore"]
y <- delta_m2$value[delta_m2$variable == "pla_zscore"]
Var.test(x, y,paired=TRUE)
# t = 11.655, df = 34, p-value = 2.007e-13

cohen.d(x, y, paired=TRUE, na.rm=TRUE, 
        hedges.correction = FALSE, conf.level = 0.95)
# d estimate: 0.6145254 (medium)

Z_ketpla <- delta_m2[delta_m2$variable == 'diff_zscore',]
pdf("Ridgeline_Z_delta_positive_flip50_V2.pdf", width=5.6,height=2.2)
ggplot(Z_ketpla, aes(value, variable, group = variable)) + geom_density_ridges(aes(fill=variable), lwd=2, rel_min_height = 0.01, scale = 2) +
  # add colours
  scale_fill_manual(values = c("#76B4C1")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) +
  coord_cartesian(xlim = c(-6, 4))
dev.off()



# PANSS NEG Z Score

# ** clear workspace ***

# panss dataset
panss <- read.csv("R_DP5_PREPOSTKET.csv", header=TRUE)
# select the 40 participants (first 40 rows)
panss <- panss[1:40,]
panss[panss == 923] <- NA

# z score
panss_zn <- panss[,c("ID", "panss_n_pre", "panss_n_psot")]

# melt
panss_zn_m <- melt(panss_zn, id=c("ID"))

# z score ket & pla values
panss_zn_m$z_score <- scale(panss_zn_m$value)

# delta
ket <- panss_zn_m[panss_zn_m$variable == 'panss_n_psot',]
colnames(ket)[3] <- "ket_value"
colnames(ket)[4] <- "ket_zscore"

pla <- panss_zn_m[panss_zn_m$variable == 'panss_n_pre',]
colnames(pla)[3] <- "pla_value"
colnames(pla)[4] <- "pla_zscore"

# combine the dataframes
delta <- cbind(ket, pla)

# subtract ket z - pla z
delta$diff_zscore <- delta$ket_zscore - delta$pla_zscore

# subset
myvars <- c("ID", "pla_zscore", "ket_zscore", "diff_zscore")
delta_s <- delta[myvars]
write.csv(delta_s, "delta_s.csv")


# melt by subject
df2 <- read.csv("delta_s.csv", header=TRUE)
df3 <- subset(df2, select = -c(X))
delta_m2 <- melt(df3, id=c("ID"))

# select delta
delta_m2$variable <- factor(delta_m2$variable,levels = c("diff_zscore","ket_zscore","pla_zscore"))

pdf("Ridgeline_Z_all_negative_V2.pdf", width=5.6,height=3.5)
ggplot(delta_m2, aes(value, variable, group = variable)) + geom_density_ridges(aes(fill=variable), lwd=2, rel_min_height = 0.01, scale = 2) +
  # add colours
  scale_fill_manual(values = c("#76B4C1","#008299","#978d85")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) +
  coord_cartesian(xlim = c(-6, 4))
dev.off()

Z_ketpla <- delta_m2[delta_m2$variable == 'ket_zscore' | delta_m2$variable == 'pla_zscore',]
pdf("Ridgeline_Z_ketpla_negative_flip50_V2.pdf", width=5.6,height=3.5)
ggplot(Z_ketpla, aes(value, variable, group = variable)) + geom_density_ridges(aes(fill=variable), lwd=1.5, rel_min_height = 0.01, scale = 2) +
  # add colours
  # add colours
  scale_fill_manual(values = c("#666666","darkgray")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) +
  coord_cartesian(xlim = c(-6, 4))
dev.off()


x <- delta_m2$value[delta_m2$variable == "ket_zscore"]
y <- delta_m2$value[delta_m2$variable == "pla_zscore"]
Var.test(x, y,paired=TRUE)
# t = 14.858, df = 33, p-value = 4.441e-16


cohen.d(x, y, paired=TRUE, na.rm=TRUE, 
        hedges.correction = FALSE, conf.level = 0.95)
# d estimate: 1.658179 (large)

Z_ketpla <- delta_m2[delta_m2$variable == 'diff_zscore',]
pdf("Ridgeline_Z_delta_negative_flip50_V2.pdf", width=5.6,height=2.2)
ggplot(Z_ketpla, aes(value, variable, group = variable)) + geom_density_ridges(aes(fill=variable), lwd=2, rel_min_height = 0.01, scale = 2) +
  # add colours
  scale_fill_manual(values = c("#76B4C1")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) +
  coord_cartesian(xlim = c(-6, 4))
dev.off()


# PANSS GENERAL Z Score

# ** clear workspace ***

# panss dataset
panss <- read.csv("R_DP5_PREPOSTKET.csv", header=TRUE)
# select the 40 participants (first 40 rows)
panss <- panss[1:40,]
panss[panss == 923] <- NA

# z score
panss_zg <- panss[,c("ID", "panss_g_pre", "panss_g_post")]

# melt
panss_zg_m <- melt(panss_zg, id=c("ID"))

# z score ket & pla values
panss_zg_m$z_score <- scale(panss_zg_m$value)

# delta
ket <- panss_zg_m[panss_zg_m$variable == 'panss_g_post',]
colnames(ket)[3] <- "ket_value"
colnames(ket)[4] <- "ket_zscore"

pla <- panss_zg_m[panss_zg_m$variable == 'panss_g_pre',]
colnames(pla)[3] <- "pla_value"
colnames(pla)[4] <- "pla_zscore"

# combine the dataframes
delta <- cbind(ket, pla)

# subtract ket z - pla z
delta$diff_zscore <- delta$ket_zscore - delta$pla_zscore

# subset
myvars <- c("ID", "pla_zscore", "ket_zscore", "diff_zscore")
delta_s <- delta[myvars]
write.csv(delta_s, "delta_s.csv")


# melt by subject
df2 <- read.csv("delta_s.csv", header=TRUE)
df3 <- subset(df2, select = -c(X))
delta_m2 <- melt(df3, id=c("ID"))

# select delta
delta_m2$variable <- factor(delta_m2$variable,levels = c("diff_zscore","ket_zscore","pla_zscore"))

pdf("Ridgeline_Z_all_general_V2.pdf", width=5.6,height=3.5)
ggplot(delta_m2, aes(value, variable, group = variable)) + geom_density_ridges(aes(fill=variable), lwd=2, rel_min_height = 0.01, scale = 2) +
  # add colours
  scale_fill_manual(values = c("#76B4C1","#008299","#978d85")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) +
  coord_cartesian(xlim = c(-6, 4))
dev.off()

Z_ketpla <- delta_m2[delta_m2$variable == 'ket_zscore' | delta_m2$variable == 'pla_zscore',]
pdf("Ridgeline_Z_ketpla_general_flip50_V2.pdf", width=5.6,height=3.5)
ggplot(Z_ketpla, aes(value, variable, group = variable)) + geom_density_ridges(aes(fill=variable), lwd=1.5, rel_min_height = 0.01, scale = 2) +
  # add colours
  # add colours
  scale_fill_manual(values = c("#666666","darkgray")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) +
  coord_cartesian(xlim = c(-6, 4))
dev.off()

x <- delta_m2$value[delta_m2$variable == "ket_zscore"]
y <- delta_m2$value[delta_m2$variable == "pla_zscore"]
Var.test(x, y,paired=TRUE)
# t = 10.34, df = 33, p-value = 6.968e-12

cohen.d(x, y, paired=TRUE, na.rm=TRUE, 
        hedges.correction = FALSE, conf.level = 0.95)
# d estimate: 2.585368 (large)


Z_ketpla <- delta_m2[delta_m2$variable == 'diff_zscore',]
pdf("Ridgeline_Z_delta_general_flip50_V2.pdf", width=5.6,height=2.2)
ggplot(Z_ketpla, aes(value, variable, group = variable)) + geom_density_ridges(aes(fill=variable), lwd=2, rel_min_height = 0.01, scale = 2) +
  # add colours
  scale_fill_manual(values = c("#76B4C1")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) +
  coord_cartesian(xlim = c(-6, 4))
dev.off()



#------------------------------------------------------->
# 3. Ridgeline PC1/2 graphs
#------------------------------------------------------->


# this is labelled pc1 but it actually has all the pcs in it:
pc <- read.table(file = 'NBRIDGE_behavdelta_neurodelta_GSR_mVWMWB1d_prep_CG_BehaviorPCARidgeline_PC1.tsv', sep = '\t', header = FALSE)

# subset only delta not controls
delta <- subset(pc, V1 == "CON-D")

# select:
# V2 for PC1
# V3 for PC2

pdf("Ridgeline_dist_delta_pc1.pdf", width=9,height=3.5)

# change last line x values to make the graph symmetrical 
ggplot(delta, aes(V2, V1, group = V1)) + geom_density_ridges(aes(fill=V1), lwd=2, scale=10) + # change scale to increase height
  # add colour
  scale_fill_manual(values = c("CON-D"="#76B4C1")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("PC1") +
  xlab("Normalised PC Score (Z)") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE)  + 
  scale_x_continuous(limits=c(-50, 50))

dev.off()


pdf("Ridgeline_dist_delta_pc2.pdf", width=9,height=3.5)

# change last line x values to make the graph symmetrical 
ggplot(delta, aes(V3, V1, group = V1)) + geom_density_ridges(aes(fill=V1), lwd=2, scale=10) + # change scale to increase height
  # add colour
  scale_fill_manual(values = c("CON-D"="#76B4C1")) +
  # make y axis start at 0
  # scale_y_discrete(expand = c(0.01, 0)) +
  # adjust x axis to start/finish and add labels
  # scale_x_continuous(expand = c(0.01, 0), limits=c(-80,90), breaks=c(-80,-60,-40,-20,0,20,40,60, 80),labels=c(-80,"",-40,"",0,"",40,"",80)) +
  ylab("PC2") +
  xlab("Normalised PC Score (Z)") +
  #geom_vline(xintercept=0, lty=2,col="grey",lwd=1) +
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  #geom_segment(data=filter(df_ketpla, V1=="Ketamine"), aes(x = 20.414365, y =.9999, xend = 20.414365, yend = 2.7), colour="grey", lwd=2) +
  #geom_segment(data=filter(df_ketpla, V1=="CON"), aes(x = 17.458831, y =2.85, xend = 17.458831, yend = 3.3), colour="grey", lwd=2) +
  theme(axis.text.x = element_text(vjust=0.5, colour="black", size=12), axis.text.y = element_text(colour="white", size=12), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line.y = element_line(color = 'gray', size = 2), axis.line.x = element_line(color = 'gray', size = 2)) + 
  #remove legend position
  guides(fill=FALSE) + 
  scale_x_continuous(limits=c(-50, 50))

dev.off()


