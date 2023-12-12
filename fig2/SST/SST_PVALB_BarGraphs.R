# install
install.packages("prevalence")
# Library
library("data.table")
library('prevalence')
library('Hmisc')
library(tidyverse)


# set working directory 
setwd("~/Dropbox/DP5_N-Bridge/Figure_Calculations/GSR/genetics_graph/gene_images_network")

# read in dataframe
df <- read.csv("gene_input.csv", sep = ",", header = TRUE)

# list column headers
colnames(df)

#---------------------------------------->
# "GABRA3"   

# mask image
df$GABRA3_sens <- df$Sensory * df$GABRA3
df$GABRA3_assoc <- df$Association * df$GABRA3
df$GABRA3_assoc_pmmdmm <- df$Association_PMM_VMM * df$GABRA3

# select columns for delta bar graph 1, 36-52
df_GABRA3 <- df %>% select(Parcels, GABRA3_sens, GABRA3_assoc)
df_GABRA3_pmmdmm <- df %>% select(Parcels, GABRA3_sens, GABRA3_assoc_pmmdmm)

# melt
df_GABRA3 <- melt(df_GABRA3, id=c("Parcels"))
df_GABRA3_pmmdmm <- melt(df_GABRA3_pmmdmm, id=c("Parcels"))

# define order 
df_GABRA3$variable <- factor(df_GABRA3$variable,levels = c(
  "GABRA3_assoc", #ED3833
  "GABRA3_sens")) #7847EB

df_GABRA3_pmmdmm$variable <- factor(df_GABRA3_pmmdmm$variable,levels = c(
  "GABRA3_assoc_pmmdmm", #ED3833
  "GABRA3_sens")) #7847EB

pdf("GABRA3_assoc_sens.pdf", width=2,height=4.5)
ggplot(df_GABRA3, aes(df_GABRA3, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()

pdf("GABRA3_assocpmmdmm_sens.pdf", width=2,height=4.5)
ggplot(df_GABRA3_pmmdmm, aes(df_GABRA3_pmmdmm, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()

#---------------------------------------->
# "GABRA4"    

# mask image
df$GABRA4_sens <- df$Sensory * df$GABRA4
df$GABRA4_assoc <- df$Association * df$GABRA4
df$GABRA4_assoc_pmmdmm <- df$Association_PMM_VMM * df$GABRA4

# select columns for delta bar graph 1, 36-52
df_GABRA4 <- df %>% select(Parcels, GABRA4_sens, GABRA4_assoc)
df_GABRA4_pmmdmm <- df %>% select(Parcels, GABRA4_sens, GABRA4_assoc_pmmdmm)

# melt
df_GABRA4 <- melt(df_GABRA4, id=c("Parcels"))
df_GABRA4_pmmdmm <- melt(df_GABRA4_pmmdmm, id=c("Parcels"))

# define order 
df_GABRA4$variable <- factor(df_GABRA4$variable,levels = c(
  "GABRA4_assoc", #ED3833
  "GABRA4_sens")) #7847EB

df_GABRA4_pmmdmm$variable <- factor(df_GABRA4_pmmdmm$variable,levels = c(
  "GABRA4_assoc_pmmdmm", #ED3833
  "GABRA4_sens")) #7847EB

pdf("GABRA4_assoc_sens.pdf", width=2,height=4.5)
ggplot(df_GABRA4, aes(df_GABRA4, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()

pdf("GABRA4_assocpmmdmm_sens.pdf", width=2,height=4.5)
ggplot(df_GABRA4_pmmdmm, aes(df_GABRA4_pmmdmm, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()
#---------------------------------------->
# "GABRA5"  

# mask image
df$GABRA5_sens <- df$Sensory * df$GABRA5
df$GABRA5_assoc <- df$Association * df$GABRA5
df$GABRA5_assoc_pmmdmm <- df$Association_PMM_VMM * df$GABRA5

# select columns for delta bar graph 1, 36-52
df_GABRA5 <- df %>% select(Parcels, GABRA5_sens, GABRA5_assoc)
df_GABRA5_pmmdmm <- df %>% select(Parcels, GABRA5_sens, GABRA5_assoc_pmmdmm)

# melt
df_GABRA5 <- melt(df_GABRA5, id=c("Parcels"))
df_GABRA5_pmmdmm <- melt(df_GABRA5_pmmdmm, id=c("Parcels"))

# define order 
df_GABRA5$variable <- factor(df_GABRA5$variable,levels = c(
  "GABRA5_assoc", #ED3833
  "GABRA5_sens")) #7847EB

df_GABRA5_pmmdmm$variable <- factor(df_GABRA5_pmmdmm$variable,levels = c(
  "GABRA5_assoc_pmmdmm", #ED3833
  "GABRA5_sens")) #7847EB

pdf("GABRA5_assoc_sens.pdf", width=2,height=4.5)
ggplot(df_GABRA5, aes(df_GABRA5, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()

pdf("GABRA5_assocpmmdmm_sens.pdf", width=2,height=4.5)
ggplot(df_GABRA5_pmmdmm, aes(df_GABRA5_pmmdmm, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()
#---------------------------------------->
# "GABRA6"  

# mask image
df$GABRA6_sens <- df$Sensory * df$GABRA6
df$GABRA6_assoc <- df$Association * df$GABRA6
df$GABRA6_assoc_pmmdmm <- df$Association_PMM_VMM * df$GABRA6

# select columns for delta bar graph 1, 36-52
df_GABRA6 <- df %>% select(Parcels, GABRA6_sens, GABRA6_assoc)
df_GABRA6_pmmdmm <- df %>% select(Parcels, GABRA6_sens, GABRA6_assoc_pmmdmm)

# melt
df_GABRA6 <- melt(df_GABRA6, id=c("Parcels"))
df_GABRA6_pmmdmm <- melt(df_GABRA6_pmmdmm, id=c("Parcels"))

# define order 
df_GABRA6$variable <- factor(df_GABRA6$variable,levels = c(
  "GABRA6_assoc", #ED3833
  "GABRA6_sens")) #7847EB

df_GABRA6_pmmdmm$variable <- factor(df_GABRA6_pmmdmm$variable,levels = c(
  "GABRA6_assoc_pmmdmm", #ED3833
  "GABRA6_sens")) #7847EB

pdf("GABRA6_assoc_sens.pdf", width=2,height=4.5)
ggplot(df_GABRA6, aes(df_GABRA6, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()

pdf("GABRA6_assocpmmdmm_sens.pdf", width=2,height=4.5)
ggplot(df_GABRA6_pmmdmm, aes(df_GABRA6_pmmdmm, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()
#---------------------------------------->
# "GABRB1" 

# mask image
df$GABRB1_sens <- df$Sensory * df$GABRB1
df$GABRB1_assoc <- df$Association * df$GABRB1
df$GABRB1_assoc_pmmdmm <- df$Association_PMM_VMM * df$GABRB1

# select columns for delta bar graph 1, 36-52
df_GABRB1 <- df %>% select(Parcels, GABRB1_sens, GABRB1_assoc)
df_GABRB1_pmmdmm <- df %>% select(Parcels, GABRB1_sens, GABRB1_assoc_pmmdmm)

# melt
df_GABRB1 <- melt(df_GABRB1, id=c("Parcels"))
df_GABRB1_pmmdmm <- melt(df_GABRB1_pmmdmm, id=c("Parcels"))

# define order 
df_GABRB1$variable <- factor(df_GABRB1$variable,levels = c(
  "GABRB1_assoc", #ED3833
  "GABRB1_sens")) #7847EB

df_GABRB1_pmmdmm$variable <- factor(df_GABRB1_pmmdmm$variable,levels = c(
  "GABRB1_assoc_pmmdmm", #ED3833
  "GABRB1_sens")) #7847EB

pdf("GABRB1_assoc_sens.pdf", width=2,height=4.5)
ggplot(df_GABRB1, aes(df_GABRB1, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()

pdf("GABRB1_assocpmmdmm_sens.pdf", width=2,height=4.5)
ggplot(df_GABRB1_pmmdmm, aes(df_GABRB1_pmmdmm, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()
#---------------------------------------->
# "GABRG3"  

# mask image
df$GABRG3_sens <- df$Sensory * df$GABRG3
df$GABRG3_assoc <- df$Association * df$GABRG3
df$GABRG3_assoc_pmmdmm <- df$Association_PMM_VMM * df$GABRG3

# select columns for delta bar graph 1, 36-52
df_GABRG3 <- df %>% select(Parcels, GABRG3_sens, GABRG3_assoc)
df_GABRG3_pmmdmm <- df %>% select(Parcels, GABRG3_sens, GABRG3_assoc_pmmdmm)

# melt
df_GABRG3 <- melt(df_GABRG3, id=c("Parcels"))
df_GABRG3_pmmdmm <- melt(df_GABRG3_pmmdmm, id=c("Parcels"))

# define order 
df_GABRG3$variable <- factor(df_GABRG3$variable,levels = c(
  "GABRG3_assoc", #ED3833
  "GABRG3_sens")) #7847EB

df_GABRG3_pmmdmm$variable <- factor(df_GABRG3_pmmdmm$variable,levels = c(
  "GABRG3_assoc_pmmdmm", #ED3833
  "GABRG3_sens")) #7847EB

pdf("GABRG3_assoc_sens.pdf", width=2,height=4.5)
ggplot(df_GABRG3, aes(df_GABRG3, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()

pdf("GABRG3_assocpmmdmm_sens.pdf", width=2,height=4.5)
ggplot(df_GABRG3_pmmdmm, aes(df_GABRG3_pmmdmm, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()
#---------------------------------------->
# "SST"    

# mask image
df$SST_sens <- df$Sensory * df$SST
df$SST_assoc <- df$Association * df$SST
df$SST_assoc_pmmdmm <- df$Association_PMM_VMM * df$SST

# select columns for delta bar graph 1, 36-52
df_SST <- df %>% select(Parcels, SST_sens, SST_assoc)
df_SST_pmmdmm <- df %>% select(Parcels, SST_sens, SST_assoc_pmmdmm)

# melt
df_SST <- melt(df_SST, id=c("Parcels"))
df_SST_pmmdmm <- melt(df_SST_pmmdmm, id=c("Parcels"))

# define order 
df_SST$variable <- factor(df_SST$variable,levels = c(
  "SST_assoc", #ED3833
  "SST_sens")) #7847EB

df_SST_pmmdmm$variable <- factor(df_SST_pmmdmm$variable,levels = c(
  "SST_assoc_pmmdmm", #ED3833
  "SST_sens")) #7847EB

pdf("SST_assoc_sens.pdf", width=2,height=4.5)
ggplot(df_SST, aes(df_SST, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()

pdf("SST_assocpmmdmm_sens.pdf", width=1.6,height=4.5)
ggplot(df_SST_pmmdmm, aes(df_SST_pmmdmm, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.2, 0.2)) 
dev.off()

#---------------------------------------->
# "GRIN3A" 

# mask image
df$GRIN3A_sens <- df$Sensory * df$GRIN3A
df$GRIN3A_assoc <- df$Association * df$GRIN3A
df$GRIN3A_assoc_pmmdmm <- df$Association_PMM_VMM * df$GRIN3A

# select columns for delta bar graph 1, 36-52
df_GRIN3A <- df %>% select(Parcels, GRIN3A_sens, GRIN3A_assoc)
df_GRIN3A_pmmdmm <- df %>% select(Parcels, GRIN3A_sens, GRIN3A_assoc_pmmdmm)

# melt
df_GRIN3A <- melt(df_GRIN3A, id=c("Parcels"))
df_GRIN3A_pmmdmm <- melt(df_GRIN3A_pmmdmm, id=c("Parcels"))

# define order 
df_GRIN3A$variable <- factor(df_GRIN3A$variable,levels = c(
  "GRIN3A_assoc", #ED3833
  "GRIN3A_sens")) #7847EB

df_GRIN3A_pmmdmm$variable <- factor(df_GRIN3A_pmmdmm$variable,levels = c(
  "GRIN3A_assoc_pmmdmm", #ED3833
  "GRIN3A_sens")) #7847EB

pdf("GRIN3A_assoc_sens.pdf", width=2,height=4.5)
ggplot(df_GRIN3A, aes(df_GRIN3A, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
 coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()

pdf("GRIN3A_assocpmmdmm_sens.pdf", width=2,height=4.5)
ggplot(df_GRIN3A_pmmdmm, aes(df_GRIN3A_pmmdmm, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
 coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()



#---------------------------------------->
# "PVALB" 

# mask image
df$PVALB_sens <- df$Sensory * df$PVALB
df$PVALB_assoc <- df$Association * df$PVALB
df$PVALB_assoc_pmmdmm <- df$Association_PMM_VMM * df$PVALB

# select columns for delta bar graph 1, 36-52
df_PVALB <- df %>% select(Parcels, PVALB_sens, PVALB_assoc)
df_PVALB_pmmdmm <- df %>% select(Parcels, PVALB_sens, PVALB_assoc_pmmdmm)

# melt
df_PVALB <- melt(df_PVALB, id=c("Parcels"))
df_PVALB_pmmdmm <- melt(df_PVALB_pmmdmm, id=c("Parcels"))

# define order 
df_PVALB$variable <- factor(df_PVALB$variable,levels = c(
  "PVALB_assoc", #ED3833
  "PVALB_sens")) #7847EB

df_PVALB_pmmdmm$variable <- factor(df_PVALB_pmmdmm$variable,levels = c(
  "PVALB_assoc_pmmdmm", #ED3833
  "PVALB_sens")) #7847EB

pdf("PVALB_assoc_sens.pdf", width=2,height=4.5)
ggplot(df_PVALB, aes(df_PVALB, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()

pdf("PVALB_assocpmmdmm_sens.pdf", width=1.6,height=4.5)
ggplot(df_PVALB_pmmdmm, aes(df_PVALB_pmmdmm, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.2, 0.2)) 
dev.off()



# # t-test
# # sign flip to match mean
# df$PC1_Association_flip <- df$PC1_Association * -1
# df$PC1_Sensory_flip <- df$PC1_Sensory * -1
# df$PC2_Association_flip <- df$PC2_Association * -1
# df$PC2_Sensory_flip <- df$PC2_Sensory * -1
# 
# t.test(df$PC1_Association_flip, df$PC1_Sensory_flip)
# t.test(df$PC2_Association_flip, df$PC2_Sensory_flip)

