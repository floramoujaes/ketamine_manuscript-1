# install
# install.packages("prevalence")
# Library
library("data.table")
library('prevalence')
library('Hmisc')
library('dplyr')

# set working directory 
setwd("~/Dropbox/DP5_N-Bridge/Figure_Calculations/GSR/N40_Neural_PCA_mVWMWB1d/association_sensory")


# PC1

# read in df
df <- read.table("PC1_ket-pla_pca_image_ptseries_ZSCORE.txt", header = FALSE)

# read in mask
mask <- read.table("network_mask.txt", header = TRUE)

# mask dataframe
df_masked <- mask * df$V1

# association/sensory graph

df2 <- df_masked[,c(1,7:18)]
colnames(df2)

df2$assoc <- rowSums(df2[ , c(5,6,7,8,10,11,12,13)],  na.rm=TRUE)
df2$sens <- rowSums(df2[ , c(2,3,4,9)], na.rm=TRUE)

t.test(df2$assoc, df2$sens)
# t = -11.393, df = 1368.2, p-value < 2.2e-16

# select columns for delta bar graph 1, 36-52
df_net <- df2 %>% select(Parcels, assoc, sens)

# melt
df_net <- melt(df_net, id=c("Parcels"))

# define order 
df_net$variable <- factor(df_net$variable,levels = c(
  "assoc", #ED3833
  "sens")) #7847EB

pdf("PC1_assoc_sens.pdf", width=2,height=4.5)
ggplot(df_net, aes(df_net, x=variable, y=value, fill=variable, width=0.65)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='black', size=2) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#000000", "#FFFFFF")) +
  coord_cartesian(ylim = c(-0.25, 0.2)) 
dev.off()


# clear
rm(list = ls())

# PC2

# read in df
df <- read.table("PC2_ket-pla_pca_image_ptseries_ZSCORE.txt", header = FALSE)

# read in mask
mask <- read.table("network_mask.txt", header = TRUE)

# mask dataframe
df_masked <- mask * df$V1

# association/sensory graph

df2 <- df_masked[,c(1,7:18)]
colnames(df2)

df2$assoc <- rowSums(df2[ , c(5,6,7,8,10,11,12,13)],  na.rm=TRUE)
df2$sens <- rowSums(df2[ , c(2,3,4,9)], na.rm=TRUE)

t.test(df2$assoc, df2$sens)
# t = 5.8858, df = 1147.3, p-value = 5.195e-09

# select columns for delta bar graph 1, 36-52
df_net <- df2 %>% select(Parcels, assoc, sens)

# melt
df_net <- melt(df_net, id=c("Parcels"))

# define order 
df_net$variable <- factor(df_net$variable,levels = c(
  "assoc", #ED3833
  "sens")) #7847EB

pdf("PC2_assoc_sens.pdf", width=2,height=4.5)
ggplot(df_net, aes(df_net, x=variable, y=value, fill=variable, width=0.65)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='black', size=2) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#000000", "#FFFFFF")) +
  coord_cartesian(ylim = c(-0.25, 0.2)) 
dev.off()

