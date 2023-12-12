# install
# install.packages("prevalence")
# Library
library("data.table")
library('prevalence')
library('Hmisc')
library('dplyr')

# set working directory 
setwd("~/Dropbox/DP5_N-Bridge/Figure_Calculations/GSR/Behaviour/cognition")

# read in mask
mask <- read.table("network_mask.txt", header = TRUE)

# read in df
df <- read.table("cognition_5000_dat_ztstat_pscalar.txt", header = FALSE)

# mask dataframe
df_masked <- mask * df$V1

# 1. network graphs
colnames(df_masked)

df_cortical <- df_masked[,c(1,7:18)]
df_subcortical <- df_masked[,c(1,2,3,4,5,6)]

# melt
df_cortical <- melt(df_cortical, id=c("Parcels"))
df_subcortical <- melt(df_subcortical, id=c("Parcels"))


# define order
#df_PC1$variable <- factor(df_PC1$variable,levels = c("DeltaPC1_Orbito.Affective", "DeltaPC1_Ventral.Multimodal", "DeltaPC1_Posterior.Multimodal", "DeltaPC1_Default", "DeltaPC1_Auditory", "DeltaPC1_Frontoparietal", "DeltaPC1_Language", "DeltaPC1_Dorsal.attention", "DeltaPC1_Cingulo.Opercular", "DeltaPC1_Somatomotor", "DeltaPC1_Visual2", "DeltaPC1_Visual1"))
# dfsubcortical_PC1$variable <- factor(dfsubcortical_PC1$variable,levels = c("DeltaPC1_Thalamus",  "DeltaPC1_Amygdala", "DeltaPC1_Striatum", "DeltaPC1_Hippocampus", "DeltaPC1_Cerebellum"))

# define order according to PC1 value
df_cortical$variable <- factor(df_cortical$variable,levels = c(
  "Default", #ED3833
  "Frontoparietal", #FFF851
  "Language", #3C9093
  "OrbitoAffective", #5C9A27
  "VentralMultimodal", #F7B042
  "PosteriorMultimodal", #C27838
  "CinguloOpercular", #B33FB1
  "Dorsalattention",  #72F23F 
  "Auditory", #EC6BFA
  "Visual1",  #2D46F6
  "Somatomotor", #6DFBFB
  "Visual2")) #7847EB


df_cortical[df_cortical ==0] <- NA
df_cortical <- df_cortical[complete.cases(df_cortical),]

pdf("cognition_cortical.pdf", width=4.5,height=4.5)
ggplot(df_cortical, aes(df_cortical, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#ED3833", "#FFF851", "#3C9093", "#5C9A27", "#F7B042", "#C27838", "#B33FB1", "#72F23F", "#EC6BFA", "#2D46F6", "#6DFBFB", "#7847EB")) +
  coord_cartesian(ylim = c(-1.3, 1.5)) 
dev.off()

df_subcortical$variable <- factor(df_subcortical$variable,levels = c(
  "Amygdala", #4E372D
  "Cerebellum", #EDCA51
  "Hippocampus", #EA6741
  "Thalamus",  #42A0B1
  "Striatum")) #CD3036
df_subcortical[df_subcortical ==0] <- NA
df_subcortical <- df_subcortical[complete.cases(df_subcortical),]

pdf("cognition_subcortical.pdf", width=2.3,height=4.5)
ggplot(df_subcortical, aes(df_subcortical, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#4E372D", "#EDCA51", "#EA6741", "#42A0B1", "#CD3036")) +
  coord_cartesian(ylim = c(-1.3, 1)) 
dev.off()




# 2. association/sensory graph

df2 <- df_masked[,c(1,7:18)]
colnames(df2)

df2$assoc <- rowSums(df2[ , c(5,6,7,8,10,11,12,13)],  na.rm=TRUE)
df2$sens <- rowSums(df2[ , c(2,3,4,9)], na.rm=TRUE)

t.test(df2$assoc, df2$sens)
# t = -9.2067, df = 1367.7, p-value < 2.2e-16

# select columns for delta bar graph 1, 36-52
df_net <- df2 %>% select(Parcels, assoc, sens)

# melt
df_net <- melt(df_net, id=c("Parcels"))

# define order 
df_net$variable <- factor(df_net$variable,levels = c(
  "assoc", #ED3833
  "sens")) #7847EB

pdf("cognition_assoc_sens.pdf", width=2,height=4.5)
ggplot(df_net, aes(df_net, x=variable, y=value, fill=variable, width=0.85)) + 
  stat_summary(fun="mean", geom="bar", position=position_dodge(width=0.95), colour='white', size=0.5) + 
  stat_summary(fun.data="mean_se", geom="linerange", fun.args = list(mult = 1), size = 2, position=position_dodge(0.9), color = "black") +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#666666", "#A9A9A9")) +
  coord_cartesian(ylim = c(-0.35, 0.2)) 
dev.off()


