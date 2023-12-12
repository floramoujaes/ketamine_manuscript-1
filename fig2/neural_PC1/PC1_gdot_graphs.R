
library(ggplot2)
library(reshape)
library(ggridges)

theme1 <- theme(panel.border = element_rect(size = 2, colour = 'NA', fill='NA'), 
               axis.line = element_line(colour = "black", size = 1),
               axis.text.x = element_text(size = 8, color="black"), 
               axis.title.x = element_text(size =8,margin=margin(15,0,0,0)),
               axis.text.y = element_text(size = 8, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
               axis.title.y = element_text(size =8,margin=margin(0,15,0,0)),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               plot.background = element_rect(fill = "transparent", colour = NA),
               panel.background = element_rect(fill = "white", colour = NA),
               legend.position="none",
               axis.ticks = element_line(colour = "black", size=1),
               axis.ticks.length = unit(.4, "cm"))

# read in dataframes
setwd("~/Dropbox/DP5_N-Bridge/Figure_Calculations/GSR/genetics/neural_pc1")
df_r = read.table("NPC1_results.csv", sep=",", header=TRUE)

# set order
df_r$gene <- factor(df_r$gene, levels = c("SST","PVALB"))

# bar chart for mean
tiff(file="NPC1_rvalues_gdot_SSTPVALB.tiff",width=1.5,height=4.5, units="in", res=1000)
ggplot(df_r, aes(x=gene, y=r.value, fill=gene)) +
  geom_bar(stat='identity') +
  theme1 +
  coord_cartesian(ylim = c(-0.25, 0.5)) +
  scale_fill_manual(values = c("SST" = "#FC2703", "PVALB" = "#2D33FF"))
dev.off()

#----------------------------->
# distribution of r plots
#----------------------------->
theme1 <- theme(panel.border = element_rect(size = 2, colour = 'NA', fill='NA'), 
                axis.line = element_line(colour = "black", size = 1),
                axis.text.x = element_text(size = 8, color="black"), 
                axis.title.x = element_text(size =8,margin=margin(15,0,0,0)),
                axis.text.y = element_text(size = 8, color="black",angle=90,hjust=0.5,margin=margin(0,10,0,0)), 
                axis.title.y = element_text(size =8,margin=margin(0,15,0,0)),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                plot.background = element_rect(fill = "transparent", colour = NA),
                panel.background = element_rect(fill = "white", colour = NA),
                legend.position="none",
                axis.ticks = element_line(colour = "black", size=1),
                axis.ticks.length = unit(.4, "cm"))


df1 = read.table("SST_NPC1_surr_r_100000.txt", header=TRUE)
df2 = read.table("PVALB_NPC1_surr_r_100000.txt", header=TRUE)

names(df1) <- c("SST")
names(df2) <- c("PVALB")

df_all <- cbind.data.frame(df1, df2)

# melt the data
df_all_melt <- melt(df_all)
names(df_all_melt) <- c("gene", "r")

colnames(df_all_melt)
sapply(df_all_melt, class)


# ridgeline 
# need to update r values for SST/PVALB
tiff(file="NPC1_rdistributions_gdot_SST_PVALB.tiff",width=10, height=6, units="in", res=1000)
ggplot(df_all_melt, aes(r, gene, group = gene)) +
  # add the ridges for each trial type
  geom_density_ridges(aes(fill=gene), lwd=1, rel_min_height = 0.01, scale = 1)+
  # add colours
  scale_fill_manual(values = c("SST"="#E5E5E5", "PVALB"="#E5E5E5")) +
  # SST
  geom_segment(aes(x = -0.466134169, y = 1, xend = -0.466134169, yend = 1.9), size = 2, color = "#FC2703") +
  # PVALB
  geom_segment(aes(x = 0.21479303, y = 2, xend = 0.21479303, yend = 2.9), size = 2, color = "#2D33FF") +
  guides(fill=FALSE) +
  scale_y_discrete(expand = expansion(add = c(0.2, 3))) +
  theme1
dev.off()


