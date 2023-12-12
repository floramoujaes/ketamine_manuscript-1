
install.packages('gplots')

library(ggplot2)
library(Hmisc)
library(gplots)
library(plyr)
library(scales)

# set working directory
setwd("~/Dropbox/DP5_N-Bridge/Figure_Calculations/GSR/N40_Neural_PCA_mVWMWB1d/correlate_5_sig_diffPCs")


# 3. correlation matrix
# 4. scatter plots mean & pc 1/2
# 5. Delta correlate pc1 & mean, pc2 & mean


# ---------------------------------------------------------------------------------------->
### Delta correlate first 5 PCS
# ---------------------------------------------------------------------------------------->

# read dataframes
pc1 = data.frame(read.table("PC1_ket-pla_pca_image_ptseries.txt"))
pc2 = data.frame(read.table("PC2_ket-pla_pca_image_ptseries.txt"))
pc3 = data.frame(read.table("PC3_ket-pla_pca_image_ptseries.txt"))
pc4 = data.frame(read.table("PC4_ket-pla_pca_image_ptseries.txt"))
pc5 = data.frame(read.table("PC5_ket-pla_pca_image_ptseries.txt"))
mean = data.frame(read.table("mean_Ket-Pla_HCS40_hpss_res-mVWMWB1d__CAB-NP_P718_BSNIP_parcels_gbc_mFz_pscalar.txt"))


# sign flip pcs to better compare with mean as sign is arbitrary
pc1 <- pc1 * -1
pc2 <- pc2 * -1
pc3 <- pc3 * -1
pc4 <- pc4 * -1
pc5 <- pc5 * -1

# rename columns headers
pc1 <- rename(pc1, c("V1"="pc1"))
pc2 <- rename(pc2, c("V1"="pc2"))
pc3 <- rename(pc3, c("V1"="pc3"))
pc4 <- rename(pc4, c("V1"="pc4"))
pc5 <- rename(pc5, c("V1"="pc5"))
mean <- rename(mean, c("V1"="mean"))

# combines dataframes
df <- cbind(pc1, pc2, pc3, pc4, pc5, mean)

# check values are numeric 
sapply(df, class) 

# create empty dataframe
df2 <- data.frame(matrix(NA, nrow = 5, ncol = 2))

# add headers
colnames(df2)
df2 <- rename(df2, c("X1"="pc", "X2"="corr_with_mean"))

pc1 <- cor.test(df$"pc1", df$"mean", use = "pairwise",method="pearson")
df2[1,] <- "pc1"
df2[1,2] <- pc1$estimate

pc2 <- cor.test(df$"pc2", df$"mean", use = "pairwise",method="pearson")
df2[2,] <- "pc2"
df2[2,2] <- pc2$estimate

pc3 <- cor.test(df$"pc3", df$"mean", use = "pairwise",method="pearson")
df2[3,] <- "pc3"
df2[3,2] <- pc3$estimate

pc4 <- cor.test(df$"pc4", df$"mean", use = "pairwise",method="pearson")
df2[4,] <- "pc4"
df2[4,2] <- pc4$estimate

pc5 <- cor.test(df$"pc5", df$"mean", use = "pairwise",method="pearson")
df2[5,] <- "pc5"
df2[5,2] <- pc5$estimate

cor.test(df$"pc1", df$"pc2", use = "pairwise",method="pearson")
# t = -0.31367, df = 716, p-value = 0.7539, r + -0.01172177 

sapply(df2, class) 
df2$corr_with_mean <- as.numeric(df2$corr_with_mean)

# round the corr values
df2$corr_with_mean <- round(df2$corr_with_mean ,digit=2)

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


pdf("corr_pcs_mean_bar_plot_v3.pdf", width=3.1,height=6)
ggplot(data=df2, aes(x=pc, y=corr_with_mean)) +
  geom_bar(stat="identity", fill="#666666", width = 0.7, size = 1.2) +
  theme1
dev.off()

p <- c(pc1$p.value, pc2$p.value, pc3$p.value, pc4$p.value, pc5$p.value)

p.adjust(p, method = "bonferroni", n = length(p))
