
# # install
# install.packages("prevalence")
# install.packages("PairedData") # pitman morgan test

# Library
library("data.table")
library('prevalence')
library('Hmisc')
library('PairedData') # pitman morgan test


# set working directory 
setwd("~/Dropbox/DP5_N-Bridge/Figure_Calculations/GSR/N-Bridge/behavdelta_neurodelta_GSR_mVWMWB1d/delta_pca_scores_graph")

# read in dataframe
# dataframe is in analysis/results: NBRIDGE_plaket_GSR_mVWMWB1d_project_PCAProjectedScoresZscores.tsv
# then need to rearrange it (see current)
df <- read.table(file = 'NBRIDGE_DP5_delta_delta_CAN-NP_V2_BehaviorPCAScores.tsv', sep = '\t', header = TRUE)


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

# -------------------------------------->
# 2 GRAPHS - indiviudal subjects bar plots & inset scatter plots
# -------------------------------------->

# invert pc1
#df$PC1 <- df$PC1 * -1

# add id numbers
df$id <- 1:nrow(df) 

# set classes
df$PC1 <- as.numeric(df$PC1)
df$id <- as.factor(df$id)

# use reorder to set order on graph
df$id = with(df, reorder(id, PC1, median))



# graph 
pdf("COMPARISON_Individual_Subject_bar_chart_pc1_gradientbluetoyellow_subnos.pdf", width=12,height=4.5)
ggplot(df, aes(y=PC1, x=id, fill=id)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("22" = "#209EFF", "21" = "#209EFF", "9" = "#37A1E4", "3" = "#37A1E4", "11" = "#37A1E4", "13" = "#37A1E4", "10" = "#42A3D7", "18" = "#42A3D7", "25" = "#4EA5C9", "16" = "#4EA5C9", "14" = "#59A6BC", "1" = "#59A6BC", "2" = "#65A8AE", "4" = "#65A8AE", "5" = "#70AAA1", "23" = "#70AAA1", "17" = "#7BAB94", "8" = "#7BAB94", "30" = "#87AD86", "24" = "#87AD86", "12" = "#92AF79", "36" = "#92AF79", "15" = "#98B072", "37" = "#9EB16B", "20" = "#A9B25E", "32" = "#A9B25E", "33" = "#B4B451", "6" = "#B4B451", "31" = "#C0B643", "26" = "#C0B643", "38" = "#CBB736", "7" = "#CBB736", "40" = "#D7B928", "27" = "#D7B928", "34" = "#E2BB1B", "29" = "#E2BB1B", "28" = "#EEBC0D", "39" = "#EEBC0D", "35" = "#F9BE00", "19" = "#F9BE00")) +
  theme1 + 
  scale_y_continuous(limits=c(-11, 5), breaks=c(-11, 0, 5)) +
  theme(text = element_text(size=50, hjust = 0.05), axis.text.x = element_text(angle = 270, color="black"), axis.text = element_text(color="black"))
dev.off()



# use reorder to set order on graph
df$id = with(df, reorder(id, PC2, median))


pdf("COMPARISON_Individual_Subject_bar_chart_pc2_gradientbluetoyellow_subnos.pdf", width=12,height=4.5)
ggplot(df, aes(y=PC2, x=id, fill=id)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("3" = "#209EFF", "37" = "#209EFF", "2" = "#37A1E4", "5" = "#37A1E4", "13" = "#37A1E4", "19" = "#37A1E4", "39" = "#42A3D7", "31" = "#42A3D7", "21" = "#4EA5C9", "15" = "#4EA5C9", "35" = "#59A6BC", "10" = "#59A6BC", "14" = "#65A8AE", "24" = "#65A8AE", "26" = "#70AAA1", "34" = "#70AAA1", "30" = "#7BAB94", "36" = "#7BAB94", "32" = "#87AD86", "27" = "#87AD86", "8" = "#92AF79", "33" = "#92AF79", "7" = "#9EB16B", "28" = "#9EB16B", "38" = "#A9B25E", "29" = "#A9B25E", "17" = "#B4B451", "40" = "#B4B451", "11" = "#C0B643", "22" = "#C0B643", "6" = "#CBB736", "1" = "#CBB736", "4" = "#D7B928", "9" = "#D7B928", "23" = "#E2BB1B", "16" = "#E2BB1B", "12" = "#EEBC0D", "25" = "#EEBC0D", "18" = "#F9BE00", "20" = "#F9BE00")) +
  theme1 + 
  scale_y_continuous(limits=c(-11, 5), breaks=c(-11, 0, 5)) +
  theme(text = element_text(size=50, hjust = 0.05), axis.text.x = element_text(angle = 270, color="black"), axis.text = element_text(color="black"))
dev.off()

