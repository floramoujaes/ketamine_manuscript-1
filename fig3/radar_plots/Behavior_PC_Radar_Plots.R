


# the code we are using is an amended version of radarchart function from fmsb library:
# install
# install.packages("fmsb")
# Library
library(fmsb)

# N.B. NEED TO RUN THE radarchartvar AT THE BOTTOM TO WORK
# tsv to read in it in /gpfs/loomis/pi/n3/Studies/Anticevic.DP5/analysis/n-bridge/behavdelta_neurodelta_GSR_mVWMWB1d/analysis/results

# set working directory to where BSNIP files are located
#setwd("~/Dropbox/DP5_N-Bridge/Figure_Calculations/GSR/n-bridge/behavdelta_neurodelta_GSR_mVWMWB1d")
setwd("~/Dropbox/DP5_N-Bridge/Figure_Calculations/GSR/N-Bridge/behavdelta_neurodelta_GSR_mVWMWB1d/behaviour_graphs")
# read in dataframe PCs & raw data
df <- read.table(file = 'NBRIDGE_DP5_behavdelta_neurodelta_GSR_mVWMWB1d_BehaviorPCALoadings.tsv', sep = '\t', header = TRUE)

# read in means for each participant that have been z scored (available in n-bridge folder, analysis, results, assembled outputs)
mean1 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-005.tsv', sep = '\t', header = TRUE)
mean2 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-006.tsv', sep = '\t', header = TRUE)
mean3 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-013.tsv', sep = '\t', header = TRUE)
mean4 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-015.tsv', sep = '\t', header = TRUE)
mean5 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-020.tsv', sep = '\t', header = TRUE)
mean6 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-022.tsv', sep = '\t', header = TRUE)
mean7 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-024.tsv', sep = '\t', header = TRUE)
mean8 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-026.tsv', sep = '\t', header = TRUE)
mean9 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-027.tsv', sep = '\t', header = TRUE)
mean10 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-031.tsv', sep = '\t', header = TRUE)
mean11 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-033.tsv', sep = '\t', header = TRUE)
mean12 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-036.tsv', sep = '\t', header = TRUE)
mean13 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-042.tsv', sep = '\t', header = TRUE)
mean14 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-043.tsv', sep = '\t', header = TRUE)
mean15 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-049.tsv', sep = '\t', header = TRUE)
mean16 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-050.tsv', sep = '\t', header = TRUE)
mean17 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-052.tsv', sep = '\t', header = TRUE)
mean18 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-068.tsv', sep = '\t', header = TRUE)
mean19 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-069.tsv', sep = '\t', header = TRUE)
mean20 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-070.tsv', sep = '\t', header = TRUE)
mean21 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-074.tsv', sep = '\t', header = TRUE)
mean22 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-075.tsv', sep = '\t', header = TRUE)
mean23 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-076.tsv', sep = '\t', header = TRUE)
mean24 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-084.tsv', sep = '\t', header = TRUE)
mean25 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-088.tsv', sep = '\t', header = TRUE)
mean26 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-091.tsv', sep = '\t', header = TRUE)
mean27 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-092.tsv', sep = '\t', header = TRUE)
mean28 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-104.tsv', sep = '\t', header = TRUE)
mean29 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-106.tsv', sep = '\t', header = TRUE)
mean30 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-110.tsv', sep = '\t', header = TRUE)
mean31 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-111.tsv', sep = '\t', header = TRUE)
mean32 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-114.tsv', sep = '\t', header = TRUE)
mean33 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-117.tsv', sep = '\t', header = TRUE)
mean34 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-120.tsv', sep = '\t', header = TRUE)
mean35 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-124.tsv', sep = '\t', header = TRUE)
mean36 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-127.tsv', sep = '\t', header = TRUE)
mean37 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-129.tsv', sep = '\t', header = TRUE)
mean38 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-133.tsv', sep = '\t', header = TRUE)
mean39 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-139.tsv', sep = '\t', header = TRUE)
mean40 <- read.table(file = 'DP5_delta_delta_CAN-NP_V2_sub-DP-5-142.tsv', sep = '\t', header = TRUE)

# create list of dataframes
dflist <- list(mean1,mean2,mean3,mean4,mean5,mean6,mean7,mean8,mean9,mean10,mean11,mean12,mean13,mean14,mean15,mean16,mean17,mean18,mean19,mean20,mean21,mean22,mean23,mean24,mean25,mean26,mean27,mean28,mean29,mean30,mean31,mean32,mean33,mean34,mean35,mean36,mean37,mean38,mean39,mean40)

# create empty matrix
data <- data.frame(matrix(NA, nrow = 31, ncol = 0))

# add in the z scores for each subject
for(i in dflist) {                                
  # new <- as.data.frame(i[, -c(1:5)]) # z score
  new <- as.data.frame(i[, 5]) # raw score
  data[ , ncol(data) + 1] <- new                
}

write.csv(data, 'data.csv')

# find mean of each column then transpose
mean_data <- as.data.frame(rowMeans(data))
mean_data_t <- as.data.frame(t(mean_data))
mean_data_scaled <- mean_data_t/10

# separate out the two PCs
df_PC1 <- df[1,]
# invert
df_PC1 <- df_PC1 * -1
# invert cognition as we have done in the bar charts
df_PC1$BACS_Comp <- df_PC1$BACS_Comp * -1

df_PC2 <- df[2,]
# invert cognition as we have done in the bar charts
df_PC2$BACS_Comp <- df_PC2$BACS_Comp * -1

# define the coloumn names
colnames(df_PC1) <- c("Cognition" , "Delusions" , "Conceptual Disorganization" , "Hallucinations" , "Excitement", "Grandiosity" , "Suspiciousness & Persecution" , "Hostility", "Blunted Affect", "Emotional Withdrawal", "Poor Rapport" , "Passive/Apapthetic Social Withdrawal" , "Abstract Thought" , "Lack of Spontaneity", "Stereotyped Thought" , "Somatic Concern" , "Anxiety", "Guilt Feelings", "Tension" , "Mannerisms & Posturing" , "Depression" , "Motor Retardation" , "Uncooperativeness", "Unusual Thought Content" , "Disorientation" , "Poor Attention", "Lack of Judgement & Insight", "Disturbance of Volition", "Impulse Control" , "Preoccupation" , "Social Avoidance")
colnames(df_PC2) <- c("Cognition" , "Delusions" , "Conceptual Disorganization" , "Hallucinations" , "Excitement", "Grandiosity" , "Suspiciousness & Persecution" , "Hostility", "Blunted Affect", "Emotional Withdrawal", "Poor Rapport" , "Passive/Apapthetic Social Withdrawal" , "Abstract Thought" , "Lack of Spontaneity", "Stereotyped Thought" , "Somatic Concern" , "Anxiety", "Guilt Feelings", "Tension" , "Mannerisms & Posturing" , "Depression" , "Motor Retardation" , "Uncooperativeness", "Unusual Thought Content" , "Disorientation" , "Poor Attention", "Lack of Judgement & Insight", "Disturbance of Volition", "Impulse Control" , "Preoccupation" , "Social Avoidance")
colnames(mean_data_scaled) <- c("Cognition" , "Delusions" , "Conceptual Disorganization" , "Hallucinations" , "Excitement", "Grandiosity" , "Suspiciousness & Persecution" , "Hostility", "Blunted Affect", "Emotional Withdrawal", "Poor Rapport" , "Passive/Apapthetic Social Withdrawal" , "Abstract Thought" , "Lack of Spontaneity", "Stereotyped Thought" , "Somatic Concern" , "Anxiety", "Guilt Feelings", "Tension" , "Mannerisms & Posturing" , "Depression" , "Motor Retardation" , "Uncooperativeness", "Unusual Thought Content" , "Disorientation" , "Poor Attention", "Lack of Judgement & Insight", "Disturbance of Volition", "Impulse Control" , "Preoccupation" , "Social Avoidance")

# To use the fmsb package, I have to add 2 lines to the dataframe: the max (1) and min (-1) of each topic to show on the plot!
df_PC1_r <- rbind(rep(0.5,31) , rep(-0.5,31) , mean_data_scaled, df_PC1)
df_PC2_r <- rbind(rep(0.5,31) , rep(-0.5,31) , mean_data_scaled, df_PC2)
mean_data_scaled_r <- rbind(rep(0.5,31) , rep(-0.5,31) , mean_data_scaled)




# define colours of each of the measures (removed first 5 green colours below as only one cog measure)
# "#009E5D","#004E2E","#005F38","#006E45","#008953",
Sxcols=c("#00B669",
         "#441978","#511D90","#6124AF","#6F26C3","#7C2CE2","#8E35F5","#A351FF",
         "#354252","#334862","#37506E","#406084","#41658F","#446C9C","#4572A6",
         "#5B013B","#6A0144","#7F0154","#950264","#AD016D","#C2017A","#CC0285","#D3038C","#DF0293","#ED039F",
         "#F702A6","#FF3DB9","#FF5ABE","#FF7FD4","#FFA1E4","#FEB5EB")

Sxcolslabels=c("Cognition" , "Delusions" , "Conceptual Disorganization" , "Hallucinations" , 
               "Excitement", "Grandiosity" , "Suspiciousness & Persecution" , "Hostility", 
               "Blunted Affect", "Emotional Withdrawal", "Poor Rapport" , "Passive/Apapthetic Social Withdrawal" , 
               "Abstract Thought" , "Lack of Spontaneity", "Stereotyped Thought" , "Somatic Concern" , "Anxiety", 
               "Guilt Feelings", "Tension" , "Mannerisms & Posturing" , "Depression" , "Motor Retardation" , 
               "Uncooperativeness", "Unusual Thought Content" , "Disorientation" , "Poor Attention", 
               "Lack of Judgement & Insight", "Disturbance of Volition", "Impulse Control" , "Preoccupation" , "Social Avoidance")

numericlabels=c("1" , "2" , "3" , "4" , "5", "6" , "7" , "8", "9", "10", "11" , "12" , 
               "13" , "14", "15" , "16" , "17", "18", "19" , "20" , "21" , "22" , 
               "23", "24" , "25" , "26", "27", "28", "29" , "30" , "31")

# RADARCHART PARAMETERS
# pcol → line color
# pfcol → fill color
# plwd → line width 
# cglcol → color of the net
# cglty → net line type (see possibilities)
# axislabcol → color of axis labels
# caxislabels → vector of axis labels to display
# cglwd → net width 
# vlcex → group labels size
# vlabcol -> label colours
# plot using radarchartvar defined below
# vlabels is where you can insert the labels but may need to input their colour as well


# save graphs

pdf("DP5_PC1_delta_invert.pdf", width=8,height=8)
radarchartvar(df_PC1_r, axistype=0, cglty=1,cglwd=2, vlabels=c(numericlabels), vlabcol =c(Sxcols), vlcex = 2, seglty = 1, seglwd = 2, segcol="grey", cglcol=c(Sxcols), plwd=c(8,8), seg=2, pcol=c("#B4B4B4","#666666"), plty=c(1,1))
dev.off()

pdf("DP5_PC2_delta.pdf", width=8,height=8)
radarchartvar(df_PC2_r, axistype=0, cglty=1,cglwd=2, vlabels=c(numericlabels), vlabcol =c(Sxcols), vlcex = 2, seglty = 1, seglwd = 2, segcol="grey", cglcol=c(Sxcols), plwd=c(8,8), seg=2, pcol=c("#B4B4B4","#666666"), plty=c(1,1))
dev.off()


# radar chart function
# to change middle 0 bar colour change latcol to grey
# to change width of the line change the following line lwd = 4.5 to lwd = 1
# polygon(xx * ((seg/2) + CGap)/(seg + CGap), yy * ((seg/2) + CGap)/(seg + CGap), lty = 1, lwd = 1, border = latcol)
# if (centerzero) {
  
radarchartvar <- function(df, latcol = "grey", latlty = 1, axistype = 0, seg = 4,
                          pty = 16, pcol = 1:8, plty = 1:6, plwd = 1,
                          pdensity = NULL, pangle = 45, pfcol = NA, cglty = 3,
                          cglwd = 1, cglcol = "navy", axislabcol = "black",
                          title = "", maxmin = TRUE, na.itp = TRUE,
                          centerzero = FALSE, vlabels = NULL, vlabcol = "black", vlcex = NULL,
                          caxislabels = NULL, calcex = NULL, paxislabels = NULL,
                          seglty = 1, seglwd = 3, segcol = "black",
                          segmaxcol = "black", segmincol = "black", palcex = NULL, ...) {
  # initial checks
  if (!is.data.frame(df)) {
    cat("The data must be given as a dataframe.\n")
    return()
  }
  if ((n <- length(df)) < 3) {
    cat("The number of variables must be 3 or more.\n")
    return()
  }
  if (maxmin == FALSE) {
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  # plot
  plot(c(-1.2, 1.2),
       c(-1.2, 1.2),
       type = "n",
       frame.plot = FALSE,
       axes = FALSE,
       xlab = "",
       ylab = "",
       main = title,
       cex.main = 2,
       asp = 1,
       ...
  )
  # create polygons
  theta <- seq(90, 450, length = n + 1) * pi / 180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  CGap <- ifelse(centerzero, 0, 1)
  for (i in 0:seg) {
    polygon(xx * (i + CGap) / (seg + CGap),
            yy * (i + CGap) / (seg + CGap),
            lty = seglty,
            lwd = seglwd,
            border = segcol
    )
     if (axistype == 1 | axistype == 3)
      CAXISLABELS <- seq(-1,1,length=seg+1)[i+1]
      CAXISLABELS <- paste(i/seg * 100, "(%)")
     if (axistype == 4 | axistype == 5)
      CAXISLABELS <- sprintf("%3.2f", i/seg)
     if (!is.null(caxislabels) & (i < length(caxislabels)))
      CAXISLABELS <- caxislabels[i + 1]
     if (axistype == 1 | axistype == 3 | axistype == 4 | axistype ==
        5) {
      if (is.null(calcex))
        text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS,
             col = axislabcol)
      else text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS,
                col = axislabcol, cex = calcex)
     }
  }
   polygon(xx * ((seg/2) + CGap)/(seg + CGap), yy * ((seg/2) + CGap)/(seg + CGap), lty = 1, lwd = 1, border = latcol)
  if (centerzero) {
    arrows(0, 0, xx * 1, yy * 1,
           lwd = cglwd, lty = cglty,
           length = 0, col = cglcol
    )
  }
  else {
    arrows(xx / (seg + CGap), yy / (seg + CGap), xx * 1, yy *
             1, lwd = cglwd, lty = cglty, length = 0, col = cglcol)
  }
  # p axis labels
  PAXISLABELS <- df[1, 1:n]
  if (!is.null(paxislabels)) {
    PAXISLABELS <- paxislabels
  }
  if (axistype == 2 | axistype == 3 | axistype == 5) {
    if (is.null(palcex)) {
      text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol)
    } else {
      text(xx[1:n], yy[1:n], PAXISLABELS,
           col = axislabcol,
           cex = palcex
      )
    }
  }
  # vlabels
  VLABELS <- colnames(df)
  if (!is.null(vlabels)) 
    VLABELS <- vlabels
  if (is.null(vlcex)) 
    text(xx * 1.2, yy * 1.2, VLABELS, col = vlabcol)
  else text(xx * 1.2, yy * 1.2, VLABELS, cex = vlcex, col = vlabcol)
  series <- length(df[[1]])
  # series
  series <- length(df[[1]])
  SX <- series - 2
  if (length(pty) < SX) {
    ptys <- rep(pty, SX)
  }
  else {
    ptys <- pty
  }
  if (length(pcol) < SX) {
    pcols <- rep(pcol, SX)
  }
  else {
    pcols <- pcol
  }
  if (length(plty) < SX) {
    pltys <- rep(plty, SX)
  }
  else {
    pltys <- plty
  }
  if (length(plwd) < SX) {
    plwds <- rep(plwd, SX)
  }
  else {
    plwds <- plwd
  }
  if (length(pdensity) < SX) {
    pdensities <- rep(pdensity, SX)
  }
  else {
    pdensities <- pdensity
  }
  if (length(pangle) < SX) {
    pangles <- rep(pangle, SX)
  }
  else {
    pangles <- pangle
  }
  if (length(pfcol) < SX) {
    pfcols <- rep(pfcol, SX)
  }
  else {
    pfcols <- pfcol
  }
  # c axis labels
  for (i in 0:seg) {
    if (axistype == 1 | axistype == 3) {
      CAXISLABELS <- seq(-1, 1, length = seg + 1)[i + 1]
    }
     CAXISLABELS <- paste(i/seg * 100, "(%)")
    if (axistype == 4 | axistype == 5) {
      CAXISLABELS <- sprintf("%3.2f", i / seg)
    }
    if (!is.null(caxislabels) & (i < length(caxislabels))) {
      CAXISLABELS <- caxislabels[i + 1]
    }
    if (axistype == 1 | axistype == 3 | axistype == 4 | axistype == 5) {
      if (is.null(calcex)) {
        text(-0.05, (i + CGap) / (seg + CGap), CAXISLABELS, col = axislabcol)
      } else {
        text(-0.05, (i + CGap) / (seg + CGap), CAXISLABELS, col = axislabcol, cex = calcex)
      }
    }
  }
  for (i in 3:series) {
    xxs <- xx
    yys <- yy
    scale <- CGap / (seg + CGap) + (df[i, ] - df[2, ]) / (df[1, ] - df[2, ]) * seg / (seg + CGap)
    if (sum(!is.na(df[i, ])) < 3) {
      cat(sprintf("[NOT ENOUGH DATA] at %d\n%g\n", i, df[i, ]))
    }
    else {
      for (j in 1:n) {
        if (is.na(df[i, j])) {
          if (na.itp) {
            left <- ifelse(j > 1, j - 1, n)
            while (is.na(df[i, left])) {
              left <- ifelse(left > 1, left - 1, n)
            }
            right <- ifelse(j < n, j + 1, 1)
            while (is.na(df[i, right])) {
              right <- ifelse(right < n, right + 1, 1)
            }
            xxleft <- xx[left] * CGap / (seg + CGap) +
              xx[left] * (df[i, left] - df[2, left]) / (df[
                1,
                left
                ] - df[2, left]) * seg / (seg + CGap)
            yyleft <- yy[left] * CGap / (seg + CGap) +
              yy[left] * (df[i, left] - df[2, left]) / (df[
                1,
                left
                ] - df[2, left]) * seg / (seg + CGap)
            xxright <- xx[right] * CGap / (seg + CGap) +
              xx[right] * (df[i, right] - df[2, right]) / (df[
                1,
                right
                ] - df[2, right]) * seg / (seg + CGap)
            yyright <- yy[right] * CGap / (seg + CGap) +
              yy[right] * (df[i, right] - df[2, right]) / (df[
                1,
                right
                ] - df[2, right]) * seg / (seg + CGap)
            if (xxleft > xxright) {
              xxtmp <- xxleft
              yytmp <- yyleft
              xxleft <- xxright
              yyleft <- yyright
              xxright <- xxtmp
              yyright <- yytmp
            }
            xxs[j] <- xx[j] * (yyleft * xxright - yyright *
                                 xxleft) / (yy[j] * (xxright - xxleft) - xx[j] *
                                              (yyright - yyleft))
            yys[j] <- (yy[j] / xx[j]) * xxs[j]
          }
          else {
            xxs[j] <- 0
            yys[j] <- 0
          }
        }
        else {
          xxs[j] <- xx[j] * CGap / (seg + CGap) + xx[j] *
            (df[i, j] - df[2, j]) / (df[1, j] - df[2, j]) *
            seg / (seg + CGap)
          yys[j] <- yy[j] * CGap / (seg + CGap) + yy[j] *
            (df[i, j] - df[2, j]) / (df[1, j] - df[2, j]) *
            seg / (seg + CGap)
        }
      }
      if (is.null(pdensities)) {
        polygon(xxs,
                yys,
                lty = pltys[i - 2],
                lwd = plwds[i - 2],
                border = pcols[i - 2],
                col = pfcols[i - 2]
        )
      }
      else {
        polygon(xxs,
                yys,
                lty = pltys[i - 2],
                lwd = plwds[i - 2],
                border = pcols[i - 2],
                density = pdensities[i - 2],
                angle = pangles[i - 2],
                col = pfcols[i - 2]
        )
      }
      points(xx * scale, yy * scale,
             pch = ptys[i - 2],
             col = pcols[i - 2]
      )
    }
  }
  # polygon(xx * (seg + CGap)/(seg + CGap), yy * (seg + CGap)/(seg + CGap), lty = seglty, lwd = seglwd, border = segmaxcol)
  # polygon(xx * (0 + CGap)/(seg + CGap), yy * (0 + CGap)/(seg + CGap), lty = seglty, lwd = seglwd, border = segmincol)
  # polygon(xx * ((seg/2) + CGap)/(seg + CGap), yy * ((seg/2) + CGap)/(seg + CGap), lty = latlty, lwd = 4.5, border = latcol)
}




