


# ----------------------------------------->
# Generate Gene Maps locally & transfer to grace
# ----------------------------------------->

# run locally via bash
cd /Users/floramoujaes/
export PYTHONPATH="/Users/floramoujaes/:$PYTHONPATH"
export PATH="/Applications/workbench/bin_macosx64:$PATH"
python
from gdot.lib.main import GeminiDot
from gdot.lib.cifti import save_dscalar, save_workbench_image
import pandas as pd
gd = GeminiDot()
import numpy as np

cortex_only=True


genes = ["SST", "PVALB"]





# saves as dscalar but if you open will look the same are parcellated as same value for each point in each parcel
for gene in genes:    
    pscalars = gd.gene_map(gene)
    dscalars = gd.upsample_parcellated(pscalars)    
    if cortex_only:
        dscalars = gd.renormalize(structures=['CORTEX_LEFT', 'CORTEX_RIGHT'], dscalars=dscalars)    
    save_dscalar(dscalars, gene + "_cortex")
    filename="/Users/floramoujaes/gdot_dp5_results/" + gene + ".txt"
    np.savetxt(filename,dscalars)

exit()

# parcellate the files & transfer to grace
cd /Users/floramoujaes/gdot/outputs
export PATH=$PATH:/Applications/workbench/bin_macosx64


CASES="SST PVALB"

for CASE in ${CASES}; do wb_command -cifti-parcellate ${CASE}_cortex.dscalar.nii CAB-NP_P718_BSNIP.dlabel.nii COLUMN ${CASE}_cortex.pscalar.nii; done
for CASE in ${CASES}; do wb_command -cifti-convert -to-text ${CASE}_cortex.pscalar.nii ${CASE}_cortex_pscalar.txt; done
for CASE in ${CASES}; do scp /Users/floramoujaes/gdot/outputs/${CASE}_cortex_pscalar.txt ffm5@transfer-grace.hpc.yale.edu:/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot; done


# ----------------------------------------->
# Correlate with neural pc1 as a sanity check
# ----------------------------------------->

# Gene files location:
/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot
# Brainsmash files location:
/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/brainsmash/neural_pc1_surrogates


# first correlate GABRA5 with the neural pc1 map
cd /gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot

# 1. read in dataframes 
R
neural_pc1 = read.table("/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot/PC1_ket-pla_pca_image_ptseries.txt", sep=",", header=FALSE)		
GABRA5 = read.table("/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot/GABRA5_cortex_pscalar.txt", sep=",", header=FALSE)
# 2. subset cortex (1-360 is cortex) 
neural_pc1_cortex <- as.data.frame(neural_pc1[1:360,])
GABRA5_cortex <- as.data.frame(GABRA5[1:360,])
# 3. rename
names(neural_pc1_cortex) <- c("neuralpc1")
names(GABRA5_cortex) <- c("gene")
# 4. correlate
cor(neural_pc1_cortex, GABRA5_cortex)
# gives a correlation equivalent to one in keynote (but in keynote results are sign flipped).


# ----------------------------------------->
# Correlate with brainsmash files and compute p value
# ----------------------------------------->

cp /gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/brainsmash/neural_pc1_surrogates/0_RL_NPC1_Surr.txt /gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/brainsmash/neural_pc1_surrogates/100000_RL_NPC1_Surr.txt


cd /gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot
screen
R

varNames <- c("SST", "PVALB")


for (i in seq_along(varNames)) {
  gene_name = varNames[i]
  file_name = paste0('/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot/', gene_name, '_cortex_pscalar.txt')
  gene = read.table(file_name, sep=",", header=FALSE)
  x <- vector(length = 100000)
  for (j in 1:100000){
	  file_name2 = paste0('/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/brainsmash/neural_pc1_surrogates/', j ,'_RL_NPC1_Surr.txt')
	  surr = read.table(file_name2, sep=",", header=FALSE)
	  r1 = cor(gene, surr)
	  x[j] <- r1
	  }
  file_name3 = paste0(gene_name, '_NPC1_surr_r_100000.txt')
  write.csv(x,file=file_name3,row.names=F)
}



# # add the last one as for some reason R makes the file name 1e+05_RL_NPC1_Surr.txt
# i = 100000
# surr = read.table('/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/brainsmash/mean_surrogates_100000/100000_RL_mean_Surr.txt', sep=",", header=FALSE)
# r1 = cor(gene, surr)
# x[i] <- r1

cd /gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot
R
# clear
rm(list = ls())

## run this first time to create dataframe to add values to
df <- data.frame(matrix(ncol = 5, nrow = 0))
x <- c("gene", "r-value", "min", "max", "p-value")
colnames(df) <- x
write.csv(df,file="NPC1_results.csv",row.names=F)





varNames <- c("SST", "PVALB")
 
#rm(list = ls())
#gene_name = "GRIN3A"

for (i in seq_along(varNames)) {
	gene_name = varNames[i]
	file_name3 = paste0(gene_name, '_NPC1_surr_r_100000.txt')
	#write.csv(x,file=file_name3,row.names=F)
	
	### Use distribution of r values to get p values
	# read in brainsmash surrogates for that gene
	file_name4 = paste0('/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot/', gene_name, '_NPC1_surr_r_100000.txt')
	df = read.table(file_name4, sep=",", header=TRUE)
	names(df) <- c("r")
	
	# calculate the R value for gene & neural pc1
	NPC1 = read.table("/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot/PC1_ket-pla_pca_image_ptseries.txt", sep=",", header=FALSE)	
	file_name5 = paste0('/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot/', gene_name, '_cortex_pscalar.txt')	
	gene_map = read.table(file_name5, sep=",", header=FALSE)
	
	# 2. subset cortex (1-360 is cortex) 
	NPC1_cortex <- as.data.frame(NPC1[1:360,])
	gene_map_cortex <- as.data.frame(gene_map[1:360,])
	
	# 3. rename
	names(NPC1_cortex) <- c("NPC1")
	names(gene_map_cortex) <- c("gene")
	
	# 4. correlate
	r = cor(NPC1_cortex, gene_map_cortex)
	
	# find max and min value
	max = max(df$r, na.rm = TRUE)
	min = min(df$r, na.rm = TRUE)
	# do these incorporate the r value?
	
	y = r[1,1]
	if (r > 0) {
	p_value <- sum(df$r > y)/sum(df$r < y)
	} else {
	p_value <- sum(df$r < y)/sum(df$r > y)
	}
	
	#print(paste0(gene_name, " p-value: ", p_value))
	
	# add values to dataframe
	# defining a row 
	row <- data.frame(gene_name, r, min, max, p_value)
	#row
	
	# sample csv name
	csv_fname = "NPC1_results.csv"
	
	# writing row in the csv file
	write.table(row, file = csv_fname, sep = ",",
            append = TRUE, quote = FALSE,
            col.names = FALSE, row.names = FALSE)
}





# fdr correct p values
R
p <- c(0.843114126, 0.004106797, 0.004399269, 0.015228426, 0.010417403, 0, 0.078911594, 0.58065281, 0.2036881, 0, 0, 0, 0, 0, 4.00E-05, 5.00E-05, 0.223915305, 0.050563627, 0.139315499, 0.00126159, 0.000260068, 0.067668852, 0.106157981, 0, 0.040571898, 2.210994445, 0, 0.014661864)
round(p.adjust(p, "fdr"), 3)


# fdr corrected p values:
# add back to excel

0.874 
0.009 
0.009 
0.025 
0.019 
0.000 
0.105 
0.625 
0.238 
0.000 
0.000 
0.000
0.000 
0.000 
0.000 
0.000 
0.251 
0.075 
0.170 
0.003 
0.001 
0.095 
0.135 
0.000
0.063 
1.000 
0.000 
0.025

