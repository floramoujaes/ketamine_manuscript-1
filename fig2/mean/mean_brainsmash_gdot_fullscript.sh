



from brainsmash.mapgen.base import Base
from brainsmash.mapgen.eval import base_fit
import numpy as np

# right
brain_map_file = "RIGHT_mean_Ket-Pla_HCS40_hpss_res-mVWMWB1d__CAB-NP_P718_BSNIP_parcels_gbc_mFz_pscalar.txt" 
dist_mat_file = "R_GeodesicParcelDistmat.txt"
base = Base(x=brain_map_file, D=dist_mat_file, resample=True)
surrogates_R = base(n=100000) # create surrogates (**HERE**)
np. save("mean_surr_100000_R.npy", surrogates_R)

# left (run at same time in new window)
brain_map_file = "LEFT_mean_Ket-Pla_HCS40_hpss_res-mVWMWB1d__CAB-NP_P718_BSNIP_parcels_gbc_mFz_pscalar.txt" 
dist_mat_file = "L_GeodesicParcelDistmat.txt"
base = Base(x=brain_map_file, D=dist_mat_file, resample=True)
surrogates_L = base(n=100000) # create surrogates
np. save("mean_surr_100000_L.npy", surrogates_L)


conda activate brainsmash
python
import numpy as np
surrogates_R = np.load('mean_surr_100000_R.npy')
surrogates_L = np.load('mean_surr_100000_L.npy')

# create dummy array of 0s for 358 subcortical parcels and append
for i in range(100000):
	x = surrogates_R[i]
	y = surrogates_L[i]
	cortex = np.append(x,y)
	d = np.arange(358) * 0
	wholebrain = np.append(cortex,d)
	np.savetxt(str(i) + "_RL_mean_Surr.txt", wholebrain)

# to check progress:
cd /gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/brainsmash/mean_surrogates_100000
ls -1t | head -5


CASES="GRIN1 GRIN2A GRIN2B GRIN2C GRIN2D GRIN3A GRIN3B GABRA1 GABRA2 GABRA3 GABRA4 GABRA5 GABRA6 GABRB1 GABRB2 GABRB3 GABRD GABRE GABRG1 GABRG2 GABRG3 GABRP GABRQ GABRR1 GABRR2 GABRR3 SST PVALB"
/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot/GABRA5_cortex_pscalar.txt


# ----------------------------------------->
# Correlate genes with brainsmash files and compute p value
# ----------------------------------------->


cd /gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot
screen
R
# varNames <- c("GRIN3A", "GABRA3", "GABRA4", "GABRA5", "GABRA6", "GABRB1", "SST")
# varNames <- c("GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3B", "GABRA1", "GABRA2", "GABRB2", "GABRB3", "GABRD", "GABRE", "GABRG1", "GABRG2", "GABRG3", "GABRP", "GABRQ", "GABRR1", "GABRR2", "GABRR3", "PVALB")

varNames <- c("GABRG3", "GABRP", "GABRQ", "GABRR1", "GABRR2", "GABRR3", "PVALB")
varNames <- c("GABRG2")

for (i in seq_along(varNames)) {
  gene_name = varNames[i]
  file_name = paste0('/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot/', gene_name, '_cortex_pscalar.txt')
  gene = read.table(file_name, sep=",", header=FALSE)
  x <- vector(length = 100000)
  for (j in 1:100000){
	  file_name2 = paste0('/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/brainsmash/mean_surrogates_100000/', j ,'_RL_mean_Surr.txt')
	  surr = read.table(file_name2, sep=",", header=FALSE)
	  r1 = cor(gene, surr)
	  x[j] <- r1
	  }
  file_name3 = paste0(gene_name, '_mean_surr_r_100000.txt')
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
write.csv(df,file="mean_results.csv",row.names=F)





varNames <- c("GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B", "GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6", "GABRB1", "GABRB2", "GABRB3", "GABRD", "GABRE", "GABRG1", "GABRG2", "GABRG3", "GABRP", "GABRQ", "GABRR1", "GABRR2", "GABRR3", "SST", "PVALB")
 
#rm(list = ls())
#gene_name = "GRIN3A"

for (i in seq_along(varNames)) {
	gene_name = varNames[i]
	file_name3 = paste0(gene_name, '_mean_surr_r_100000.txt')
	#write.csv(x,file=file_name3,row.names=F)
	
	### Use distribution of r values to get p values
	# read in brainsmash surrogates for that gene
	file_name4 = paste0('/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot/', gene_name, '_mean_surr_r_100000.txt')
	df = read.table(file_name4, sep=",", header=TRUE)
	names(df) <- c("r")
	
	# calculate the R value for gene & neural pc1
	mean = read.table("/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot/mean_Ket-Pla_HCS40_hpss_res-mVWMWB1d__CAB-NP_P718_BSNIP_parcels_gbc_mFz_pscalar.txt", sep=",", header=FALSE)	
	file_name5 = paste0('/gpfs/project/fas/n3/Studies/Anticevic.DP5/analysis/gdot/', gene_name, '_cortex_pscalar.txt')	
	gene_map = read.table(file_name5, sep=",", header=FALSE)
	
	# 2. subset cortex (1-360 is cortex) 
	mean_cortex <- as.data.frame(mean[1:360,])
	gene_map_cortex <- as.data.frame(gene_map[1:360,])
	
	# 3. rename
	names(mean_cortex) <- c("mean")
	names(gene_map_cortex) <- c("gene")
	
	# 4. correlate
	r = cor(mean_cortex, gene_map_cortex)
	
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
	csv_fname = "mean_results.csv"
	
	# writing row in the csv file
	write.table(row, file = csv_fname, sep = ",",
            append = TRUE, quote = FALSE,
            col.names = FALSE, row.names = FALSE)
}





# fdr correct p values
R
p <- c(0.033784063, 0.451863467, 1.072667731, 0.965099827, 0.933749734, 0.18044243, 0.704448611, 0.002988907, 0.945298214, 0.429040971, 0.259239671, 0.138329842, 0.22797323, 0.037387443, 0.541972491, 0.089823231, 0.781864186, 0.71977918, 0.706484642, 0.358474162, 0, 0.04565322, 0.855253149, 0.070950469, 0.250093757, 0.073410547, 0.116582365, 0.355491094)
round(p.adjust(p, "fdr"), 3)

# fdr corrected p values:
# add back to excel

0.256 
0.703 
1.000 
1.000 
1.000 
0.459 
0.916 
0.042 
1.000 
0.703 
0.518 
0.387
0.518 
0.256 
0.799 
0.314 
0.952 
0.916 
0.916 
0.627 
0.000 
0.256 
0.998 
0.294
0.518 
0.294 
0.363 
0.627









#####


