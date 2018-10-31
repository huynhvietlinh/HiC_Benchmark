library(hicrep)
load_matrix <- function(mode, filename, resolution, max_length) {
	if (mode == 0)	{ 	# raw data
	        hi_c_mat <- read.matrix(filename, header = FALSE, stringsAsFactors = FALSE)		
	}
	else { 			# model
		parameter_mat <- read.table(filename, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
		beta = parameter_mat[1,2]
		alpha = parameter_mat[2:nrow(parameter_mat) - 1,1]
		insulator = parameter_mat[2:nrow(parameter_mat) - 1,2]
		matrix_size = length(alpha)
		hi_c_mat <- matrix(nrow = matrix_size, ncol = matrix_size);
		
		for (i in 1:matrix_size) {
			for (j in 1:matrix_size) {
				if ((i != j) && (abs(j-i)*resolution <= max_length)) {
					alpha_i <- alpha[i]
					alpha_j <- alpha[j]
                		        alpha_ij_2 <- (alpha_i + alpha_j)/2;
					log_pred <- alpha_ij_2 + beta*log(abs(j-i))
				
					i_min <- min(i,j)
					i_max <- max(i,j)
					mid_insulation <- 0;
					if (i_max > i_min + 1) {
	        	        		for (k in (i_min + 1):(i_max - 1)) {
							mid_insulation = mid_insulation + insulator[k]
						}
					}	              
					hi_c_mat[i,j] <- exp(log_pred - mid_insulation)
				}
			}		
		}
	}
	matrix_size = nrow(hi_c_mat)
	final_hic_mat <- matrix(nrow = matrix_size, ncol = matrix_size + 3)
	final_hic_mat[,1] <- 1	# chr1
	final_hic_mat[,2] <- seq(0,matrix_size - 1)*resolution
	final_hic_mat[,3] <- seq(1, matrix_size)*resolution
	
	return final_hic_mat
}

SCC <- function(mode, filename_1, filename_2, resolution, max_length) {
	
	# load raw data matrix
        HiC_raw_data <- read.table(data_filename, header = FALSE, stringsAsFactors = FALSE)
	#print(dim(HiC_raw_data))
	# load model parameters
	parameter_mat <- read.table(model_filename, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
	beta = parameter_mat[1,2]
	alpha = parameter_mat[2:nrow(parameter_mat) - 1,1]
	insulator = parameter_mat[2:nrow(parameter_mat) - 1,2]
	# compute and format two matrices
	matrix_size = length(alpha)
	#matrix_size = 100

	#print(matrix_size)
	HiC_model_mat <- matrix(nrow = matrix_size, ncol = matrix_size + 3)
	HiC_raw_mat <- matrix(nrow = matrix_size, ncol = matrix_size + 3)
	HiC_model_mat[,1] <- HiC_raw_mat[,1] <- 1	# chr1
	HiC_model_mat[,2] <- HiC_raw_mat[,2] <- seq(0,matrix_size - 1)*resolution
	HiC_model_mat[,3] <- HiC_raw_mat[,3] <- seq(1, matrix_size)*resolution

	for (j in 1:matrix_size) {
		HiC_model_mat[,j+3] <- HiC_raw_mat[,j+3] <- as.numeric(HiC_raw_data[,j])
	}
	for (i in 1:matrix_size) {
		for (j in 1:matrix_size) {
			if ((i != j) && (abs(j-i)*resolution <= max_length)) {
				alpha_i <- alpha[i]
				alpha_j <- alpha[j]
                	        alpha_ij_2 <- (alpha_i + alpha_j)/2;
	                        log_pred <- alpha_ij_2 + beta*log(abs(j-i))
				
				i_min <- min(i,j)
				i_max <- max(i,j)
				mid_insulation <- 0;
				if (i_max > i_min + 1) {
	        	        	for (k in (i_min + 1):(i_max - 1)) {
						mid_insulation = mid_insulation + insulator[k]
					}
				}	              
				HiC_model_mat[i,j+3] <- exp(log_pred - mid_insulation)
			}		
		}
	}
	#print("Pre-process now!")
	x <- data.frame(HiC_raw_mat, stringsAsFactors = FALSE)
	y <- data.frame(HiC_model_mat, stringsAsFactors = FALSE)
	# compare two matrices	
	processed <- prep(x, y, resolution, 2, max_length)
	#print("Compare now!") 
	scc.out <- get.scc(processed, resolution, max_length)
	print(data_filename)
	print(model_filename)
	print(scc.out$scc)
	#print(scc.out$std)
}

SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_mes_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_mes_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_msc_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_msc_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)

SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_mes_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_mes_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_msc_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_msc_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)

SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/mes_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/mes_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/mes_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_mes_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/mes_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_mes_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/mes_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_msc_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/mes_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_msc_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)

SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/mes_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/mes_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/mes_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_mes_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/mes_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_mes_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/mes_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_msc_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/mes_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_msc_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)

SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/msc_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/msc_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/msc_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_mes_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/msc_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_mes_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/msc_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_msc_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/msc_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_msc_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)

SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/msc_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/msc_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/msc_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_mes_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/msc_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_mes_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/msc_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_msc_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/msc_rep2.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_msc_rep2.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)

#SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep1.40Kb.raw.chr1.mat_-1_-1_1_10", 4e4, 2e6)
#SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep1.40Kb.raw.chr1.mat_-1_-1_1_20", 4e4, 2e6)
#SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep1.40Kb.raw.chr1.mat_-1_-1_1_30", 4e4, 2e6)
#SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep1.40Kb.raw.chr1.mat_-1_-1_1_40", 4e4, 2e6)
#SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep1.40Kb.raw.chr1.mat_-1_-1_5_20", 4e4, 2e6)
#SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep1.40Kb.raw.chr1.mat_-1_-1_10_20", 4e4, 2e6)
#SCC("/share/hormozdiarilab/Data/HiC/Schmitt_2016/contact_maps/RAW/primary_cohort/h1_rep1.40Kb.raw.chr1.mat", "../Model/_share_hormozdiarilab_Data_HiC_Schmitt_2016_contact_maps_RAW_primary_cohort_h1_rep1.40Kb.raw.chr1.mat_-1_-1_20_30", 4e4, 2e6)

