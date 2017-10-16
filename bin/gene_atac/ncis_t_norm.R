library(LSD)
library(changepoint)

### get parameters
args = commandArgs(trailingOnly=TRUE)
ncis_table_list_file = args[1]
atac_TSS_matrix = args[2]
atac_TSS_matrix_norm_output = args[3]
input_folder = args[4]


### read ncis table file names
ncis_table_list = read.table(ncis_table_list_file, header=FALSE)

### read TSS matrix
atac_TSS_matrix_all = read.table(atac_TSS_matrix, header=TRUE, sep='\t')
### extract signal matrix
atac_TSS_matrix_sig = atac_TSS_matrix_all[,c(-1,-4,-5,-6)]
atac_TSS_matrix_sig = as.matrix(atac_TSS_matrix_sig)
head(atac_TSS_matrix_sig)

### standarize the interval length
atac_TSS_matrix_sig_200 = t(apply(atac_TSS_matrix_sig, 1, FUN=function(x) x/(x[2]-x[1])*200))
### remove coordinates
atac_TSS_matrix_sig_200 = atac_TSS_matrix_sig_200[,c(-1,-2)]

### get info matrix for output matrix
atac_TSS_matrix_norm = atac_TSS_matrix_all[,c(1,2,3,4,5,6)]
print(dim(atac_TSS_matrix_sig_200))

for (i in c(1:dim(ncis_table_list)[1])){
	### get input & output files
	input_file = paste(input_folder, ncis_table_list[i,3], sep='')
	output_name = paste(ncis_table_list[i,1], ncis_table_list[i,2], sep='_')
	print(output_name)

	### read r vs t table of r vs t: r=sum(x1)/sum(x2); t=x1+x2
	data = read.table(input_file,header = F)
	data_m = data[data[,1]>10,]

	### log2 transform r
	t_od = (data_m[,1])
	r_od = log2(data_m[,2])-log2(data_m[,3])

	### remove x1 OR x2 equals 0
	t = t_od[as.logical((!is.na(r_od)) * (is.finite(r_od)) ) ]
	r = r_od[as.logical((!is.na(r_od)) * (is.finite(r_od)) ) ]
	print(summary(r))

	### polynomial regression fit the log(r) vs log(t) pattern
	print('fit polynomial regression model')
	lo = lm(r~poly(log(t),50, raw=TRUE))
	### get polynomial regression fitted model predicted value
	lo_fit_value = predict(lo, data.frame(x=t))

	###### read TSS atac-signal
	sig1 = atac_TSS_matrix_sig_200[ , ncis_table_list[i,4]-6]
	sig2 = atac_TSS_matrix_sig_200[ , ncis_table_list[i,5]-6]


	### read t value
	t_tss = sig1+sig2
	print(length(t_tss))
	### calculate r value for each gene based on the fitted t-r model
	lo_fit_value_tss = predict(lo, data.frame(x=t_tss))

	### nor sig1 based on t & r  ### r = log2(sig1/sig2) => scale factor of sig1 = 1/(2^r)
	sig1_norm = apply(cbind(sig1, lo_fit_value_tss), 1, FUN=function(x) x[1] / (2^x[2]) )

	### add to norm matrix
	atac_TSS_matrix_norm = cbind(atac_TSS_matrix_norm, sig1_norm)

}

colnames(atac_TSS_matrix_norm) = colnames(atac_TSS_matrix_all)
### write normed table
write.table(atac_TSS_matrix_norm, atac_TSS_matrix_norm_output, quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')

