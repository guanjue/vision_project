library(LSD)

### get parameters
args = commandArgs(trailingOnly=TRUE)
ncis_table_list_file = args[1]
scale_factor_output_file = args[2]

input_folder = args[3]
output_folder = args[4]

sampling_num = as.numeric(args[5])
log = args[6]

### read ncis table file names
ncis_table_list = read.table(ncis_table_list_file, header=FALSE)

### initialize scale factor matrix
sf_all_matrix = c()

for (i in c(1:dim(ncis_table_list)[1])){
	print(ncis_table_list[i,])
	xais_variable_file = paste(input_folder, ncis_table_list[i, 4], sep='')
	yais_variable_file = paste(input_folder, ncis_table_list[i, 5], sep='')
	output_file_name = paste(ncis_table_list[i, 1], ncis_table_list[i, 2], sep='_')
	t_threshold = as.numeric(ncis_table_list[i, 6])

	### read input reads table
	data_x = read.table(xais_variable_file,header = T)
	data_x = as.matrix(data_x[,2])
	#print(head(data_x))
	data_y = read.table(yais_variable_file,header = T)
	data_y = as.matrix(data_y[,2])
	#print(head(data_y))

	### sampling calculate scale factor & plotting
	used_id = sample(length(data_y)[1],sampling_num)
	data_x = data_x[used_id,]
	data_y = data_y[used_id,]

	### get input t value
	data_t = data_x+data_y
	#print(head(data_t))
	#print(head(data_t<=t_threshold))

	### get t value info & t threshold
	print('t value info & t threshold')
	print(sum(data_t>t_threshold))
	print(t_threshold)


	data_x_low_vs_sig = data_x[data_t<=t_threshold]
	data_y_low_vs_sig = data_y[data_t<=t_threshold]
	### remove zero for plotting
	data_x_low_vs_sig_non0 = data_x_low_vs_sig[as.logical((data_x_low_vs_sig!=0) * (data_y_low_vs_sig!=0))]
	data_y_low_vs_sig_non0 = data_y_low_vs_sig[as.logical((data_x_low_vs_sig!=0) * (data_y_low_vs_sig!=0))]
	print('length(data_y_low_vs_sig_non0):')
	print(length(data_y_low_vs_sig_non0))

	### get sig_vs_sig regions
	data_x_sig_vs_sig = data_x[data_t>t_threshold]
	data_y_sig_vs_sig = data_y[data_t>t_threshold]
	### remove zero for plotting
	data_x_sig_vs_sig_non0 = data_x_sig_vs_sig[as.logical((data_x_sig_vs_sig!=0) * (data_y_sig_vs_sig!=0))]
	data_y_sig_vs_sig_non0 = data_y_sig_vs_sig[as.logical((data_x_sig_vs_sig!=0) * (data_y_sig_vs_sig!=0))]
	print('length(data_y_sig_vs_sig_non0):')
	print(length(data_y_sig_vs_sig_non0))


	###### scale factors
	scale_factor = function(data_x, data_y, method){
		### get mean
		if (method=='mean'){
				merge_x = mean(data_x)
				merge_y = mean(data_y)
			} else{
				merge_x = median(data_x)
				merge_y = median(data_y)
			}

		### get scale factor log scale
		sf = merge_y / merge_x

		### normalize y value
		data_x = data_x * sf

		### get std
		std_x = sd(data_x)
		std_y = sd(data_y)

		print(std_y)

		### get M & A for MA plot
		M = log(data_y/1) - log(data_x/1)
		A = 0.5*(log(data_y/1) + log(data_x/1) )

		### remove 0s for plotting
		for_plotting_id = as.logical( (!is.na(M)) * (!is.na(A)) )
		M = M[for_plotting_id]
		A = A[for_plotting_id]

		#print(length(M))
		#print(length(A))
		output = list('sf'=sf, 'M'=M, 'A'=A, 'merge_x'=merge_x, 'merge_y'=merge_y, 'std_x'=std_x, 'std_y'=std_y)
		return(output)
	}

	### original signal M & A
	od_M = log(data_y) - log(data_x) 
	od_A = 0.5*(log(data_y) + log(data_x) )
	### remove 0s for plotting
	for_plotting_id = as.logical( (!is.na(od_M)) * (!is.na(od_A)) )
	od_M = od_M[for_plotting_id]
	od_A = od_A[for_plotting_id]

	### sf od mean
	sf_info = scale_factor(data_x, data_y, 'mean')

	### sf od median
	sf_info_median = scale_factor(data_x, data_y, 'median')

	### scale low_vs_sig
	sf_info_low_vs_sig = scale_factor(data_x_low_vs_sig, data_y_low_vs_sig, 'mean')

	### scale sig_vs_sig
	sf_info_sig_vs_sig = scale_factor(data_x_sig_vs_sig, data_y_sig_vs_sig, 'mean')

	if (t_threshold > 0){
		### plot sig_vs_sig only
		print('scatter plot')
		png(paste(output_folder, output_file_name, '.sig_vs_sig.hs.png', sep =''))
		if (log == 'T'){
			heatscatter(data_x_sig_vs_sig_non0, data_y_sig_vs_sig_non0, log='xy', pch = 20, ylim=c(1,10000), xlim=c(1,10000))
		} else{
			heatscatter(data_x_sig_vs_sig, data_y_sig_vs_sig, pch = 20, ylim=c(-10,5), xlim=c(20,10000))
		}
		abline(0,1,col = 'blue')
		dev.off()

		### plot low_vs_sig only
		print('scatter plot')
		png(paste(output_folder, output_file_name, '.low_vs_sig.hs.png', sep =''))
		if (log == 'T'){
			heatscatter(data_x_low_vs_sig_non0, data_y_low_vs_sig_non0, log='xy', pch = 20, ylim=c(1,10000), xlim=c(1,10000))
		} else{
			heatscatter(data_x_low_vs_sig, data_y_low_vs_sig, pch = 20, ylim=c(-10,5), xlim=c(20,10000))
		}
		abline(0,1,col = 'blue')
		dev.off()


		###### MA plot for different strategy
		png(paste(output_folder, output_file_name, '.4.MA.png', sep =''))
		par(mfrow=c(2,2))
		###
		heatscatter(od_A, od_M, pch = 20, ylim=c(-10,10), main='original Signal')#, xlim=c(20,10000))
		abline(h=0,col = 'red')
		###
		heatscatter(sf_info$A, sf_info$M, pch = 20, ylim=c(-10,10), main='original Signal mean')#, xlim=c(20,10000))
		abline(h=0,col = 'red')
		###
		heatscatter(sf_info_low_vs_sig$A, sf_info_low_vs_sig$M, pch = 20, ylim=c(-10,10), main='low vs sig mean')#, xlim=c(20,10000))
		abline(h=0,col = 'red')
		###
		heatscatter(sf_info_sig_vs_sig$A, sf_info_sig_vs_sig$M, pch = 20, ylim=c(-10,10), main='sig vs sig mean')#, xlim=c(20,10000))
		abline(h=0,col = 'red')
		###
		dev.off()	

		### get scale factor vector
		sf_vector = c(sf_info$merge_x, sf_info_median$merge_x, sf_info_low_vs_sig$merge_x, sf_info_sig_vs_sig$merge_x,   sf_info$merge_y, sf_info_median$merge_y, sf_info_low_vs_sig$merge_y, sf_info_sig_vs_sig$merge_y,   sf_info$std_x, sf_info_median$std_x, sf_info_low_vs_sig$std_x, sf_info_sig_vs_sig$std_x,   sf_info$std_y, sf_info_median$std_y, sf_info_low_vs_sig$std_y, sf_info_sig_vs_sig$std_y )
	} else{
		sf_vector = c(sf_info$merge_x, sf_info_median$merge_x, sf_info_low_vs_sig$merge_x, sf_info_sig_vs_sig$merge_x,   sf_info$merge_y, sf_info_median$merge_y, sf_info_low_vs_sig$merge_y, sf_info_sig_vs_sig$merge_y,   sf_info$std_x, sf_info_median$std_x, sf_info_low_vs_sig$std_x, sf_info_sig_vs_sig$std_x,   sf_info$std_y, sf_info_median$std_y, sf_info_low_vs_sig$std_y, sf_info_sig_vs_sig$std_y )
	}

	### get low_vs_sig regions

	### append sf to sf_all_matrix
	sf_all_matrix = rbind(sf_all_matrix, sf_vector)
}

### write sf_all_matrix 
print('sf_all_matrix')
sf_all_matrix = cbind(ncis_table_list, sf_all_matrix )
print('write output')
write.table(sf_all_matrix, scale_factor_output_file, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

