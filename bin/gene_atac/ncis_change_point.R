library(LSD)
library(changepoint)

### get parameters
args = commandArgs(trailingOnly=TRUE)
ncis_table_list_file = args[1]
t_thresh_output_file = args[2]

input_folder = args[3]
output_folder = args[4]

changepoint_method = args[5]

### read ncis table file names
ncis_table_list = read.table(ncis_table_list_file, header=FALSE)

### initialize t table
t_thresh = c()

for (i in c(1:dim(ncis_table_list)[1])){
	### get input & output files
	input_file = paste(input_folder, ncis_table_list[i,3], sep='')
	output_name = paste(ncis_table_list[i,1], ncis_table_list[i,2], sep='_')
	print(output_name)

	### read r vs t table of r vs t: r=sum(x1)/sum(x2); t=x1+x2
	data = read.table(input_file,header = F)
	data_m = data#[data[,1]>10,]

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

	### variance change-point without polynomial regression norm
	print('find variance change-point without polynomial regression norm')
	ansvar=cpt.var(r, class=FALSE, method = changepoint_method, penalty = 'BIC', Q=5)

	### variance change-point with polynomial regression norm
	print('find variance change-point with polynomial regression norm')
	if (max(r)!=0){
		ansvar_norm=cpt.var(r-lo_fit_value, class=FALSE, method = changepoint_method, penalty = 'BIC', Q=5)
		print(ansvar_norm)
		}	else{
			ansvar_norm=c(0)
		}
	

	### plot r vs t pattern without polynomial regression norm & and variance change point
	png(paste(output_folder, output_name, '.png', sep =''))
	heatscatter(t, r, pch = 20, ylim=c(-0.5,0.5), xlim=c(1,10000), log='x', main=paste(toString(ansvar[1]), 'VS', toString(ansvar_norm[1]), sep=' ') )
	lines(t, lo_fit_value, col = 'blue', lty=2)
	abline(v = t[ansvar[1]], col = 'gray', lty=2)
	abline(v = t[ansvar_norm[1]], col = 'red', lty=2)
	dev.off()

	### plot r vs t pattern with polynomial regression norm & and variance change point
	png(paste(output_folder, output_name, '.pn.png', sep =''))
	heatscatter(t, r-lo_fit_value, pch = 20, ylim=c(-1,1), xlim=c(1,10000), log='x', main=toString(ansvar_norm[1]))
	abline(v = t[ansvar_norm[1]], col = 'red', lty=2)
	dev.off()


	if (ansvar_norm[1] != 0){
		t_thresh[i] = t[ansvar_norm[1]]
	} else{
		t_thresh[i] = 0
	}
	
}

### write t_thresh table
t_thresh = cbind(ncis_table_list, t_thresh )
write.table(t_thresh, t_thresh_output_file, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

