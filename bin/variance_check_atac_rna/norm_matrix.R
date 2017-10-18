library(LSD)

### get parameters
args = commandArgs(trailingOnly=TRUE)
sig_matrix_inputfile = args[1]
scale_factor_list = args[2]
output_name = args[3]

### read signal matrix
data_atac = read.table(sig_matrix_inputfile, header = T, sep='\t')
print(dim(data_atac))
print(head(data_atac))

norm_info = read.table(scale_factor_list)
print(head(norm_info))
### get scale factor
scale_factor = norm_info[,14] / norm_info[,10]
#scale_factor = norm_info[,11] / norm_info[,7]
print(scale_factor)
scale_factor_reorder = c(1, 2, 4, 3, 5, 6, 11, 12, 7, 8, 9, 10, 13, 14, 15, 16)

scale_factor = scale_factor[scale_factor_reorder]
print(scale_factor)
print(length(scale_factor))
### atac matrices
data_atac_sig = data_atac[,-1]
print(head(data_atac_sig))
### scale factor norm
data_atac_sig_sf = t(apply(data_atac_sig, 1, FUN=function(x) x*scale_factor))
print(head(data_atac_sig_sf))

### get output matrix
header = colnames(data_atac)
cREs_id = data_atac[,1]
sig_matrix = cbind(cREs_id, data_atac_sig_sf)
output_matrix = rbind(header, sig_matrix)

### write output file 
write.table(output_matrix, output_name, quote=F, sep='\t', row.names = FALSE, col.names = FALSE)
