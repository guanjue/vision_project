### get parameters
args = commandArgs(trailingOnly=TRUE)
sig_matrix_inputfile = args[1]
sample_info_list = args[2]
output_name = args[3]

###### Load DEseq2 
library(DESeq2)

### reads signal matrix
print('reads data')
data = read.table(sig_matrix_inputfile, header = F, sep = '\t')
data_sample_info = read.table(sample_info_list, header = F)

### read column name & condition
col_label = data_sample_info[,1]
condition = data_sample_info[,2]

### extract signal infomation
data_sig = data[,c(-1,-2,-3,-4,-5)]
data_info = data[,c(1:5)]
rownames(data_sig) = data[,4]
colnames(data_sig) = col_label

### get cell types
condition = factor(condition)

### convert to counts data
countdata = round(data_sig*100, 0)

### convert to data.frame
coldata = data.frame(row.names=colnames(countdata), condition)

### reads signal for DEseq2
dds = DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~condition)

### Regularized log transformation for clustering/heatmaps, etc
print('rlog')
rld = rlog(dds)
### extract rlog matrix
rld_matrix = assay(rld)

### write rlog transformed results
write.table(cbind(data_info, rld_matrix), paste(output_name, 'rld_matrix.txt', sep='.'), quote=F, sep='\t', row.names = FALSE, col.names = FALSE)


### DEseq analysis
print('DEseq scale factor')
dds = DESeq(dds)
### get normalized signal matrix
norm_matrix = counts(dds, normalized=T)
### write scale factor normalized results
write.table(cbind(data_info, log2(norm_matrix+1)), paste(output_name, 'log2_norm_matrix_plus1.txt', sep='.'), quote=F, sep='\t', row.names = FALSE, col.names = FALSE)



