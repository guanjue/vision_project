library(edgeR)

data = read.table('/Volumes/MAC_Data/data/labs/hardison_lab/vision/gene_rnaseq_atac/rsem_matrix_folder/rsem_matrix.txt', header = F, sep = '\t')
data_sample_info = read.table('/Volumes/MAC_Data/data/labs/hardison_lab/vision/gene_rnaseq_atac/input_folder/rsem_list_sample_col.txt', header = F)

data = read.table('/Volumes/MAC_Data/data/labs/hardison_lab/vision/gene_rnaseq_atac/atac_reads_table/reads_count_matrix_wg_5end_sample_1.txt', header = F, sep = '\t')

### read column name & condition
col_label = data_sample_info[,1]
condition = data_sample_info[,2]

### extract signal infomation
data_sig = data[,c(-1,-2,-3,-4,-5)]
data_sig = data[,c(-1,-2,-3)]

data_info = data[,c(1:5)]
rownames(data_sig) = data[,4]

countdata = as.matrix(round(data_sig*100, 0))


group = factor(c(1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,15,16))

y = DGEList(counts=countdata,group=group)

y_cpm = cpm(y)*1000000

y = calcNormFactors(y, method='TMM')

y_norm_cpm = cpm(y)*1000000


write.table(y$samples$norm.factors, 'edgeR_scale_factor.txt', quote=F, sep='\t', row.names = FALSE, col.names = FALSE)


