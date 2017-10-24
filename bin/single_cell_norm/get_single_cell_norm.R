library(scran)
library(LSD)
library(org.Mm.eg.db)
library(DESeq2)

data = read.table('/Volumes/MAC_Data/data/labs/hardison_lab/vision/gene_rnaseq_atac/rsem_matrix_folder/rsem_matrix.txt', header = F, sep = '\t')
data_sample_info = read.table('/Volumes/MAC_Data/data/labs/hardison_lab/vision/gene_rnaseq_atac/input_folder/rsem_list_sample_col.txt', header = F)

### read column name & condition
col_label = data_sample_info[,1]
condition = data_sample_info[,2]

### extract signal infomation
data_sig = data[,c(-1,-2,-3,-4,-5)]
data_info = data[,c(1:5)]
rownames(data_sig) = data[,4]

countdata = as.matrix(round(data_sig*100, 0))

### scran scale factor
sf = computeSumFactors(countdata, size=5)
write.table(t(sf), 'scran_scale_factor.txt', quote=F, sep='\t', row.names = FALSE, col.names = FALSE)

### DESeq2 size factor
coldata = data.frame(row.names=colnames(countdata), condition)
### reads signal for DEseq2
dds = DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~condition)
print('DEseq scale factor')
dds = DESeq(dds)
sf_deseq = dds$sizeFactor

png('compare_deseq_scran.png')
plot(sf, sf_deseq, log="xy", ylab="DESeq2 size factor", xlab="Scran scale factor", xlim=c(0,2.3), ylim=c(0,2.3))
abline(0,1, col='red')
dev.off()

