library(LSD)
data_atac = read.table('gene_atac/gencode_pc_sort.atac.celltype.matched.txt', header = T, sep='\t')
data_rna = read.table('rsem_matrix_folder/rsem_matrix.norm.rld_matrix.celltype.txt', header = T, sep='\t')

data_atac_cmp = data_rna[,8]
data_rna_cmp = data_rna[,7]

data_atac_sig = data_atac[,c(-1,-2,-3,-4,-5,-6)]
data_rna_sig = data_rna[,c(-1,-2,-3,-4,-5)]

data_atac_sig_var = apply(log2(data_atac_sig+1), 1, var)
data_rna_sig_var = apply(data_rna_sig, 1, var)

plot(data_atac_sig_var, data_rna_sig_var)
heatscatter(data_atac_sig_var, data_rna_sig_var, pch = 20)#, ylim=c(-10,5), xlim=c(20,10000))
heatscatter(log(data_atac_sig_var), log(data_rna_sig_var), pch = 20)#, ylim=c(-10,5), xlim=c(20,10000))


data_atac_sig_rank = t(apply(data_atac_sig, 1, rank))
data_rna_sig_rank = t(apply(data_rna_sig, 1, rank))


data_atac_rna_sig_cor = apply(cbind(data_atac_sig, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method='spearman') )


plot(data_atac_sig_rank_var, data_rna_sig_var)



heatscatter(data_atac_cmp, data_rna_cmp, pch = 20)

