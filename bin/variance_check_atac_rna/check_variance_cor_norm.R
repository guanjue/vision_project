library(LSD)
library(seewave)

data_atac = read.table('gene_atac/gencode_pc_sort.TSSexp10kb.atac.celltype.matched.txt', header = T, sep='\t')

#data_atac = read.table('gene_atac/gencode_pc_sort.TSSexp10kb.atac.scran.celltype.matched.txt', header = T, sep='\t')
data_atac = read.table('gene_atac/gencode_pc_sort.TSSexp10kb.atac.edger.celltype.matched.txt', header = T, sep='\t')

data_atac_1kb = read.table('gene_atac/gencode_pc_sort.TSSexp1kb.atac.celltype.matched.txt', header = T, sep='\t')
data_atac_5kb = read.table('gene_atac/gencode_pc_sort.TSSexp5kb.atac.celltype.matched.txt', header = T, sep='\t')
data_atac_1mb = read.table('gene_atac/gencode_pc_sort.TSSup1000kb.atac.celltype.matched.txt', header = T, sep='\t')
data_atac_tr = read.table('gene_atac/gencode_pc_sort.TSSexp10kb.atac.celltype.matched.TRnormed.txt', header = T, sep='\t')

data_atac1kb = read.table('gene_atac/gencode_pc_sort.TSSexp1kb.atac.celltype.matched.txt', header = T, sep='\t')
data_atac5kb = read.table('gene_atac/gencode_pc_sort.TSSexp5kb.atac.celltype.matched.txt', header = T, sep='\t')
data_atac10kb = read.table('gene_atac/gencode_pc_sort.TSSexp10kb.atac.celltype.matched.txt', header = T, sep='\t')
data_atac1000kbkb = read.table('gene_atac/gencode_pc_sort.TSSup1000kb.atac.celltype.matched.txt', header = T, sep='\t')

data_atac_tr_1mb = read.table('gene_atac/gencode_pc_sort.TSSup1000kb.atac.celltype.matched.TRnormed.txt', header = T, sep='\t')
data_atac_tr_1kb = read.table('gene_atac/gencode_pc_sort.TSSexp1kb.atac.celltype.matched.TRnormed.txt', header = T, sep='\t')
data_atac_tr_5kb = read.table('gene_atac/gencode_pc_sort.TSSexp5kb.atac.celltype.matched.TRnormed.txt', header = T, sep='\t')
data_atac_tr_10kb = read.table('gene_atac/gencode_pc_sort.TSSexp10kb.atac.celltype.matched.TRnormed.txt', header = T, sep='\t')
print('dim(data_atac_tr_10kb)')
print(dim(data_atac_tr_10kb))

norm_info = read.table('scale_factor_matrix/ncis_table_list.sf.txt')
scale_factor_allmean = norm_info[c(1:12),11] / norm_info[c(1:12),7]
scale_factor = norm_info[c(1:12),14] / norm_info[c(1:12),10]


### rna matrix
print('rna matrix')
data_rna = read.table('rsem_matrix_folder/rsem_matrix.norm.rld_matrix.celltype.matched.txt', header = F, sep='\t')
print(dim(data_rna))
data_rna_sig = data_rna[,c(-1,-2,-3,-4,-5,-6)]
print(dim(data_rna))

data_rna = read.table('rsem_matrix_folder/rsem_matrix.norm.log2_norm_matrix_plus1.celltype.matched.txt', header = F, sep='\t')
print(dim(data_rna))
data_rna_sig = data_rna[,c(-1,-2,-3,-4,-5,-6)]
print(dim(data_rna))


### atac matrices
data_atac_sig = data_atac[,c(-1,-2,-3,-4,-5,-6)]
data_atac_sig_1kb = data_atac_1kb[,c(-1,-2,-3,-4,-5,-6)]
data_atac_sig_5kb = data_atac_5kb[,c(-1,-2,-3,-4,-5,-6)]
data_atac_sig_1mb = data_atac_1mb[,c(-1,-2,-3,-4,-5,-6)]

data_atac_sig_tr = data_atac_tr[,c(-1,-2,-3,-4,-5,-6)]
data_atac_sig_tr_1mb = data_atac_tr_1mb[,c(-1,-2,-3,-4,-5,-6)]
data_atac_sig_tr_1kb = data_atac_tr_1kb[,c(-1,-2,-3,-4,-5,-6)]
data_atac_sig_tr_5kb = data_atac_tr_5kb[,c(-1,-2,-3,-4,-5,-6)]
data_atac_sig_tr_10kb = data_atac_tr_10kb[,c(-1,-2,-3,-4,-5,-6)]

data_atac_sig1 = data_atac1kb[,c(-1,-2,-3,-4,-5,-6)]
data_atac_sig5 = data_atac5kb[,c(-1,-2,-3,-4,-5,-6)]
data_atac_sig10 = data_atac10kb[,c(-1,-2,-3,-4,-5,-6)]
data_atac_sig1000kb = data_atac1000kbkb[,c(-1,-2,-3,-4,-5,-6)]

### scale factor norm
data_atac_sig_sf = t(apply(data_atac_sig, 1, FUN=function(x) x*scale_factor))
data_atac_sig_sf_am = t(apply(data_atac_sig, 1, FUN=function(x) x*scale_factor_allmean))

data_atac_sig_sf1 = t(apply(data_atac_sig1, 1, FUN=function(x) x*scale_factor))
data_atac_sig_sf_am1 = t(apply(data_atac_sig1, 1, FUN=function(x) x*scale_factor_allmean))

data_atac_sig_sf5 = t(apply(data_atac_sig5, 1, FUN=function(x) x*scale_factor))
data_atac_sig_sf_am5 = t(apply(data_atac_sig5, 1, FUN=function(x) x*scale_factor_allmean))

data_atac_sig_sf10 = t(apply(data_atac_sig10, 1, FUN=function(x) x*scale_factor))
data_atac_sig_sf_am10 = t(apply(data_atac_sig10, 1, FUN=function(x) x*scale_factor_allmean))

data_atac_sig_sf1000kb = t(apply(data_atac_sig1000kb, 1, FUN=function(x) x*scale_factor))
data_atac_sig_sf_am1000kb = t(apply(data_atac_sig1000kb, 1, FUN=function(x) x*scale_factor_allmean))

data_rna_sig = as.matrix(data_rna_sig)
class(data_rna_sig) = 'numeric'

png('norm_compare/rna_hist.png')
plot( density(data_rna_sig) )
dev.off()


###### input matrix
### od log2
data_atac_sig_log2 = log2((data_atac_sig+1)/10000*1000)
print(head(data_atac_sig_log2))
print('get used id')
used_id = as.logical( apply((data_atac_sig_log2), 1, FUN=function(x) (min((x)) > log2(0.1)) ) * apply((data_rna_sig), 1, FUN=function(x) (max((x)) > 4) ) )
used_id = as.logical(  apply((data_rna_sig), 1, FUN=function(x) (max((x)) > 4) ) )

print(length(used_id))
print(sum(used_id))
data_atac_sig_log2 = data_atac_sig_log2[used_id, ]
print('dim(data_atac_sig_log2)')
print(dim(data_atac_sig_log2))

data_atac_sig_log2_1kb = log2((data_atac_sig_1kb+1)/10000*1000)[used_id, ]
data_atac_sig_log2_5kb = log2((data_atac_sig_5kb+1)/10000*1000)[used_id, ]
data_atac_sig_log2_1mb = log2((data_atac_sig_1mb+1)/10000*1000)[used_id, ]

data_atac_sig_log2_tr = log2((data_atac_sig_tr+1)/10000*1000)[used_id, ]
data_atac_sig_log2_tr_1mb = log2((data_atac_sig_tr_1mb+1)/10000*1000)[used_id, ]
data_atac_sig_log2_tr_1kb = log2((data_atac_sig_tr_1kb+1)/1000*1000)[used_id, ]
data_atac_sig_log2_tr_5kb = log2((data_atac_sig_tr_5kb+1)/5000*1000)[used_id, ]
data_atac_sig_log2_tr_10kb = log2((data_atac_sig_tr_10kb+1)/10000*1000)[used_id, ]

### scale log2
#data_atac_sig_sf_am_log2 = log2(data_atac_sig_sf_am+1)
### sf_sd_norm log2
#data_atac_sig_sf_log2 = log2(data_atac_sig_sf+1)

### scale log2
print(dim(data_rna_sig))
data_rna_sig=data_rna_sig[used_id, ]
data_rna = data_rna[used_id, ]


data_atac_sig_sf_am_log2_1 = log2((data_atac_sig_sf_am1+1)/1000*1000)[used_id, ]
data_atac_sig_sf_log2_1 = log2((data_atac_sig_sf1+1)/1000*1000)[used_id, ]

data_atac_sig_sf_am_log2_5 = log2((data_atac_sig_sf_am5+1)/5000*1000)[used_id, ]
data_atac_sig_sf_log2_5 = log2((data_atac_sig_sf5+1)/5000*1000)[used_id, ]

data_atac_sig_sf_am_log2_10 = log2((data_atac_sig_sf_am10+1)/10000*1000)[used_id, ]
data_atac_sig_sf_log2_10 = log2((data_atac_sig_sf10+1)/10000*1000)[used_id, ]

data_atac_sig_sf_am_log2_1000kb = log2((data_atac_sig_sf_am1000kb+1)/10000*1000)[used_id, ]
data_atac_sig_sf_log2_1000kb = log2((data_atac_sig_sf1000kb+1)/10000*1000)[used_id, ]


for (method in c('pearson', 'spearman')){
	### spearman correlation
	data_atac_rna_sig_cor = apply(cbind(data_atac_sig_log2, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_rna_sig_cor_1kb = apply(cbind(data_atac_sig_log2_1kb, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_rna_sig_cor_5kb = apply(cbind(data_atac_sig_log2_5kb, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_rna_sig_cor_1mb = apply(cbind(data_atac_sig_log2_1mb, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )


	data_atac_rna_sig_cor_tr = apply(cbind(data_atac_sig_log2_tr, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_rna_sig_cor_tr_1mb = apply(cbind(data_atac_sig_log2_tr_1mb, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_rna_sig_cor_tr_1kb = apply(cbind(data_atac_sig_log2_tr_1kb, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_rna_sig_cor_tr_5kb = apply(cbind(data_atac_sig_log2_tr_5kb, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )

	data_atac_sig_log2_tr_1mb_cor1000kb_shuffle = apply(cbind( data_atac_sig_log2_tr_1mb[sample(nrow(data_atac_sig_log2_tr_1mb)),], data_rna_sig ), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )


	#print(head(data_atac_sig_log2))
	#print(head(data_rna_sig))
	test_data = (data_rna_sig[as.logical((data_atac_rna_sig_cor<=(-0.4)) * (data_atac_rna_sig_cor>=(-0.6))),1])
	print(summary(test_data))
	print(summary(data_rna_sig[, 1]))
	print('dif rna 0')
	print(sum(test_data[!is.na(test_data)]==0)/length(test_data[!is.na(test_data)]))
	print(sum((data_rna_sig[, 1]==0))/dim(data_rna_sig)[1])#*(data_rna_sig[, 2]==0)*(data_rna_sig[, 3]==0)*(data_rna_sig[, 4]==0)*(data_rna_sig[, 5]==0)*(data_rna_sig[, 6]==0)*(data_rna_sig[, 7]==0)*(data_rna_sig[, 8]==0)*(data_rna_sig[, 9]==0)*(data_rna_sig[, 10]==0)*(data_rna_sig[, 11]==0)*(data_rna_sig[, 12]==0)))
	png('norm_compare/scatterplot_check.png')
	heatscatter(data_atac_sig_log2[as.logical((data_atac_rna_sig_cor<=(-0.4)) * (data_atac_rna_sig_cor>=(-0.6))), 1], data_rna_sig[as.logical((data_atac_rna_sig_cor<=(-0.4)) * (data_atac_rna_sig_cor>=(-0.6))), 1], pch = 20)
	dev.off()
	png('norm_compare/scatterplot_check_all.png')
	heatscatter(data_atac_sig_log2[, 1], data_rna_sig[, 1], pch = 20)
	dev.off()

	print(dim(data_atac_sig_sf_am_log2_1))
	print(dim(data_rna_sig))
	data_atac_sf_sd_am_rna_sig_cor1 = apply(cbind(data_atac_sig_sf_am_log2_1, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_sf_sd_rna_sig_cor1 = apply(cbind(data_atac_sig_sf_log2_1, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )

	data_atac_sf_sd_am_rna_sig_cor5 = apply(cbind(data_atac_sig_sf_am_log2_5, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_sf_sd_rna_sig_cor5 = apply(cbind(data_atac_sig_sf_log2_5, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )

	data_atac_sf_sd_am_rna_sig_cor10 = apply(cbind(data_atac_sig_sf_am_log2_10, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_sf_sd_rna_sig_cor10 = apply(cbind(data_atac_sig_sf_log2_10, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )

	data_atac_sf_sd_am_rna_sig_cor1000kb = apply(cbind(data_atac_sig_sf_am_log2_1000kb, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_sf_sd_rna_sig_cor1000kb = apply(cbind(data_atac_sig_sf_log2_1000kb, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )


	data_atac_sig_log2_tr_1mb_cor1000kb_shuffle = apply(cbind( data_atac_sig_log2_tr_1mb[sample(nrow(data_atac_sig_log2_tr_1mb)),], data_rna_sig ), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_sig_log2_tr_1mb_cor1000kb_shuffle = apply(cbind(data_atac_sig_log2_tr_1mb[sample(nrow(data_atac_sig_log2_tr_1mb)),], data_rna_sig[sample(nrow(data_rna_sig)),]), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )

	data_atac_sig_log2_1mb_shuffle = apply(cbind(data_atac_sig_log2_1mb[sample(nrow(data_atac_sig_log2_1mb)),], data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )

	data_atac_sf_sd_am_rna_sig_cor1000kb_shuffle = apply(cbind(data_atac_sig_sf_am_log2_1000kb[sample(dim(data_atac_sig_sf_am_log2_1000kb)[1]),], data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_sf_sd_rna_sig_cor1000kb_shuffle = apply(cbind(data_atac_sig_sf_log2_1000kb[sample(dim(data_atac_sig_sf_log2_1000kb)[1]),], data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )


	data_rna = cbind(data_rna, data_atac_rna_sig_cor, data_atac_rna_sig_cor_tr, data_atac_sf_sd_am_rna_sig_cor10, data_atac_sf_sd_rna_sig_cor10)


	print('kl.dist')
	dist_tss_sig_cor_10kb = as.matrix(cbind( density(data_atac_rna_sig_cor, na.rm = T)$x, density(data_atac_rna_sig_cor, na.rm = T)$y ))
	dist_tss_sig_cor_1000kb = as.matrix(cbind( density(data_atac_rna_sig_cor_1mb, na.rm = T)$x, density(data_atac_rna_sig_cor_1mb, na.rm = T)$y ))
	print('raw')
	print(kl.dist(dist_tss_sig_cor_10kb, dist_tss_sig_cor_1000kb))

	print('kl.dist')
	dist_tss_sig_cor_10kb = as.matrix(cbind( density(data_atac_rna_sig_cor_tr, na.rm = T)$x, density(data_atac_rna_sig_cor_tr, na.rm = T)$y ))
	dist_tss_sig_cor_1000kb = as.matrix(cbind( density(data_atac_rna_sig_cor_tr_1mb, na.rm = T)$x, density(data_atac_rna_sig_cor_tr_1mb, na.rm = T)$y ))
	print('TR')
	print(kl.dist(dist_tss_sig_cor_10kb, dist_tss_sig_cor_1000kb))

	print('kl.dist')
	dist_tss_sig_cor_10kb = as.matrix(cbind( density(data_atac_sf_sd_am_rna_sig_cor10, na.rm = T)$x, density(data_atac_sf_sd_am_rna_sig_cor10, na.rm = T)$y ))
	dist_tss_sig_cor_1000kb = as.matrix(cbind( density(data_atac_sf_sd_am_rna_sig_cor1000kb, na.rm = T)$x, density(data_atac_sf_sd_am_rna_sig_cor1000kb, na.rm = T)$y ))
	print('all mean')
	print(kl.dist(dist_tss_sig_cor_10kb, dist_tss_sig_cor_1000kb))

	print('kl.dist')
	dist_tss_sig_cor_10kb = as.matrix(cbind( density(data_atac_sf_sd_rna_sig_cor10, na.rm = T)$x, density(data_atac_sf_sd_rna_sig_cor10, na.rm = T)$y ))
	dist_tss_sig_cor_1000kb = as.matrix(cbind( density(data_atac_sf_sd_rna_sig_cor1000kb, na.rm = T)$x, density(data_atac_sf_sd_rna_sig_cor1000kb, na.rm = T)$y ))
	print('ncis sig sig')
	print(kl.dist(dist_tss_sig_cor_10kb, dist_tss_sig_cor_1000kb))


	print('correlation summary: ')
	print(length(data_atac_rna_sig_cor))
	print(summary((data_atac_rna_sig_cor)))
	print(summary((data_atac_rna_sig_cor_tr)))

	cor_summary_matrix = c()
	print('1kb')
	print(summary((data_atac_rna_sig_cor_1kb)))
	print(summary((data_atac_rna_sig_cor_tr_1kb)))
	print(summary((data_atac_sf_sd_am_rna_sig_cor1)))
	print(summary((data_atac_sf_sd_rna_sig_cor1)))

	print('5kb')
	print(summary((data_atac_rna_sig_cor_5kb)))
	print(summary((data_atac_rna_sig_cor_tr_5kb)))
	print(summary((data_atac_sf_sd_am_rna_sig_cor5)))
	print(summary((data_atac_sf_sd_rna_sig_cor5)))

	print('10kb')
	print(summary((data_atac_rna_sig_cor)))
	print(summary((data_atac_rna_sig_cor_tr)))
	print(summary((data_atac_sf_sd_am_rna_sig_cor10)))
	print(summary((data_atac_sf_sd_rna_sig_cor10)))

	print('1000kb')
	print(summary((data_atac_rna_sig_cor_1mb)))
	print(summary((data_atac_rna_sig_cor_tr_1mb)))
	print(summary((data_atac_sf_sd_am_rna_sig_cor1000kb)))
	print(summary((data_atac_sf_sd_rna_sig_cor1000kb)))

	print('shuffle')
	print(summary((data_atac_sig_log2_1mb_shuffle)))
	print(summary((data_atac_sig_log2_tr_1mb_cor1000kb_shuffle)))
	print(summary((data_atac_sf_sd_am_rna_sig_cor1000kb_shuffle)))
	print(summary((data_atac_sf_sd_rna_sig_cor1000kb_shuffle)))


	cor_summary_matrix = cbind(summary((data_atac_rna_sig_cor)), summary((data_atac_rna_sig_cor_tr)), summary((data_atac_sf_sd_am_rna_sig_cor10)), summary((data_atac_sf_sd_rna_sig_cor10)), summary((data_atac_rna_sig_cor_1mb)), summary((data_atac_rna_sig_cor_tr_1mb)), summary((data_atac_sf_sd_am_rna_sig_cor1000kb)), summary((data_atac_sf_sd_rna_sig_cor1000kb)))
	write.table(cor_summary_matrix, paste('norm_compare/data_rna.cor.', method,'.txt', sep=''), quote=F, sep='\t', row.names = FALSE, col.names = FALSE)

	png(paste('norm_compare/', method, '_corr.png', sep=''), width = 8, height = 8, units = 'in', res = 300)
	par(mfrow=c(2,2))
	plot(density((data_atac_rna_sig_cor[!is.na(data_atac_rna_sig_cor)])), ylim=c(0,2.5), main='original atac signal TSS 10 kb', col='blue')
	lines(density((data_atac_rna_sig_cor_1kb[!is.na(data_atac_rna_sig_cor_1kb)])), ylim=c(0,2.5), col='yellow')
	lines(density((data_atac_rna_sig_cor_5kb[!is.na(data_atac_rna_sig_cor_5kb)])), ylim=c(0,2.5), col='green')
	lines(density((data_atac_rna_sig_cor_1mb[!is.na(data_atac_rna_sig_cor_1mb)])), ylim=c(0,2.5), col='gray')
	lines(density((data_atac_sig_log2_1mb_shuffle[!is.na(data_atac_sig_log2_1mb_shuffle)])), ylim=c(0,2.5), col='black')

	abline(v=0, lty=2)

	plot(density((data_atac_rna_sig_cor_tr[!is.na(data_atac_rna_sig_cor_tr)])), ylim=c(0,2.5), main='TR norm atac signal TSS 10 kb', col='blue')
	lines(density((data_atac_rna_sig_cor_tr_1kb[!is.na(data_atac_rna_sig_cor_tr_1kb)])), ylim=c(0,2.5), col='yellow')
	lines(density((data_atac_rna_sig_cor_tr_5kb[!is.na(data_atac_rna_sig_cor_tr_5kb)])), ylim=c(0,2.5), col='green')
	lines(density((data_atac_rna_sig_cor_tr_1mb[!is.na(data_atac_rna_sig_cor_tr_1mb)])), ylim=c(0,2.5), col='gray')
	lines(density((data_atac_sig_log2_tr_1mb_cor1000kb_shuffle[!is.na(data_atac_sig_log2_tr_1mb_cor1000kb_shuffle)])), ylim=c(0,2.5), col='black')

	abline(v=0, lty=2)

	plot(density((data_atac_sf_sd_am_rna_sig_cor1[!is.na(data_atac_sf_sd_am_rna_sig_cor1)])), ylim=c(0,2.5), col='yellow', main='Total mean norm atac signal TSS 10 kb')
	lines(density((data_atac_sf_sd_am_rna_sig_cor5[!is.na(data_atac_sf_sd_am_rna_sig_cor5)])), ylim=c(0,2.5), col='green')
	lines(density((data_atac_sf_sd_am_rna_sig_cor10[!is.na(data_atac_sf_sd_am_rna_sig_cor10)])), ylim=c(0,2.5), col='blue')
	lines(density((data_atac_sf_sd_am_rna_sig_cor1000kb[!is.na(data_atac_sf_sd_am_rna_sig_cor1000kb)])), ylim=c(0,2.5), col='gray')
	lines(density((data_atac_sf_sd_am_rna_sig_cor1000kb_shuffle[!is.na(data_atac_sf_sd_am_rna_sig_cor1000kb_shuffle)])), ylim=c(0,2.5), col='black')

	abline(v=0, lty=2)

	plot(density((data_atac_sf_sd_rna_sig_cor1[!is.na(data_atac_sf_sd_rna_sig_cor1)])), ylim=c(0,2.5), col='yellow', main='Both signal mean norm atac signal TSS 10 kb')
	lines(density((data_atac_sf_sd_rna_sig_cor5[!is.na(data_atac_sf_sd_rna_sig_cor5)])), ylim=c(0,2.5), col='green')
	lines(density((data_atac_sf_sd_rna_sig_cor10[!is.na(data_atac_sf_sd_rna_sig_cor10)])), ylim=c(0,2.5), col='blue')
	lines(density((data_atac_sf_sd_rna_sig_cor1000kb[!is.na(data_atac_sf_sd_rna_sig_cor1000kb)])), ylim=c(0,2.5), col='gray')
	lines(density((data_atac_sf_sd_rna_sig_cor1000kb_shuffle[!is.na(data_atac_sf_sd_rna_sig_cor1000kb_shuffle)])), ylim=c(0,2.5), col='black')

	abline(v=0, lty=2)

	dev.off()

	print(colnames(data_atac_sig_log2))
	print(colnames(data_rna_sig))
}


write.table(data_rna, 'norm_compare/data_rna.cor.txt', quote=F, sep='\t', row.names = FALSE, col.names = FALSE)


data_atac_sig_var = apply(data_atac_sig_log2, 1, var)
data_atac_sig_var_tr = apply(data_atac_sig_log2_tr, 1, var)
data_atac_sig_sf_am_log2_10_var = apply(data_atac_sig_sf_am_log2_10, 1, var)
data_atac_sig_sf_log2_10_var = apply(data_atac_sig_sf_log2_10, 1, var)
data_rna_sig_var = apply(data_rna_sig, 1, var)

#print(dim(data_rna_sig))
#print(length(data_atac_sig_log2))
#print(length(data_rna_sig_var))
png('norm_compare/variance.png')
par(mfrow=c(2,2))
rna_log=log(data_rna_sig_var)
rm_rna=as.logical( (!is.na(rna_log)) * (is.finite(rna_log)) )

summary(rna_log[rm_rna])
heatscatter((data_atac_sig_var), rna_log, pch = 20, xlim=c(0, 2), ylim=c(-5, 5))
abline(lm(rna_log[rm_rna]~data_atac_sig_var[rm_rna]), col="red") # regression line (y~x) 
lines(lowess(data_atac_sig_var[rm_rna],rna_log[rm_rna]), col="blue") # lowess line (x,y)

heatscatter((data_atac_sig_var_tr), rna_log, pch = 20, xlim=c(0, 2), ylim=c(-5, 5))
abline(lm(rna_log[rm_rna]~data_atac_sig_var_tr[rm_rna]), col="red") # regression line (y~x) 
lines(lowess(data_atac_sig_var_tr[rm_rna],rna_log[rm_rna]), col="blue") # lowess line (x,y)

heatscatter((data_atac_sig_sf_am_log2_10_var), rna_log, pch = 20, xlim=c(0, 2), ylim=c(-5, 5))
abline(lm(rna_log[rm_rna]~data_atac_sig_sf_am_log2_10_var[rm_rna]), col="red") # regression line (y~x) 
lines(lowess(data_atac_sig_sf_am_log2_10_var[rm_rna],rna_log[rm_rna]), col="blue") # lowess line (x,y)

heatscatter((data_atac_sig_sf_log2_10_var), rna_log, pch = 20, xlim=c(0, 2), ylim=c(-5, 5))
abline(lm(rna_log[rm_rna]~data_atac_sig_sf_log2_10_var[rm_rna]), col="red") # regression line (y~x) 
lines(lowess(data_atac_sig_sf_log2_10_var[rm_rna],rna_log[rm_rna]), col="blue") # lowess line (x,y)

dev.off()

png('norm_compare/variance_high.png')
par(mfrow=c(2,2))
heatscatter((data_atac_sig_var)[log(data_rna_sig_var)>-10], log(data_rna_sig_var)[log(data_rna_sig_var)>-10], pch = 20, xlim=c(0, 2), ylim=c(-5, 5))
heatscatter((data_atac_sig_var_tr)[log(data_rna_sig_var)>-10], log(data_rna_sig_var)[log(data_rna_sig_var)>-10], pch = 20, xlim=c(0, 2), ylim=c(-5, 5))
heatscatter((data_atac_sig_sf_am_log2_10_var)[log(data_rna_sig_var)>-10], log(data_rna_sig_var)[log(data_rna_sig_var)>-10], pch = 20, xlim=c(0, 2), ylim=c(-5, 5))
heatscatter((data_atac_sig_sf_log2_10_var)[log(data_rna_sig_var)>-10], log(data_rna_sig_var)[log(data_rna_sig_var)>-10], pch = 20, xlim=c(0, 2), ylim=c(-5, 5))
dev.off()

print('start scatterplot!!!!!!')
for (i in c(1:12)){
	png(paste('norm_compare/log2/scatterplot_log2_',toString(i),'.png',sep=''))
	par(mfrow=c(2,2))
	accuracy = sum( ((data_atac_sig_log2[,i])-data_atac_sig_log2[,2])*(data_rna_sig[,i]-data_rna_sig[,2])>0 ) / length(data_atac_sig_log2[,2])
	heatscatter(data_atac_sig_log2[,i]-data_atac_sig_log2[,2], data_rna_sig[,i]-data_rna_sig[,2], pch = 20, main=paste(colnames(data_atac_sig_log2)[i], toString(accuracy)) )
	abline(0,1)
	abline(h=0)
	abline(v=0)
	accuracy = sum( ((data_atac_sig_log2_tr[,i])-data_atac_sig_log2_tr[,2])*(data_rna_sig[,i]-data_rna_sig[,2])>0 ) / length(data_atac_sig_log2[,2])
	heatscatter(data_atac_sig_log2_tr[,i]-data_atac_sig_log2_tr[,2], data_rna_sig[,i]-data_rna_sig[,2], pch = 20, main=paste(colnames(data_atac_sig_log2)[i], toString(accuracy)) )
	abline(0,1)
	abline(h=0)
	abline(v=0)
	accuracy = sum( ((data_atac_sig_sf_am_log2_10[,i])-data_atac_sig_sf_am_log2_10[,2])*(data_rna_sig[,i]-data_rna_sig[,2])>0 ) / length(data_atac_sig_log2[,2])
	heatscatter((data_atac_sig_sf_am_log2_10[,i])-data_atac_sig_sf_am_log2_10[,2], data_rna_sig[,i]-data_rna_sig[,2], pch = 20, main=paste(colnames(data_atac_sig_log2)[i], toString(accuracy)) )
	abline(0,1)
	abline(h=0)
	abline(v=0)
	accuracy = sum( ((data_atac_sig_sf_log2_10[,i])-data_atac_sig_sf_log2_10[,2])*(data_rna_sig[,i]-data_rna_sig[,2])>0 ) / length(data_atac_sig_log2[,2])
	heatscatter((data_atac_sig_sf_log2_10[,i])-data_atac_sig_sf_log2_10[,2], data_rna_sig[,i]-data_rna_sig[,2], pch = 20, main=paste(colnames(data_atac_sig_log2)[i], toString(accuracy)) )
	abline(0,1)
	abline(h=0)
	abline(v=0)
	dev.off()
}


print('start scatterplot!!!!!!')
for (i in c(1:12)){
	png(paste('norm_compare/log2/scatterplot_log2_bg_',toString(i),'.png',sep=''))
	par(mfrow=c(2,2))
	accuracy = sum( ((data_atac_sig_log2[,i])-data_atac_sig_log2[,2])*(data_rna_sig[,i]-data_rna_sig[,2])>0 ) / length(data_atac_sig_log2[,2])
	heatscatter(data_atac_sig_log2[,i]-data_atac_sig_log2[,2], data_rna_sig[,i]-data_rna_sig[,2], pch = 20, main=paste(colnames(data_atac_sig_log2)[i], toString(accuracy)) )
	abline(0,1)
	abline(h=0)
	abline(v=0)
	accuracy = sum( ((data_atac_sig_log2_tr_1mb[,i])-data_atac_sig_log2_tr_1mb[,2])*(data_rna_sig[,i]-data_rna_sig[,2])>0 ) / length(data_atac_sig_log2[,2])
	heatscatter(data_atac_sig_log2_tr_1mb[,i]-data_atac_sig_log2_tr_1mb[,2], data_rna_sig[,i]-data_rna_sig[,2], pch = 20, main=paste(colnames(data_atac_sig_log2)[i], toString(accuracy)) )
	abline(0,1)
	abline(h=0)
	abline(v=0)
	accuracy = sum( ((data_atac_sig_sf_am_log2_1000kb[,i])-data_atac_sig_sf_am_log2_1000kb[,2])*(data_rna_sig[,i]-data_rna_sig[,2])>0 ) / length(data_atac_sig_log2[,2])
	heatscatter((data_atac_sig_sf_am_log2_1000kb[,i])-data_atac_sig_sf_am_log2_1000kb[,2], data_rna_sig[,i]-data_rna_sig[,2], pch = 20, main=paste(colnames(data_atac_sig_log2)[i], toString(accuracy)) )
	abline(0,1)
	abline(h=0)
	abline(v=0)
	accuracy = sum( ((data_atac_sig_sf_log2_1000kb[,i])-data_atac_sig_sf_log2_1000kb[,2])*(data_rna_sig[,i]-data_rna_sig[,2])>0 ) / length(data_atac_sig_log2[,2])
	heatscatter((data_atac_sig_sf_log2_1000kb[,i])-data_atac_sig_sf_log2_1000kb[,2], data_rna_sig[,i]-data_rna_sig[,2], pch = 20, main=paste(colnames(data_atac_sig_log2)[i], toString(accuracy)) )
	abline(0,1)
	abline(h=0)
	abline(v=0)
	dev.off()
}



for (i in c(1:12)){
	png(paste('norm_compare/linear/scatterplot_log2_',toString(i),'.png',sep=''))
	par(mfrow=c(2,2))
	heatscatter(data_atac_sig_log2[,i], data_rna_sig[,i], pch = 20, ylim=c(-5,15), xlim=c(-5,15), main=colnames(data_atac_sig_log2)[i])
	abline(0,1)
	heatscatter(data_atac_sig_log2_tr[,i], data_rna_sig[,i], pch = 20, ylim=c(-5,15), xlim=c(-5,15), main=colnames(data_atac_sig_log2)[i])
	abline(0,1)
	heatscatter((data_atac_sig_sf_am_log2_10[,i]), data_rna_sig[,i], pch = 20, ylim=c(-5,15), xlim=c(-5,15), main=colnames(data_atac_sig_log2)[i])
	abline(0,1)
	heatscatter((data_atac_sig_sf_log2_10[,i]), data_rna_sig[,i], pch = 20, ylim=c(-5,15), xlim=c(-5,15), main=colnames(data_atac_sig_log2)[i])
	abline(0,1)
	dev.off()
}


png('norm_compare/rna_scatterplot1.png')
par(mfrow=c(2,2))
heatscatter(data_rna_sig[,2], data_rna_sig[,1], pch = 20, xlim=c(-5, 20), ylim=c(-5, 20), main='lsk vs cmp')
abline(0,1, col='black')
heatscatter(data_rna_sig[,2], data_rna_sig[,2], pch = 20, xlim=c(-5, 20), ylim=c(-5, 20), main='cmp vs cmp')
abline(0,1, col='black')
heatscatter(data_rna_sig[,2], data_rna_sig[,3], pch = 20, xlim=c(-5, 20), ylim=c(-5, 20), main='gmp vs cmp')
abline(0,1, col='black')
heatscatter(data_rna_sig[,2], data_rna_sig[,4], pch = 20, xlim=c(-5, 20), ylim=c(-5, 20), main='mep vs cmp')
abline(0,1, col='black')
dev.off()
png('norm_compare/rna_scatterplot2.png')
par(mfrow=c(2,2))
heatscatter(data_rna_sig[,2], data_rna_sig[,5], pch = 20, xlim=c(-5, 20), ylim=c(-5, 20), main='cfue vs cmp')
abline(0,1, col='black')
heatscatter(data_rna_sig[,2], data_rna_sig[,6], pch = 20, xlim=c(-5, 20), ylim=c(-5, 20), main='ery vs cmp')
abline(0,1, col='black')
heatscatter(data_rna_sig[,2], data_rna_sig[,7], pch = 20, xlim=c(-5, 20), ylim=c(-5, 20), main='cfumk vs cmp')
abline(0,1, col='black')
heatscatter(data_rna_sig[,2], data_rna_sig[,8], pch = 20, xlim=c(-5, 20), ylim=c(-5, 20), main='meg vs cmp')
abline(0,1, col='black')
dev.off()
png('norm_compare/rna_scatterplot3.png')
par(mfrow=c(2,2))
heatscatter(data_rna_sig[,2], data_rna_sig[,9], pch = 20, xlim=c(-5, 20), ylim=c(-5, 20), main='mono vs cmp')
abline(0,1, col='black')
heatscatter(data_rna_sig[,2], data_rna_sig[,10], pch = 20, xlim=c(-5, 20), ylim=c(-5, 20), main='neu vs cmp')
abline(0,1, col='black')
heatscatter(data_rna_sig[,2], data_rna_sig[,11], pch = 20, xlim=c(-5, 20), ylim=c(-5, 20), main='g1e vs cmp')
abline(0,1, col='black')
heatscatter(data_rna_sig[,2], data_rna_sig[,12], pch = 20, xlim=c(-5, 20), ylim=c(-5, 20), main='er4 vs cmp')
abline(0,1, col='black')
dev.off()

print('start scatterplot DONE!!!!!!')

