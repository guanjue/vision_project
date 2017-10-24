library(LSD)
library(seewave)

data_atac = read.table('gene_atac/gencode_pc_sort.TSSexp10kb.atac.celltype.matched.txt', header = T, sep='\t')
data_atac_1kb = read.table('gene_atac/gencode_pc_sort.TSSexp1kb.atac.celltype.matched.txt', header = T, sep='\t')
data_atac_5kb = read.table('gene_atac/gencode_pc_sort.TSSexp5kb.atac.celltype.matched.txt', header = T, sep='\t')
data_atac_1mb = read.table('gene_atac/gencode_pc_sort.TSSup1000kb.atac.celltype.matched.txt', header = T, sep='\t')
data_atac_tr = read.table('gene_atac/gencode_pc_sort.TSSexp10kb.atac.celltype.matched.TRnormed.txt', header = T, sep='\t')
data_atac1kb = read.table('gene_atac/gencode_pc_sort.TSSexp1kb.atac.celltype.matched.txt', header = T, sep='\t')
data_atac5kb = read.table('gene_atac/gencode_pc_sort.TSSexp5kb.atac.celltype.matched.txt', header = T, sep='\t')
data_atac10kb = read.table('gene_atac/gencode_pc_sort.TSSexp10kb.atac.celltype.matched.txt', header = T, sep='\t')
data_atac110kb = read.table('gene_atac/gencode_pc_sort.TSSup1000kb.atac.celltype.matched.txt', header = T, sep='\t')
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
data_atac_sig110 = data_atac110kb[,c(-1,-2,-3,-4,-5,-6)]

### scale factor norm
data_atac_sig_sf = t(apply(data_atac_sig, 1, FUN=function(x) x*scale_factor))
data_atac_sig_sf_am = t(apply(data_atac_sig, 1, FUN=function(x) x*scale_factor_allmean))

data_atac_sig_sf1 = t(apply(data_atac_sig1, 1, FUN=function(x) x*scale_factor))
data_atac_sig_sf_am1 = t(apply(data_atac_sig1, 1, FUN=function(x) x*scale_factor_allmean))

data_atac_sig_sf5 = t(apply(data_atac_sig5, 1, FUN=function(x) x*scale_factor))
data_atac_sig_sf_am5 = t(apply(data_atac_sig5, 1, FUN=function(x) x*scale_factor_allmean))

data_atac_sig_sf10 = t(apply(data_atac_sig10, 1, FUN=function(x) x*scale_factor))
data_atac_sig_sf_am10 = t(apply(data_atac_sig10, 1, FUN=function(x) x*scale_factor_allmean))

data_atac_sig_sf110 = t(apply(data_atac_sig110, 1, FUN=function(x) x*scale_factor))
data_atac_sig_sf_am110 = t(apply(data_atac_sig110, 1, FUN=function(x) x*scale_factor_allmean))


###### input matrix
### od log2
data_atac_sig_log2 = log2((data_atac_sig+1)/10000*1000)
print(head(data_atac_sig_log2))
used_id = apply((data_atac_sig_log2), 1, FUN=function(x) (min((x)) > log2(0.1)) )
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
data_atac_sig_sf_am_log2_1 = log2((data_atac_sig_sf_am1+1)/1000*1000)[used_id, ]
data_atac_sig_sf_log2_1 = log2((data_atac_sig_sf1+1)/1000*1000)[used_id, ]

data_atac_sig_sf_am_log2_5 = log2((data_atac_sig_sf_am5+1)/5000*1000)[used_id, ]
data_atac_sig_sf_log2_5 = log2((data_atac_sig_sf5+1)/5000*1000)[used_id, ]

data_atac_sig_sf_am_log2_10 = log2((data_atac_sig_sf_am10+1)/10000*1000)[used_id, ]
data_atac_sig_sf_log2_10 = log2((data_atac_sig_sf10+1)/10000*1000)[used_id, ]

data_atac_sig_sf_am_log2_110 = log2((data_atac_sig_sf_am110+1)/10000*1000)[used_id, ]
data_atac_sig_sf_log2_110 = log2((data_atac_sig_sf110+1)/10000*1000)[used_id, ]

data_rna_sig=data_rna_sig[used_id, ]


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

	data_atac_sig_log2_tr_1mb_cor110_shuffle = apply(cbind( data_atac_sig_log2_tr_1mb[sample(nrow(data_atac_sig_log2_tr_1mb)),], data_rna_sig ), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )

	print('kl.dist')
	print(cor(data_atac_rna_sig_cor_tr, data_atac_rna_sig_cor_tr_1mb, use ='complete.obs'))
	print(cor(data_atac_rna_sig_cor_tr, data_atac_sig_log2_tr_1mb_cor110_shuffle, use ='complete.obs'))

	dist_tss = as.matrix(cbind( density(data_atac_rna_sig_cor_tr, na.rm = T)$x, density(data_atac_rna_sig_cor_tr, na.rm = T)$y ))
	dist_tss_5kb = as.matrix(cbind( density(data_atac_rna_sig_cor_tr_5kb, na.rm = T)$x, density(data_atac_rna_sig_cor_tr_5kb, na.rm = T)$y ))
	dist_tss_1mb = as.matrix(cbind( density(data_atac_rna_sig_cor_tr_1mb, na.rm = T)$x, density(data_atac_rna_sig_cor_tr_1mb, na.rm = T)$y ))
	dist_tss_1mb_shuffle = as.matrix(cbind( density(data_atac_sig_log2_tr_1mb_cor110_shuffle, na.rm = T)$x, density(data_atac_sig_log2_tr_1mb_cor110_shuffle, na.rm = T)$y ))

	print('1mb')
	print(kl.dist(dist_tss, dist_tss_1mb))
	print('1mb shuffle')
	print(kl.dist(dist_tss, dist_tss_1mb_shuffle))
	print('5kb')
	print(kl.dist(dist_tss, dist_tss_5kb))
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


	data_atac_sf_sd_am_rna_sig_cor1 = apply(cbind(data_atac_sig_sf_am_log2_1, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_sf_sd_rna_sig_cor1 = apply(cbind(data_atac_sig_sf_log2_1, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )

	data_atac_sf_sd_am_rna_sig_cor5 = apply(cbind(data_atac_sig_sf_am_log2_5, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_sf_sd_rna_sig_cor5 = apply(cbind(data_atac_sig_sf_log2_5, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )

	data_atac_sf_sd_am_rna_sig_cor10 = apply(cbind(data_atac_sig_sf_am_log2_10, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_sf_sd_rna_sig_cor10 = apply(cbind(data_atac_sig_sf_log2_10, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )

	data_atac_sf_sd_am_rna_sig_cor110 = apply(cbind(data_atac_sig_sf_am_log2_110, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_sf_sd_rna_sig_cor110 = apply(cbind(data_atac_sig_sf_log2_110, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )



	data_atac_sig_log2_tr_1mb_cor110_shuffle = apply(cbind( data_atac_sig_log2_tr_1mb[sample(nrow(data_atac_sig_log2_tr_1mb)),], data_rna_sig ), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_sig_log2_tr_1mb_cor110_shuffle = apply(cbind(data_atac_sig_log2_tr_1mb[sample(nrow(data_atac_sig_log2_tr_1mb), 1000),], data_rna_sig[sample(nrow(data_rna_sig), 10000),]), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )

	data_atac_sig_log2_1mb_shuffle = apply(cbind(data_atac_sig_log2_1mb[sample(nrow(data_atac_sig_log2_1mb)),], data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )

	data_atac_sf_sd_am_rna_sig_cor110_shuffle = apply(cbind(data_atac_sig_sf_am_log2_110[sample(nrow(data_atac_sig_sf_am_log2_110)),], data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )
	data_atac_sf_sd_rna_sig_cor110_shuffle = apply(cbind(data_atac_sig_sf_log2_110[sample(nrow(data_atac_sig_sf_log2_110)),], data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method=method) )


	print('correlation summary: ')
	print(length(data_atac_rna_sig_cor))
	print(summary((data_atac_rna_sig_cor)))
	print(summary((data_atac_rna_sig_cor_tr)))

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
	print(summary((data_atac_sf_sd_am_rna_sig_cor110)))
	print(summary((data_atac_sf_sd_rna_sig_cor110)))

	print('shuffle')
	print(summary((data_atac_sig_log2_1mb_shuffle)))
	print(summary((data_atac_sig_log2_tr_1mb_cor110_shuffle)))
	print(summary((data_atac_sf_sd_am_rna_sig_cor110_shuffle)))
	print(summary((data_atac_sf_sd_rna_sig_cor110_shuffle)))

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
	lines(density((data_atac_sig_log2_tr_1mb_cor110_shuffle[!is.na(data_atac_sig_log2_tr_1mb_cor110_shuffle)])), ylim=c(0,2.5), col='black')

	abline(v=0, lty=2)

	plot(density((data_atac_sf_sd_am_rna_sig_cor1[!is.na(data_atac_sf_sd_am_rna_sig_cor1)])), ylim=c(0,2.5), col='yellow', main='Total mean norm atac signal TSS 10 kb')
	lines(density((data_atac_sf_sd_am_rna_sig_cor5[!is.na(data_atac_sf_sd_am_rna_sig_cor5)])), ylim=c(0,2.5), col='green')
	lines(density((data_atac_sf_sd_am_rna_sig_cor10[!is.na(data_atac_sf_sd_am_rna_sig_cor10)])), ylim=c(0,2.5), col='blue')
	lines(density((data_atac_sf_sd_am_rna_sig_cor110[!is.na(data_atac_sf_sd_am_rna_sig_cor110)])), ylim=c(0,2.5), col='gray')
	lines(density((data_atac_sf_sd_am_rna_sig_cor110_shuffle[!is.na(data_atac_sf_sd_am_rna_sig_cor110_shuffle)])), ylim=c(0,2.5), col='black')

	abline(v=0, lty=2)

	plot(density((data_atac_sf_sd_rna_sig_cor1[!is.na(data_atac_sf_sd_rna_sig_cor1)])), ylim=c(0,2.5), col='yellow', main='Both signal mean norm atac signal TSS 10 kb')
	lines(density((data_atac_sf_sd_rna_sig_cor5[!is.na(data_atac_sf_sd_rna_sig_cor5)])), ylim=c(0,2.5), col='green')
	lines(density((data_atac_sf_sd_rna_sig_cor10[!is.na(data_atac_sf_sd_rna_sig_cor10)])), ylim=c(0,2.5), col='blue')
	lines(density((data_atac_sf_sd_rna_sig_cor110[!is.na(data_atac_sf_sd_rna_sig_cor110)])), ylim=c(0,2.5), col='gray')
	lines(density((data_atac_sf_sd_rna_sig_cor110_shuffle[!is.na(data_atac_sf_sd_rna_sig_cor110_shuffle)])), ylim=c(0,2.5), col='black')

	abline(v=0, lty=2)

	dev.off()

	print(colnames(data_atac_sig_log2))
	print(colnames(data_rna_sig))
}





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
rna_log=(data_rna_sig_var)
rm_rna=as.logical( (!is.na(rna_log)) * (is.finite(rna_log)) )

summary(rna_log[rm_rna])
heatscatter((data_atac_sig_var), rna_log, pch = 20, xlim=c(0, 3), ylim=c(-10, 20))
abline(lm(rna_log[rm_rna]~data_atac_sig_var[rm_rna]), col="red") # regression line (y~x) 
lines(lowess(data_atac_sig_var[rm_rna],rna_log[rm_rna]), col="blue") # lowess line (x,y)

heatscatter((data_atac_sig_var_tr), rna_log, pch = 20, xlim=c(0, 3), ylim=c(-10, 20))
abline(lm(rna_log[rm_rna]~data_atac_sig_var_tr[rm_rna]), col="red") # regression line (y~x) 
lines(lowess(data_atac_sig_var_tr[rm_rna],rna_log[rm_rna]), col="blue") # lowess line (x,y)

heatscatter((data_atac_sig_sf_am_log2_10_var), rna_log, pch = 20, xlim=c(0, 3), ylim=c(-10, 20))
abline(lm(rna_log[rm_rna]~data_atac_sig_sf_am_log2_10_var[rm_rna]), col="red") # regression line (y~x) 
lines(lowess(data_atac_sig_sf_am_log2_10_var[rm_rna],rna_log[rm_rna]), col="blue") # lowess line (x,y)

heatscatter((data_atac_sig_sf_log2_10_var), rna_log, pch = 20, xlim=c(0, 3), ylim=c(-10, 20))
abline(lm(rna_log[rm_rna]~data_atac_sig_sf_log2_10_var[rm_rna]), col="red") # regression line (y~x) 
lines(lowess(data_atac_sig_sf_log2_10_var[rm_rna],rna_log[rm_rna]), col="blue") # lowess line (x,y)

dev.off()

png('norm_compare/variance_high.png')
par(mfrow=c(2,2))
heatscatter((data_atac_sig_var)[log(data_rna_sig_var)>-10], log(data_rna_sig_var)[log(data_rna_sig_var)>-10], pch = 20, xlim=c(0, 3), ylim=c(-10, 5))
heatscatter((data_atac_sig_var_tr)[log(data_rna_sig_var)>-10], log(data_rna_sig_var)[log(data_rna_sig_var)>-10], pch = 20, xlim=c(0, 3), ylim=c(-10, 5))
heatscatter((data_atac_sig_sf_am_log2_10_var)[log(data_rna_sig_var)>-10], log(data_rna_sig_var)[log(data_rna_sig_var)>-10], pch = 20, xlim=c(0, 3), ylim=c(-10, 5))
heatscatter((data_atac_sig_sf_log2_10_var)[log(data_rna_sig_var)>-10], log(data_rna_sig_var)[log(data_rna_sig_var)>-10], pch = 20, xlim=c(0, 3), ylim=c(-10, 5))
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
	accuracy = sum( ((data_atac_sig_sf_am_log2_110[,i])-data_atac_sig_sf_am_log2_110[,2])*(data_rna_sig[,i]-data_rna_sig[,2])>0 ) / length(data_atac_sig_log2[,2])
	heatscatter((data_atac_sig_sf_am_log2_110[,i])-data_atac_sig_sf_am_log2_110[,2], data_rna_sig[,i]-data_rna_sig[,2], pch = 20, main=paste(colnames(data_atac_sig_log2)[i], toString(accuracy)) )
	abline(0,1)
	abline(h=0)
	abline(v=0)
	accuracy = sum( ((data_atac_sig_sf_log2_110[,i])-data_atac_sig_sf_log2_110[,2])*(data_rna_sig[,i]-data_rna_sig[,2])>0 ) / length(data_atac_sig_log2[,2])
	heatscatter((data_atac_sig_sf_log2_110[,i])-data_atac_sig_sf_log2_110[,2], data_rna_sig[,i]-data_rna_sig[,2], pch = 20, main=paste(colnames(data_atac_sig_log2)[i], toString(accuracy)) )
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
'''
data_atac_rna_sig_cor = apply(cbind(data_atac_sig_log2_sd_norm, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method='pearson') )
data_atac_sf_sd_am_rna_sig_cor = apply(cbind(data_atac_sig_sf_am_log2_sd_norm, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method='pearson') )
data_atac_sf_sd_rna_sig_cor = apply(cbind(data_atac_sig_sf_log2_sd_norm, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method='pearson') )



summary((data_atac_rna_sig_cor))
summary((data_atac_sf_sd_am_rna_sig_cor))
summary((data_atac_sf_sd_rna_sig_cor))



par(mfrow=c(2,2))
plot(density((data_atac_rna_sig_cor[!is.na(data_atac_rna_sig_cor)])), ylim=c(0,2.5))
plot(density((data_atac_sf_sd_am_rna_sig_cor[!is.na(data_atac_sf_sd_am_rna_sig_cor)])), ylim=c(0,2.5))
plot(density((data_atac_sf_sd_rna_sig_cor[!is.na(data_atac_sf_sd_rna_sig_cor)])), ylim=c(0,2.5))























data_atac_sig_log2_cmp = data_atac_sig_log2_sd_norm[,2]
data_atac_sig_log2_gmp = data_atac_sig_log2_sd_norm[,3]
data_atac_sig_log2_ery = data_atac_sig_log2_sd_norm[,6]
data_atac_sig_log2_g1e = data_atac_sig_log2_sd_norm[,11]

data_atac_sig_scale_log2_cmp = data_atac_sig_sf_am_log2_sd_norm[,2]
data_atac_sig_scale_log2_gmp = data_atac_sig_sf_am_log2_sd_norm[,3]
data_atac_sig_scale_log2_ery = data_atac_sig_sf_am_log2_sd_norm[,6]
data_atac_sig_scale_log2_g1e = data_atac_sig_sf_am_log2_sd_norm[,11]

data_atac_sig_sf_sd_log2_cmp = data_atac_sig_sf_log2_sd_norm[,2]
data_atac_sig_sf_sd_log2_gmp = data_atac_sig_sf_log2_sd_norm[,3]
data_atac_sig_sf_sd_log2_ery = data_atac_sig_sf_log2_sd_norm[,6]
data_atac_sig_sf_sd_log2_g1e = data_atac_sig_sf_log2_sd_norm[,11]

data_rna_sig_cmp = data_rna_sig[,2]
data_rna_sig_gmp = data_rna_sig[,3]
data_rna_sig_ery = data_rna_sig[,6]
data_rna_sig_g1e = data_rna_sig[,11]


data_atac_sig_var = apply(data_atac_sig_log2_sd_norm, 1, var)
data_atac_sig_scale_log2_var = apply(data_atac_sig_sf_am_log2_sd_norm, 1, var)
data_atac_sig_sf_sd_log2_var = apply(data_atac_sig_sf_log2_sd_norm, 1, var)
data_rna_sig_var = apply(data_rna_sig, 1, var)


par(mfrow=c(3,3))
heatscatter(log(data_atac_sig_var), log(data_rna_sig_var), pch = 20, xlim=c(-8, 2))
heatscatter(log(data_atac_sig_scale_log2_var), log(data_rna_sig_var), pch = 20, xlim=c(-8, 2))
heatscatter(log(data_atac_sig_sf_sd_log2_var), log(data_rna_sig_var), pch = 20, xlim=c(-8, 2), ylim=c(-8, 2))
heatscatter(log(data_atac_sig_sf_sd_log2_var), log(data_atac_sig_var), pch = 20, xlim=c(-8, 2), ylim=c(-8, 2))
heatscatter(log(data_atac_sig_scale_log2_var), log(data_atac_sig_var), pch = 20, xlim=c(-8, 2), ylim=c(-8, 2))



par(mfrow=c(2,2))
heatscatter(data_atac_sig_log2_gmp, data_rna_sig_gmp, pch = 20)
abline(0,1)
heatscatter(data_atac_sig_scale_log2_gmp, data_rna_sig_gmp, pch = 20)
abline(0,1)
heatscatter(data_atac_sig_sf_sd_log2_gmp, data_rna_sig_gmp, pch = 20)
abline(0,1)

par(mfrow=c(3,3))
heatscatter(data_atac_sig_sf_sd_log2_gmp, data_atac_sig_scale_log2_gmp, pch = 20)
abline(0,1)

heatscatter(data_atac_sig_sf_sd_log2_gmp, data_atac_sig_log2_gmp, pch = 20)
abline(0,1)


heatscatter(data_atac_sig_log2_gmp, data_atac_sig_scale_log2_gmp, pch = 20)
abline(0,1)



heatscatter(data_atac_sig_sf_sd_log2_cmp, data_atac_sig_sf_sd_log2_gmp, pch = 20)
abline(0,1)

heatscatter(data_atac_sig_scale_log2_cmp, data_atac_sig_scale_log2_gmp, pch = 20)
abline(0,1)

heatscatter(data_atac_sig_log2_cmp, data_atac_sig_log2_gmp, pch = 20)
abline(0,1)


par(mfrow=c(3,3))
heatscatter(data_atac_sig_sf_sd_log2_cmp, data_atac_sig_sf_sd_log2_ery, pch = 20)
df = as.data.frame(cbind(data_atac_sig_sf_sd_log2_ery, data_atac_sig_sf_sd_log2_cmp))
lm_fit = lm(data_atac_sig_sf_sd_log2_ery~data_atac_sig_sf_sd_log2_cmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

heatscatter(data_atac_sig_scale_log2_cmp, data_atac_sig_scale_log2_ery, pch = 20)
df = as.data.frame(cbind(data_atac_sig_scale_log2_ery, data_atac_sig_scale_log2_cmp))
lm_fit = lm(data_atac_sig_scale_log2_ery~data_atac_sig_scale_log2_cmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

heatscatter(data_atac_sig_log2_cmp, data_atac_sig_log2_ery, pch = 20)
df = as.data.frame(cbind(data_atac_sig_log2_ery, data_atac_sig_log2_cmp))
lm_fit = lm(data_atac_sig_log2_ery~data_atac_sig_log2_cmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)


heatscatter(data_atac_sig_sf_sd_log2_cmp, data_atac_sig_sf_sd_log2_gmp, pch = 20)
df = as.data.frame(cbind(data_atac_sig_sf_sd_log2_gmp, data_atac_sig_sf_sd_log2_cmp))
lm_fit = lm(data_atac_sig_sf_sd_log2_gmp~data_atac_sig_sf_sd_log2_cmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

heatscatter(data_atac_sig_scale_log2_cmp, data_atac_sig_scale_log2_gmp, pch = 20)
df = as.data.frame(cbind(data_atac_sig_scale_log2_gmp, data_atac_sig_scale_log2_cmp))
lm_fit = lm(data_atac_sig_scale_log2_gmp~data_atac_sig_scale_log2_cmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

heatscatter(data_atac_sig_log2_cmp, data_atac_sig_log2_gmp, pch = 20)
df = as.data.frame(cbind(data_atac_sig_log2_gmp, data_atac_sig_log2_cmp))
lm_fit = lm(data_atac_sig_log2_gmp~data_atac_sig_log2_cmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)



heatscatter(data_atac_sig_sf_sd_log2_cmp, data_atac_sig_sf_sd_log2_g1e, pch = 20)
df = as.data.frame(cbind(data_atac_sig_sf_sd_log2_g1e, data_atac_sig_sf_sd_log2_cmp))
lm_fit = lm(data_atac_sig_sf_sd_log2_g1e~data_atac_sig_sf_sd_log2_cmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

heatscatter(data_atac_sig_scale_log2_cmp, data_atac_sig_scale_log2_g1e, pch = 20)
df = as.data.frame(cbind(data_atac_sig_scale_log2_g1e, data_atac_sig_scale_log2_cmp))
lm_fit = lm(data_atac_sig_scale_log2_g1e~data_atac_sig_scale_log2_cmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

heatscatter(data_atac_sig_log2_cmp, data_atac_sig_log2_g1e, pch = 20)
df = as.data.frame(cbind(data_atac_sig_log2_g1e, data_atac_sig_log2_cmp))
lm_fit = lm(data_atac_sig_log2_g1e~data_atac_sig_log2_cmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)




par(mfrow=c(2,2))
heatscatter(data_rna_sig_cmp-data_rna_sig_ery, data_atac_sig_log2_cmp-data_atac_sig_log2_ery, pch = 20)
df = as.data.frame(cbind(data_rna_sig_cmp-data_rna_sig_ery, data_atac_sig_log2_cmp-data_atac_sig_log2_ery))
lm_fit = lm(data_atac_sig_log2_cmp-data_atac_sig_log2_ery~data_rna_sig_cmp-data_rna_sig_ery, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

heatscatter(data_rna_sig_cmp-data_rna_sig_ery, data_atac_sig_sf_sd_log2_cmp-data_atac_sig_sf_sd_log2_ery, pch = 20)
df = as.data.frame(cbind(data_rna_sig_cmp-data_rna_sig_ery, data_atac_sig_sf_sd_log2_cmp-data_atac_sig_sf_sd_log2_ery))
lm_fit = lm(data_atac_sig_sf_sd_log2_cmp-data_atac_sig_sf_sd_log2_ery~data_rna_sig_cmp-data_rna_sig_ery, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

heatscatter(data_rna_sig_cmp-data_rna_sig_ery, data_atac_sig_scale_log2_cmp-data_atac_sig_scale_log2_ery, pch = 20)
df = as.data.frame(cbind(data_rna_sig_cmp-data_rna_sig_ery, data_atac_sig_scale_log2_cmp-data_atac_sig_scale_log2_ery))
lm_fit = lm(data_atac_sig_scale_log2_cmp-data_atac_sig_scale_log2_ery~data_rna_sig_cmp-data_rna_sig_ery, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

heatscatter(data_rna_sig_cmp-data_rna_sig_gmp, data_atac_sig_log2_cmp-data_atac_sig_log2_gmp, pch = 20)
df = as.data.frame(cbind(data_rna_sig_cmp-data_rna_sig_gmp, data_atac_sig_log2_cmp-data_atac_sig_log2_gmp))
lm_fit = lm(data_atac_sig_log2_cmp-data_atac_sig_log2_gmp~data_rna_sig_cmp-data_rna_sig_gmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)


heatscatter(data_rna_sig_cmp-data_rna_sig_g1e, data_atac_sig_log2_cmp-data_atac_sig_log2_g1e, pch = 20)
df = as.data.frame(cbind(data_rna_sig_cmp-data_rna_sig_g1e, data_atac_sig_log2_cmp-data_atac_sig_log2_g1e))
lm_fit = lm(data_atac_sig_log2_cmp-data_atac_sig_log2_g1e~data_rna_sig_cmp-data_rna_sig_g1e, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)


heatscatter(data_rna_sig_cmp-data_rna_sig_ery, data_atac_sig_log2_cmp-data_atac_sig_log2_ery, pch = 20)
df = as.data.frame(cbind(data_rna_sig_cmp-data_rna_sig_ery, data_atac_sig_log2_cmp-data_atac_sig_log2_ery))
lm_fit = lm(data_atac_sig_log2_cmp-data_atac_sig_log2_ery~data_rna_sig_cmp-data_rna_sig_ery, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)


###
pdf('scatter.pdf', width = 9, height = 12)
par(mfrow=c(4,3))
heatscatter(data_rna_sig_cmp, data_atac_sig_log2_cmp, pch = 20)
df = as.data.frame(cbind(data_atac_sig_log2_cmp, data_rna_sig_cmp))
lm_fit = lm(data_atac_sig_log2_cmp~data_rna_sig_cmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)


heatscatter(data_rna_sig_cmp, data_atac_sig_scale_log2_cmp, pch = 20)
df = as.data.frame(cbind(data_atac_sig_scale_log2_cmp, data_rna_sig_cmp))
lm_fit = lm(data_atac_sig_scale_log2_cmp~data_rna_sig_cmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

heatscatter(data_rna_sig_cmp, data_atac_sig_sf_sd_log2_cmp, pch = 20)
df = as.data.frame(cbind(data_atac_sig_sf_sd_log2_cmp, data_rna_sig_cmp))
lm_fit = lm(data_atac_sig_sf_sd_log2_cmp~data_rna_sig_cmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)


heatscatter(data_rna_sig_gmp, data_atac_sig_sf_sd_log2_gmp, pch = 20)
df = as.data.frame(cbind(data_atac_sig_sf_sd_log2_gmp, data_rna_sig_gmp))
lm_fit = lm(data_atac_sig_sf_sd_log2_gmp~data_rna_sig_gmp, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)


heatscatter(data_rna_sig_g1e, data_atac_sig_scale_log2_g1e, pch = 20)
df = as.data.frame(cbind(data_atac_sig_scale_log2_g1e, data_rna_sig_g1e))
lm_fit = lm(data_atac_sig_scale_log2_g1e~data_rna_sig_g1e, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

heatscatter(data_rna_sig_g1e, data_atac_sig_sf_sd_log2_g1e, pch = 20)
df = as.data.frame(cbind(data_atac_sig_sf_sd_log2_g1e, data_rna_sig_g1e))
lm_fit = lm(data_atac_sig_sf_sd_log2_g1e~data_rna_sig_g1e, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)



heatscatter(data_rna_sig_ery, data_atac_sig_log2_ery, pch = 20)
df = as.data.frame(cbind(data_atac_sig_log2_ery, data_rna_sig_ery))
lm_fit = lm(data_atac_sig_log2_ery~data_rna_sig_ery, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)


heatscatter(data_rna_sig_ery, data_atac_sig_scale_log2_ery, pch = 20)
df = as.data.frame(cbind(data_atac_sig_scale_log2_ery, data_rna_sig_ery))
lm_fit = lm(data_atac_sig_scale_log2_ery~data_rna_sig_ery, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

heatscatter(data_rna_sig_ery, data_atac_sig_sf_sd_log2_ery, pch = 20)
df = as.data.frame(cbind(data_atac_sig_sf_sd_log2_ery, data_rna_sig_ery))
lm_fit = lm(data_atac_sig_sf_sd_log2_ery~data_rna_sig_ery, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)


heatscatter(data_rna_sig_g1e, data_atac_sig_log2_g1e, pch = 20)
df = as.data.frame(cbind(data_atac_sig_log2_g1e, data_rna_sig_g1e))
lm_fit = lm(data_atac_sig_log2_g1e~data_rna_sig_g1e, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)


heatscatter(data_rna_sig_g1e, data_atac_sig_scale_log2_g1e, pch = 20)
df = as.data.frame(cbind(data_atac_sig_scale_log2_g1e, data_rna_sig_g1e))
lm_fit = lm(data_atac_sig_scale_log2_g1e~data_rna_sig_g1e, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

heatscatter(data_rna_sig_g1e, data_atac_sig_sf_sd_log2_g1e, pch = 20)
df = as.data.frame(cbind(data_atac_sig_sf_sd_log2_g1e, data_rna_sig_g1e))
lm_fit = lm(data_atac_sig_sf_sd_log2_g1e~data_rna_sig_g1e, data=df)
abline(lm_fit$coefficients[1],lm_fit$coefficients[2], col='blue')
abline(0,1)

dev.off()









data_atac_sf_sd_rna_sig_cor_cmp = cor(data_atac_sig_sf_sd_log2_cmp, data_rna_sig_cmp, method='spearman')
summary((data_atac_sf_sd_rna_sig_cor_cmp))
data_atac_sf_sd_rna_sig_cor_cmp = cor(data_atac_sig_sf_sd_log2_cmp, data_rna_sig_cmp, method='pearson')
summary((data_atac_sf_sd_rna_sig_cor_cmp))

data_atac_sf_sd_rna_sig_cor_ery = cor(data_atac_sig_sf_sd_log2_ery, data_rna_sig_ery, method='spearman')
summary((data_atac_sf_sd_rna_sig_cor_ery))
data_atac_sf_sd_rna_sig_cor_ery = cor(data_atac_sig_sf_sd_log2_ery, data_rna_sig_ery, method='pearson')
summary((data_atac_sf_sd_rna_sig_cor_ery))

data_atac_sf_sd_rna_sig_cor_gmp = cor(data_atac_sig_sf_sd_log2_gmp, data_rna_sig_gmp, method='spearman')
summary((data_atac_sf_sd_rna_sig_cor_gmp))
data_atac_sf_sd_rna_sig_cor_gmp = cor(data_atac_sig_sf_sd_log2_gmp, data_rna_sig_gmp, method='pearson')
summary((data_atac_sf_sd_rna_sig_cor_gmp))




sd1=apply(data_atac_sig_sf,2,sd)


data_atac_sig_sf_sd = t(apply(data_atac_sig_sf, 1, FUN=function(x) x/norm_info[c(1:12),18]*norm_info[c(1:12),14])) #data_atac_sig_sf#
data_atac_sig_sf_sd_am = data_atac_sig_sf_am#t(apply(data_atac_sig_sf_am, 1, FUN=function(x) x/norm_info[c(1:12),15]*norm_info[c(1:12),11])) #data_atac_sig_sf_am#
colnames(data_atac_sig_sf_sd) = colnames(data_atac_sig)
colnames(data_atac_sig_sf_sd_am) = colnames(data_atac_sig)


###### input matrix
### od log2
data_atac_sig_log2 = log2(data_atac_sig+1)
### scale log2
data_atac_sig_sf_sd_am_log2 = log2(data_atac_sig_sf_sd_am+1)
### sf_sd_norm log2
data_atac_sig_sf_sd_log2 = log2(data_atac_sig_sf_sd+1)
### rna matrix
data_rna_sig = data_rna[,c(-1,-2,-3,-4,-5)]



###### scale cor
data_atac_sf_sd_rna_sig_cor_cmp = cor(data_atac_sig_scale_log2_cmp, data_rna_sig_cmp, method='spearman')
summary((data_atac_sf_sd_rna_sig_cor_cmp))
data_atac_sf_sd_rna_sig_cor_cmp = cor(data_atac_sig_scale_log2_cmp, data_rna_sig_cmp, method='pearson')
summary((data_atac_sf_sd_rna_sig_cor_cmp))

data_atac_sf_sd_rna_sig_cor_cmp = cor(data_atac_sig_scale_log2_ery, data_rna_sig_ery, method='spearman')
summary((data_atac_sf_sd_rna_sig_cor_cmp))
data_atac_sf_sd_rna_sig_cor_cmp = cor(data_atac_sig_scale_log2_ery, data_rna_sig_ery, method='pearson')
summary((data_atac_sf_sd_rna_sig_cor_cmp))


data_atac_sf_sd_rna_sig_cor_cmp = cor(data_atac_sig_scale_log2_gmp, data_rna_sig_gmp, method='spearman')
summary((data_atac_sf_sd_rna_sig_cor_cmp))
data_atac_sf_sd_rna_sig_cor_cmp = cor(data_atac_sig_scale_log2_gmp, data_rna_sig_gmp, method='pearson')
summary((data_atac_sf_sd_rna_sig_cor_cmp))



data_atac_sf_sd_rna_sig_cor_cmp = apply(cbind(data_atac_sig_sf_sd_log2_cmp, data_rna_sig_cmp), MARGIN=1, FUN=function(x) cor(x[1], x[2], method='spearman') )
summary((data_atac_sf_sd_rna_sig_cor_cmp))


par(mfrow=c(2,2))
heatscatter(data_atac_sig_log2_ery, data_rna_sig_ery, pch = 20)
heatscatter(data_atac_sig_scale_log2_ery, data_rna_sig_ery, pch = 20)
heatscatter(data_atac_sig_sf_sd_log2_ery, data_rna_sig_ery, pch = 20)

heatscatter(data_atac_sig_log2_cmp, data_atac_sig_log2_gmp, pch = 20)
heatscatter(data_atac_sig_sf_sd_log2_cmp, data_atac_sig_sf_sd_log2_gmp, pch = 20)

heatscatter(data_rna_sig_cmp, data_rna_sig_gmp, pch = 20)

### plot variance
par(mfrow=c(2,2))
heatscatter(log(data_atac_sig_var), log(data_rna_sig_var), pch = 20, xlim=c(-8, 2))
heatscatter(log(data_atac_sig_scale_log2_var), log(data_rna_sig_var), pch = 20, xlim=c(-8, 2))
heatscatter(log(data_atac_sig_sf_sd_log2_var), log(data_rna_sig_var), pch = 20, xlim=c(-8, 2), ylim=c(-8, 2))
heatscatter(log(data_atac_sig_sf_sd_log2_var), log(data_atac_sig_var), pch = 20, xlim=c(-8, 2), ylim=c(-8, 2))

heatscatter(log(data_atac_sig_sf_sd_log2_var), log(data_rna_sig_var), pch = 20, xlim=c(-2,2))
heatscatter(log(data_atac_sig_var), log(data_rna_sig_var), pch = 20, xlim=c(-2,2))


par(mfrow=c(2,2))
heatscatter(log(data_atac_sig_var)[log(data_rna_sig_var)<(-4)], log(data_rna_sig_var)[log(data_rna_sig_var)<(-4)], pch = 20, xlim=c(-6.5, 1.5))
heatscatter(log(data_atac_sig_scale_log2_var)[log(data_rna_sig_var)<(-4)], log(data_rna_sig_var)[log(data_rna_sig_var)<(-4)], pch = 20, xlim=c(-6.5, 1.5))
heatscatter(log(data_atac_sig_sf_sd_log2_var)[log(data_rna_sig_var)<(-4)], log(data_rna_sig_var)[log(data_rna_sig_var)<(-4)], pch = 20, xlim=c(-6.5, 1.5))
heatscatter(log(data_atac_sig_sf_sd_log2_var), log(data_atac_sig_var), pch = 20)

par(mfrow=c(2,2))
heatscatter(log(data_atac_sig_var)[log(data_rna_sig_var)>=(-4)], log(data_rna_sig_var)[log(data_rna_sig_var)>=(-4)], pch = 20, xlim=c(-8, 2))
heatscatter(log(data_atac_sig_scale_log2_var)[log(data_rna_sig_var)>=(-4)], log(data_rna_sig_var)[log(data_rna_sig_var)>=(-4)], pch = 20, xlim=c(-8, 2))
heatscatter(log(data_atac_sig_sf_sd_log2_var)[log(data_rna_sig_var)>=(-4)], log(data_rna_sig_var)[log(data_rna_sig_var)>=(-4)], pch = 20, xlim=c(-8, 2))
heatscatter(log(data_atac_sig_sf_sd_log2_var), log(data_atac_sig_var), pch = 20)

### plot sd
heatscatter(log(data_atac_sig_sd), log(data_rna_sig_sd), pch = 20)
heatscatter(log(data_atac_sig_sf_sd_log2_sd), log(data_rna_sig_sd), pch = 20)
heatscatter(log(data_atac_sig_sf_sd_log2_sd), log(data_rna_sig_sd), pch = 20)
heatscatter(log(data_atac_sig_sf_sd_log2_sd), log(data_atac_sig_sd), pch = 20)

par(mfrow=c(1,1))
hist(log(data_rna_sig_var), breaks=50)

plot(density(log(data_rna_sig_var)))

plot(density((data_rna_sig_var)))

hist((data_rna_sig_var), breaks=50)


data_atac_cmp = data_atac[,8]
data_atac_gmp = data_atac[,9]

data_rna_cmp = data_rna[,7]
data_rna_gmp = data_rna[,8]

### covert to rank
data_atac_sig_rank = t(apply(data_atac_sig, 1, rank))
data_atac_sig_sf_sd_log2_rank = t(apply(data_atac_sig_sf_sd_log2, 1, rank))
data_rna_sig_rank = t(apply(data_rna_sig, 1, rank))

### spearman correlation
data_atac_rna_sig_cor = apply(cbind(data_atac_sig, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method='spearman') )
data_atac_scale_rna_sig_cor = apply(cbind(data_atac_sig_scale_log2, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method='spearman') )
data_atac_sf_sd_rna_sig_cor = apply(cbind(data_atac_sig_sf_sd_log2, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method='spearman') )

summary((data_atac_rna_sig_cor))
summary((data_atac_scale_rna_sig_cor))
summary((data_atac_sf_sd_rna_sig_cor))

par(mfrow=c(2,2))
hist(abs(data_atac_rna_sig_cor), breaks = 50)
hist(abs(data_atac_sf_sd_rna_sig_cor), breaks = 50)


par(mfrow=c(2,2))
plot(density(abs(data_atac_rna_sig_cor[!is.na(data_atac_rna_sig_cor)])), ylim=c(0,2.5))
plot(density(abs(data_atac_scale_rna_sig_cor[!is.na(data_atac_scale_rna_sig_cor)])), ylim=c(0,2.5))
plot(density(abs(data_atac_sf_sd_rna_sig_cor[!is.na(data_atac_sf_sd_rna_sig_cor)])), ylim=c(0,2.5))

png('spearman_corr.png')
par(mfrow=c(2,2))
plot(density((data_atac_rna_sig_cor[!is.na(data_atac_rna_sig_cor)])), ylim=c(0,2.5))
plot(density((data_atac_scale_rna_sig_cor[!is.na(data_atac_scale_rna_sig_cor)])), ylim=c(0,2.5))
plot(density((data_atac_sf_sd_rna_sig_cor[!is.na(data_atac_sf_sd_rna_sig_cor)])), ylim=c(0,2.5))
dev.off()

par(mfrow=c(2,2))
plot(density((data_atac_rna_sig_cor[log(data_rna_sig_var)<(-4)][!is.na(data_atac_rna_sig_cor[log(data_rna_sig_var)<(-4)])])), ylim=c(0,4), xlim=c(-1,1))
plot(density((data_atac_scale_rna_sig_cor[log(data_rna_sig_var)<(-4)][!is.na(data_atac_scale_rna_sig_cor[log(data_rna_sig_var)<(-4)])])), ylim=c(0,4), xlim=c(-1,1))
plot(density((data_atac_sf_sd_rna_sig_cor[log(data_rna_sig_var)<(-4)][!is.na(data_atac_sf_sd_rna_sig_cor[log(data_rna_sig_var)<(-4)])])), ylim=c(0,4), xlim=c(-1,1))


par(mfrow=c(2,2))
plot(density((data_atac_rna_sig_cor[log(data_rna_sig_var)>=(-4)][!is.na(data_atac_rna_sig_cor[log(data_rna_sig_var)>=(-4)])])), ylim=c(0,1.5))
plot(density((data_atac_scale_rna_sig_cor[log(data_rna_sig_var)>=(-4)][!is.na(data_atac_scale_rna_sig_cor[log(data_rna_sig_var)>=(-4)])])), ylim=c(0,1.5))
plot(density((data_atac_sf_sd_rna_sig_cor[log(data_rna_sig_var)>=(-4)][!is.na(data_atac_sf_sd_rna_sig_cor[log(data_rna_sig_var)>=(-4)])])), ylim=c(0,1.5))


### pearson correlation
data_atac_rna_sig_cor_pearson = apply(cbind(data_atac_sig, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method='pearson') )
data_atac_sf_sd_rna_sig_cor_pearson = apply(cbind(data_atac_sig_sf_sd, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method='pearson') )

summary(abs(data_atac_rna_sig_cor_pearson))
summary(abs(data_atac_sf_sd_rna_sig_cor_pearson))

hist(abs(data_atac_rna_sig_cor_pearson), breaks = 50)
hist(abs(data_atac_sf_sd_rna_sig_cor_pearson), breaks = 50)




plot(data_atac_sig_var, data_rna_sig_var)
heatscatter(data_atac_sig_var, data_rna_sig_var, pch = 20)#, ylim=c(-10,5), xlim=c(20,10000))
heatscatter(log(data_atac_sig_var), log(data_rna_sig_var), pch = 20)#, ylim=c(-10,5), xlim=c(20,10000))
plot(data_atac_sig_rank_var, data_rna_sig_var)
heatscatter(data_atac_cmp, data_atac_cmp, pch = 20)
'''
