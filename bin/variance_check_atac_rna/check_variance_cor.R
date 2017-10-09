library(LSD)
data_atac = read.table('gene_atac/gencode_pc_sort.atac.celltype.matched.txt', header = T, sep='\t')

norm_info = read.table('scale_factor_matrix/ncis_table_list.sf.txt')

scale_factor_allmean = norm_info[c(1:12),7]

scale_factor = norm_info[c(1:12),10]


data_rna = read.table('rsem_matrix_folder/rsem_matrix.norm.rld_matrix.celltype.txt', header = T, sep='\t')
#data_rna = read.table('rsem_matrix_folder/rsem_matrix.norm.log2_norm_matrix_plus1.celltype.txt', header = T, sep='\t')

data_atac_sig = data_atac[,c(-1,-2,-3,-4,-5,-6)]

### scale factor norm
data_atac_sig_sf = t(apply(data_atac_sig, 1, FUN=function(x) x*scale_factor))
data_atac_sig_sf_am = t(apply(data_atac_sig, 1, FUN=function(x) x*scale_factor_allmean))

### sd norm
#col_sd = apply(data_atac_sig_sf, 2, sd)
#data_atac_sig_sf_sd = t(apply(data_atac_sig_sf, 1, FUN=function(x) x/col_sd))
data_atac_sig_sf_sd = data_atac_sig_sf_sd#*5000
data_atac_sig_sf_sd_am = data_atac_sig_sf_am#*5000
colnames(data_atac_sig_sf_sd) = colnames(data_atac_sig)
colnames(data_atac_sig_sf_sd_am) = colnames(data_atac_sig)


###### input matrix
### od log2
data_atac_sig_log2 = log2(data_atac_sig+1)
### scale log2
data_atac_sig_scale_log2 = scale(log2(data_atac_sig_sf_sd_am+1))
### sf_sd_norm log2
data_atac_sig_sf_sd_log2 = log2(data_atac_sig_sf_sd+1)
### rna matrix
data_rna_sig = data_rna[,c(-1,-2,-3,-4,-5)]





data_atac_sig_log2_cmp = data_atac_sig_log2[,2]
data_atac_sig_log2_gmp = data_atac_sig_log2[,3]
data_atac_sig_log2_ery = data_atac_sig_log2[,6]

data_atac_sig_scale_log2_cmp = data_atac_sig_scale_log2[,2]
data_atac_sig_scale_log2_gmp = data_atac_sig_scale_log2[,3]
data_atac_sig_scale_log2_ery = data_atac_sig_scale_log2[,6]

data_atac_sig_sf_sd_log2_cmp = data_atac_sig_sf_sd_log2[,2]
data_atac_sig_sf_sd_log2_gmp = data_atac_sig_sf_sd_log2[,3]
data_atac_sig_sf_sd_log2_ery = data_atac_sig_sf_sd_log2[,6]

data_rna_sig_cmp = data_rna_sig[,2]
data_rna_sig_gmp = data_rna_sig[,3]
data_rna_sig_ery = data_rna_sig[,6]


data_atac_sig_var = apply(data_atac_sig_log2, 1, var)
data_atac_sig_scale_log2_var = apply(data_atac_sig_scale_log2, 1, var)

data_atac_sig_sf_sd_log2_var = apply(data_atac_sig_sf_sd_log2, 1, var)
data_rna_sig_var = apply(data_rna_sig, 1, var)

data_atac_sig_sd = apply(data_atac_sig_log2, 1, sd)
data_atac_sig_sf_sd_log2_sd = apply(data_atac_sig_sf_sd_log2, 1, sd)
data_rna_sig_sd = apply(data_rna_sig, 1, sd)


par(mfrow=c(2,2))
heatscatter(data_atac_sig_log2_cmp, data_rna_sig_cmp, pch = 20)
abline(0,1)
heatscatter(data_atac_sig_scale_log2_cmp, data_rna_sig_cmp, pch = 20)
abline(0,1)
heatscatter(data_atac_sig_sf_sd_log2_cmp, data_rna_sig_cmp, pch = 20)
abline(0,1)

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
data_atac_sf_sd_rna_sig_cor = apply(cbind(data_atac_sig_sf_sd, data_rna_sig), MARGIN=1, FUN=function(x) cor(x[1:12], x[13:24], method='spearman') )

summary((data_atac_rna_sig_cor))
summary((data_atac_scale_rna_sig_cor))
summary((data_atac_sf_sd_rna_sig_cor))

hist(abs(data_atac_rna_sig_cor), breaks = 50)
hist(abs(data_atac_sf_sd_rna_sig_cor), breaks = 50)

par(mfrow=c(2,2))
plot(density(abs(data_atac_rna_sig_cor[!is.na(data_atac_rna_sig_cor)])), ylim=c(0,2.5))
plot(density(abs(data_atac_scale_rna_sig_cor[!is.na(data_atac_scale_rna_sig_cor)])), ylim=c(0,2.5))
plot(density(abs(data_atac_sf_sd_rna_sig_cor[!is.na(data_atac_sf_sd_rna_sig_cor)])), ylim=c(0,2.5))

par(mfrow=c(2,2))
plot(density((data_atac_rna_sig_cor[!is.na(data_atac_rna_sig_cor)])), ylim=c(0,2.5))
plot(density((data_atac_scale_rna_sig_cor[!is.na(data_atac_scale_rna_sig_cor)])), ylim=c(0,2.5))
plot(density((data_atac_sf_sd_rna_sig_cor[!is.na(data_atac_sf_sd_rna_sig_cor)])), ylim=c(0,2.5))


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

