### get parameters
args = commandArgs(trailingOnly=TRUE)
index_matrix_ideas_state_inputfile = args[1]
ideas_state_color = args[2]
cREs_IDEASpro_outfile = args[3]

### read index set matrix
read_matrix = function(inputfile){
	data_index_set = as.matrix(read.table(inputfile, header=T))
	rownames(data_index_set) = data_index_set[,1]
	data_index_set = data_index_set[,-1]
	class(data_index_set) = "numeric" 
	return(data_index_set)
}

###########
ideas_state_matrix = read_matrix(index_matrix_ideas_state_inputfile)

### get unique elements from the matrix
ideas_state_matrix_flatten = as.vector(ideas_state_matrix)
ideas_state_matrix_uniq = unique(ideas_state_matrix_flatten)
### sort index
ideas_state_matrix_uniq = sort(ideas_state_matrix_uniq)

###### extract counts matrix
counts_matrix = c()
for (i in c(1: dim(ideas_state_matrix)[2]) ){
	### extract ith cell type data
	ideas_state_matrix_table_tmp = as.matrix(ideas_state_matrix[,i])

	table_tmp = c()
	for (j in c( 1: length(ideas_state_matrix_uniq)) ){
		### count the number of cREs have jth IDEAS state
		table_tmp[j] = sum(ideas_state_matrix_table_tmp==(j-1))
	}

	### vector to matrix
	counts_matrix = rbind(counts_matrix, table_tmp)
}

### transpose matrix
counts_matrix_t = t( counts_matrix )
### add colnames
colnames(counts_matrix_t) = colnames(ideas_state_matrix)

### set heatmap colors
rgb_col_num=read.table(ideas_state_color,header=F)
rgb_col_num=rbind(rgb_col_num, c(255,255,255))
print(rgb_col_num)
rgb_col=apply(rgb_col_num,1,function(x) rgb(x[1],x[2],x[3],max=255))

my_colorbar=colorRampPalette(rgb_col)(n = 18)
col_breaks = c(
        seq(0.1, 1,length=2),
        seq(1.1, 2,length=2),
        seq(2.1, 3,length=2),
        seq(3.1, 4,length=2),
        seq(4.1, 5,length=2),
        seq(5.1, 6,length=2),
        seq(6.1, 7,length=2),
        seq(7.1, 8,length=2),
        seq(8.1, 9,length=2),
        seq(9.1, 10,length=2),
        seq(10.1, 11,length=2),
        seq(11.1, 12,length=2),
        seq(12.1, 13,length=2),
        seq(13.1, 14,length=2),
        seq(14.1, 15,length=2),
        seq(15.1, 16,length=2),
        seq(16.1, 17,length=2),
        seq(17.1, 18,length=2)
)

### save figure
pdf(cREs_IDEASpro_outfile, 16, 16)
barplot(counts_matrix_t, col=my_colorbar)
dev.off()


# Rscript plot_index_set_region_hist.R celltype.index.sorted.txt cRE_celltype_number.pdf

