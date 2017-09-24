library(pheatmap)
### get parameters
args = commandArgs(trailingOnly=TRUE)
index_set_inputfile = args[1]
index_inputfile = args[2]
index_set_inputfile_binary = args[3]

index_set_all_heatmap = args[4]
index_set_thresh_heatmap = args[5]
index_all_heatmap = args[6]

color = args[7]
threshold = as.numeric(args[8])
small_num_for_log2 = as.numeric(args[9])
tranform = args[10]

### set heatmap colors
my_colorbar=colorRampPalette(c('white',color))(n = 128)
col_breaks = c(seq(0, 2000,length=33))

### read index set signal matrix
read_index_set_sig_matrix = function(inputfile){
	data_index_set = as.matrix(read.table(inputfile, header=T))
	rownames(data_index_set) = data_index_set[,1]
	data_index_set = data_index_set[,-1]
	class(data_index_set) = "numeric" 
	return(data_index_set)
}

### read index set matrix
read_matrix = function(inputfile){
	data_index_set = as.matrix(read.table(inputfile, header=T))
	rownames(data_index_set) = data_index_set[,1]
	data_index_set = data_index_set[,-1]
	class(data_index_set) = "numeric" 
	return(data_index_set)
}

###########
data_index_set_sig = read_index_set_sig_matrix(index_set_inputfile)
if (tranform == 'log2'){
	data_index_set_sig = log2(data_index_set_sig + small_num_for_log2)
}
data_index_set_01 = read_matrix(index_set_inputfile_binary)

### log2 scale
if (tranform == 'log2'){
	data_index_set_01 = log2(data_index_set_01+small_num_for_log2)
}
### plot index set heatmap
pheatmap(data_index_set_sig, color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE, filename = index_set_all_heatmap)
### plot index set heatmap (filter by DNA region number threshold)
data_index_set_filtered = data_index_set_sig[apply(data_index_set_01, 1, max) >= log2(threshold), ]
pheatmap(data_index_set_filtered, color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE, filename = index_set_thresh_heatmap)
###########
### read index matrix
data_index = read_matrix(index_inputfile)
if (tranform == 'log2'){
	data_index = log2(data_index + small_num_for_log2)
}

### plot index heatmap
pheatmap(data_index, color=my_colorbar, cluster_cols = FALSE, cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE, filename = index_all_heatmap)
###########

# Rscript plot_index_set_module.R celltype.index_set.sorted.txt celltype.binary_pattern.sorted.txt black 200 index_set_all.pdf index_set_thresh.pdf index.png