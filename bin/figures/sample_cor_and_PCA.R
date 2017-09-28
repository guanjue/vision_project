###### run correlation and cell type distance & PCA
### get parameters
args = commandArgs(trailingOnly = TRUE)
input_file = args[1]
output_file = args[2]
start_col = as.numeric(args[3])
celltype_distfun = args[4]

### read index set signal matrix
read_sig_matrix = function(inputfile){
	data_index_set = as.matrix(read.table(inputfile, header=T))
	rownames(data_index_set) = data_index_set[,4]
	print(dim(data_index_set))
	data_index_set = data_index_set[,c(start_col:(start_col+27))]
	class(data_index_set) = "numeric" 
	return(data_index_set)
}

###########
###### read input file
### read signal matrix
sig_matrix = read_sig_matrix(input_file)

### transpose matrix
sig_matrix_t = t(sig_matrix)

### run PCA
print('run PCA:')
pca=prcomp(sig_matrix_t, scale=TRUE)
print('run PCA: DONE')

###### plot PCA
library(ggplot2)
scores=data.frame(pca$x)
cols=ncol(scores)

### get point labels
cells=gsub("_.*", "", rownames(scores))
#print(cells)

### get coordinates
scores=data.frame(scores,cells)
colnames(scores)[cols+1]="cell"

### plot PCA plot
pdf(paste(output_file, 'PCA.pdf', sep='_'))
par(mfrow=c(2,2))
ggplot(scores, aes(x=PC1, y=PC2))+geom_point(aes(color=cell), size=6)+geom_text(aes(label=rownames(scores), hjust=-.2), size=3)
ggplot(scores, aes(x=PC1, y=PC3))+geom_point(aes(color=cell), size=6)+geom_text(aes(label=rownames(scores), hjust=-.2), size=3)
ggplot(scores, aes(x=PC1, y=PC4))+geom_point(aes(color=cell), size=6)+geom_text(aes(label=rownames(scores), hjust=-.2), size=3)
ggplot(scores, aes(x=PC2, y=PC3))+geom_point(aes(color=cell), size=6)+geom_text(aes(label=rownames(scores), hjust=-.2), size=3)
ggplot(scores, aes(x=PC2, y=PC4))+geom_point(aes(color=cell), size=6)+geom_text(aes(label=rownames(scores), hjust=-.2), size=3)
ggplot(scores, aes(x=PC3, y=PC4))+geom_point(aes(color=cell), size=6)+geom_text(aes(label=rownames(scores), hjust=-.2), size=3)
dev.off()

#variance plot
#http://strata.uga.edu/6370/lecturenotes/principalComponents.html
sd=pca$sdev
var=sd^2
var.percent=var/sum(var) * 100
pdf(paste(output_file, 'Percent_Variance.pdf', sep='_'))
barplot(var.percent, xlab="PC", ylab="Percent Variance", names.arg=1:length(var.percent), las=1, ylim=c(0,max(var.percent)), col="blue")
abline(h=1/ncol(pca)*100, col="red")
dev.off()



### calculate cell type spearman correlation
c=cor(sig_matrix, method=celltype_distfun) # use spearman to match ENCODE
print(head(sig_matrix))
### heatmap labels
cr=round(c, digits=2)

### get color vector
my_palette <- c("#1000FFFF","#0047FFFF","#00C8FFFF","#00FFE0FF","#00FF89FF","#00FF33FF","#99FF00FF","#AFFF00FF","#FFD800FF","#FFB800FF","#FFAD00FF","#FF8100FF","#FF5600FF","#FF2B00FF","#FF0000FF")
breaks = c(seq(-0.1, 1.0, length=16))

### plot spearman correlation
library(gplots)
pdf(paste(output_file, 'spearman_correlation.pdf', sep='_'))
heatmap.2(c, dendrogram="row", col=my_palette, breaks=breaks, trace="none", cellnote=cr, notecol="black", revC=TRUE, notecex=.5, symm=TRUE, distfun=dist, symkey=FALSE)
dev.off()






