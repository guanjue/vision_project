##################################
script_folder='/Volumes/MAC_Data/data/labs/hardison_lab/vision/bin/'
analysis_folder='/Volumes/MAC_Data/data/labs/hardison_lab/vision/rnaseq/'
##################################
	###### initiate folders
	gene_list='gene_list/'
	rsem_folder='rsem/'
	rsem_sort_folder='rsem_sort/'

	if [ -d "$rsem_sort" ]; then  
		rm -r $rsem_sort
	fi
	mkdir $rsem_sort

### go to folder
cd $analysis_folder

###### convert gencode.vM4.annotation.gtf to bed 
### extract wanted info (extract only genes)
cat $gene_list'gencode.vM4.annotation.gtf' | awk -F '\t' -v OFS='\t' '{if ($3=="gene") print $1,$4-1,$5,$3,$7,$9}' | awk -F ';' -v OFS='\t' '{print $1, $3, $5}' > $gene_list'gencode.txt'

### convert to bed format 
cat $gene_list'gencode.txt' | awk -F '"' -v OFS='\t' '{print $1,$2,$4,$6}' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3, $7";"$8, $9, $5}' > $gene_list'gencode.bed'
### convert to bed format (protein_coding only)
cat $gene_list'gencode.txt' | awk -F '"' -v OFS='\t' '{print $1,$2,$4,$6}' | awk -F '\t' -v OFS='\t' '{if ($8=="protein_coding") print $1,$2,$3, $7";"$8, $9, $5}' > $gene_list'gencode_protein_coding.bed'

### sort the protein_coding bed
sort -k1,1 -k2,2n $gene_list'gencode_protein_coding.bed' > $gene_list'gencode_pc_sort.bed'

### get expand gene body
cat $gene_list'gencode_pc_sort.bed' | awk -F '\t' -v OFS='\t' '{if ($2-10000 > 0) print $1,$2-10000,$3+10000,$4,$5,$6}' > $gene_list'gencode_pc_sort.exp10kb.bed'
### get TSS
cat $gene_list'gencode_pc_sort.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+" && $2-10000>0) print $1,$2-10000,$2+10000,$4,$5,$6; else if ($6=="-" && $3-10000>0) print $1,$3-10000,$3+10000,$4,$5,$6}' > $gene_list'gencode_pc_sort.TSSexp10kb.bed'

### scp $gene_list'gencode_pc_sort'*'.bed' gzx103@biostar.psu.edu:/gpfs/home/gzx103/scratch/vision_clustering/gene_atac_sig/


###### row sort RSEM matrix for each sample
for rsem_bed in $(cat rsem_list.txt)
do
	echo $rsem_bed
	python $script_folder'rna_matrix/vlookup_bed.py' -t $rsem_folder$rsem_bed -s $gene_list'gencode_pc_sort.bed' -o $rsem_sort_folder$rsem_bed'.sort.bed' -k T
done

### merge to rsem matrix
cat $rsem_sort_folder'LSK.rsem.1007.bed.sort.bed' | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $16";"$17, $6}' > $rsem_sort_folder'rsem_matrix.txt'
for rsem_bed in $(cat rsem_list.txt)
do
	echo $rsem_bed
	### extract signal column
	cat $rsem_sort_folder$rsem_bed'.sort.bed' | awk -F '\t' -v OFS='\t' '{print $4}' > tmp_sig.txt
	### merge matrix
	paste $rsem_sort_folder'rsem_matrix.txt' tmp_sig.txt > $rsem_sort_folder'rsem_matrix.tmp.txt'
	### change names
	mv $rsem_sort_folder'rsem_matrix.tmp.txt' $rsem_sort_folder'rsem_matrix.txt'
	### rm tmp
	rm tmp_sig.txt
done


