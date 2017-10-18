##################################
script_folder='/Volumes/MAC_Data/data/labs/hardison_lab/vision/bin/'
analysis_folder='/Volumes/MAC_Data/data/labs/hardison_lab/vision/gene_rnaseq_atac/'
##################################
	###### initiate folders
	gene_list='gene_list/'
	rsem_folder='rsem/'
	rsem_sort_folder='rsem_sort/'
	rsem_matrix_folder='rsem_matrix_folder/'
	input_folder='input_folder/'
	gene_atac='gene_atac/'
	atac_ncis_table='atac_ncis_table/'
	atac_ncis_table_plot='atac_ncis_table_plot/'
	scale_factor_matrix='scale_factor_matrix/'
	atac_reads_table='atac_reads_table/'


	if [ -d "$rsem_sort" ]; then  
		rm -r $rsem_sort
	fi
	mkdir $rsem_sort

	if [ -d "$rsem_matrix_folder" ]; then  
		rm -r $rsem_matrix_folder
	fi
	mkdir $rsem_matrix_folder

	if [ -d "$atac_ncis_table_plot" ]; then  
		rm -r $atac_ncis_table_plot
	fi
	mkdir $atac_ncis_table_plot

	if [ -d "$scale_factor_matrix" ]; then  
		rm -r $scale_factor_matrix
	fi
	mkdir $scale_factor_matrix

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
	cat $gene_list'gencode_pc_sort.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+" && $2-5000>0) print $1,$2-10000,$2+5000,$4,$5,$6; else if ($6=="-" && $3-5000>0) print $1,$3-5000,$3+5000,$4,$5,$6}' > $gene_list'gencode_pc_sort.TSSexp5kb.bed'
	cat $gene_list'gencode_pc_sort.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+" && $2-1000>0) print $1,$2-10000,$2+1000,$4,$5,$6; else if ($6=="-" && $3-1000>0) print $1,$3-1000,$3+1000,$4,$5,$6}' > $gene_list'gencode_pc_sort.TSSexp1kb.bed'

	cat $gene_list'gencode_pc_sort.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+" && $2-100000>0) print $1,$2-100000,$2+100000,$4,$5,$6; else if ($6=="-" && $3-100000>0) print $1,$3-100000,$3+100000,$4,$5,$6}' > $gene_list'gencode_pc_sort.TSSexp100kb.bed'
	cat $gene_list'gencode_pc_sort.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+" && $2-110000>0) print $1,$2-110000,$2+110000,$4,$5,$6; else if ($6=="-" && $3-110000>0) print $1,$3-110000,$3+110000,$4,$5,$6}' > $gene_list'gencode_pc_sort.TSSexp110kb.bed'

	cat $gene_list'gencode_pc_sort.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+" && $2-100000-5000>0) print $1,$2-100000-5000,$2-100000,$4,$5,$6; else if ($6=="-" && $3-100000-5000>0) print $1,$3+100000,$3+100000+5000,$4,$5,$6}' > $gene_list'gencode_pc_sort.TSSup100kb.bed'
	cat $gene_list'gencode_pc_sort.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+" && $2-100000-5000>0) print $1,$2+100000,$2+100000+5000,$4,$5,$6; else if ($6=="-" && $3-100000-5000>0) print $1,$3-100000-5000,$3-100000,$4,$5,$6}' > $gene_list'gencode_pc_sort.TSSdown100kb.bed'

	cat $gene_list'gencode_pc_sort.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+" && $2-1000000-5000>0) print $1,$2-1000000-5000,$2-1000000,$4,$5,$6; else if ($6=="-" && $3-1000000-5000>0) print $1,$3+1000000,$3+1000000+5000,$4,$5,$6}' > $gene_list'gencode_pc_sort.TSSup1000kb.bed'
	cat $gene_list'gencode_pc_sort.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+" && $2-1000000-5000>0) print $1,$2+1000000,$2+1000000+5000,$4,$5,$6; else if ($6=="-" && $3-1000000-5000>0) print $1,$3-1000000-5000,$3-1000000,$4,$5,$6}' > $gene_list'gencode_pc_sort.TSSdown1000kb.bed'

	bedtools intersect -a $gene_list'gencode_pc_sort.TSSexp110kb.bed' -b $gene_list'gencode_pc_sort.TSSexp100kb.bed' > $gene_list'gencode_pc_sort.TSSexp110kb_sub_100kb.bed'


### scp $gene_list'gencode_pc_sort'*'.bed' gzx103@biostar.psu.edu:/gpfs/home/gzx103/scratch/vision_clustering/gene_atac_sig/


###### get rsem cell type matrix for each gene
	### row sort RSEM matrix for each sample (the gtf gene coordinates and rsem signal matrix gene coordinates are not the same)
	while read rsem_bed
	do
		echo $rsem_bed
		python $script_folder'rna_matrix/vlookup_bed.py' -t $rsem_folder$rsem_bed -s $gene_list'gencode_pc_sort.bed' -o $rsem_sort_folder$rsem_bed'.sort.bed' -k T
	done < $input_folder'rsem_list.txt'

	### merge to rsem sample matrix
	cat $rsem_sort_folder'LSK.rsem.1007.bed.sort.bed' | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $16";"$17, $6}' > $rsem_matrix_folder'rsem_matrix.txt'
	while read LINE
	do
		rsem_bed=$(echo "$LINE" | awk '{print $1}')
		cell_order=$(echo "$LINE" | awk '{print $2}')
		echo $rsem_bed
		echo $cell_order
		### extract signal column
		cat $rsem_sort_folder$rsem_bed'.sort.bed' | awk -F '\t' -v OFS='\t' '{print $4}' > tmp_sig.txt
		### merge matrix
		paste $rsem_matrix_folder'rsem_matrix.txt' tmp_sig.txt > $rsem_matrix_folder'rsem_matrix.tmp.txt'
		### change names
		mv $rsem_matrix_folder'rsem_matrix.tmp.txt' $rsem_matrix_folder'rsem_matrix.txt'
		### rm tmp
		rm tmp_sig.txt
	done < $input_folder'rsem_list.txt'

	### merge rsem sample matrix to cell type matrix
	time python $script_folder'rna_matrix/merge_cell_type_data_rsem.py' -i $rsem_matrix_folder'rsem_matrix.txt' -m $input_folder'rsem_list_sample2celltype.txt' -n 6 -o $rsem_matrix_folder'rsem_matrix.celltype.txt'

	### sort columns without merge cell types
	time python $script_folder'rna_matrix/merge_cell_type_data_rsem.py' -i $rsem_matrix_folder'rsem_matrix.txt' -m $input_folder'rsem_list_sample2sortsample.txt' -n 6 -o $rsem_matrix_folder'rsem_matrix.sortsample.txt'

###### normalize the rsem relative to CMP
	### quantile normalization
	time python $script_folder'rna_matrix/quantile_normalization.py' -i $rsem_matrix_folder'rsem_matrix.celltype.txt' -n 6 -o $rsem_matrix_folder'rsem_matrix.celltype.qn.txt'

	### DEseq2 normalization
	time Rscript $script_folder'rna_matrix/DEseq2_norm.R' $rsem_matrix_folder'rsem_matrix.txt' $input_folder'rsem_list_sample_col.txt' $rsem_matrix_folder'rsem_matrix.norm'

	### merge normed rsem sample matrix to cell type matrix
	time python $script_folder'rna_matrix/merge_cell_type_data_rsem.py' -i $rsem_matrix_folder'rsem_matrix.norm.rld_matrix.txt' -m $input_folder'rsem_list_sample2celltype.txt' -n 6 -o $rsem_matrix_folder'rsem_matrix.norm.rld_matrix.celltype.txt'
	time python $script_folder'rna_matrix/merge_cell_type_data_rsem.py' -i $rsem_matrix_folder'rsem_matrix.norm.log2_norm_matrix_plus1.txt' -m $input_folder'rsem_list_sample2celltype.txt' -n 6 -o $rsem_matrix_folder'rsem_matrix.norm.log2_norm_matrix_plus1.celltype.txt'


###### gene atac
	### merge normed rsem sample matrix to cell type matrix
	time python $script_folder'rna_matrix/merge_cell_type_data_rsem.py' -i $gene_atac'gencode_pc_sort.atac.txt' -m $input_folder'gene_atac_list_sample2celltype.txt' -n 7 -o $gene_atac'gencode_pc_sort.atac.celltype.txt'
	time python $script_folder'rna_matrix/merge_cell_type_data_rsem.py' -i $gene_atac'gencode_pc_sort.TSSexp1kb.atac.txt' -m $input_folder'gene_atac_list_sample2celltype.txt' -n 7 -o $gene_atac'gencode_pc_sort.TSSexp1kb.atac.celltype.txt'
	time python $script_folder'rna_matrix/merge_cell_type_data_rsem.py' -i $gene_atac'gencode_pc_sort.TSSexp5kb.atac.txt' -m $input_folder'gene_atac_list_sample2celltype.txt' -n 7 -o $gene_atac'gencode_pc_sort.TSSexp5kb.atac.celltype.txt'
	time python $script_folder'rna_matrix/merge_cell_type_data_rsem.py' -i $gene_atac'gencode_pc_sort.TSSexp10kb.atac.txt' -m $input_folder'gene_atac_list_sample2celltype.txt' -n 7 -o $gene_atac'gencode_pc_sort.TSSexp10kb.atac.celltype.txt'
	time python $script_folder'rna_matrix/merge_cell_type_data_rsem.py' -i $gene_atac'gencode_pc_sort.TSSup100kb.atac.txt' -m $input_folder'gene_atac_list_sample2celltype.txt' -n 7 -o $gene_atac'gencode_pc_sort.TSSup100kb.atac.celltype.txt'
	time python $script_folder'rna_matrix/merge_cell_type_data_rsem.py' -i $gene_atac'gencode_pc_sort.TSSdown100kb.atac.txt' -m $input_folder'gene_atac_list_sample2celltype.txt' -n 7 -o $gene_atac'gencode_pc_sort.TSSdown100kb.atac.celltype.txt'
	time python $script_folder'rna_matrix/merge_cell_type_data_rsem.py' -i $gene_atac'gencode_pc_sort.TSSup1000kb.atac.txt' -m $input_folder'gene_atac_list_sample2celltype.txt' -n 7 -o $gene_atac'gencode_pc_sort.TSSup1000kb.atac.celltype.txt'
	time python $script_folder'rna_matrix/merge_cell_type_data_rsem.py' -i $gene_atac'gencode_pc_sort.TSSdown1000kb.atac.txt' -m $input_folder'gene_atac_list_sample2celltype.txt' -n 7 -o $gene_atac'gencode_pc_sort.TSSdown1000kb.atac.celltype.txt'
	time python $script_folder'rna_matrix/merge_cell_type_data_rsem.py' -i $gene_atac'gencode_pc_sort.TSSexp110kb_sub_100kb.atac.txt' -m $input_folder'gene_atac_list_sample2celltype.txt' -n 7 -o $gene_atac'gencode_pc_sort.TSSexp110kb_sub_100kb.atac.celltype.txt'

	### merge normed rsem sample matrix to cell type matrix
	time python $script_folder'rna_matrix/vlookup_bed.py' -t $gene_atac'gencode_pc_sort.atac.celltype.txt' -s $rsem_matrix_folder'rsem_matrix.norm.rld_matrix.celltype.txt' -o $gene_atac'gencode_pc_sort.atac.celltype.matched.txt' -k F

	###### get TSS expand atac-signal
	### get rsem id table
	echo 3 > $rsem_matrix_folder'rsem_matrix_id.txt'
	tail -n+2 $rsem_matrix_folder'rsem_matrix.norm.rld_matrix.celltype.txt' | awk -F ';' -v OFS=';' '{print $1,$2}' | awk -F '\t' -v OFS='\t' '{print $4}' >> $rsem_matrix_folder'rsem_matrix_id.txt'

	cat $rsem_matrix_folder'rsem_matrix.norm.rld_matrix.celltype.txt' | awk -F ';' -v OFS='\t' '{print $1,$2,$3}' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4";"$5,$6,$7,  $8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19 }' > $rsem_matrix_folder'rsem_matrix.norm.rld_matrix.celltype.notmatched.txt'

	### sort TSS atac-signal matrix
	time python $script_folder'rna_matrix/vlookup_uniq.py' -t $gene_atac'gencode_pc_sort.TSSup1000kb.atac.celltype.txt' -m 4 -s $rsem_matrix_folder'rsem_matrix_id.txt' -n 1 -o $gene_atac'gencode_pc_sort.TSSup1000kb.atac.celltype.matched.txt' -k F

	time python $script_folder'rna_matrix/vlookup_uniq.py' -t $gene_atac'gencode_pc_sort.TSSdown1000kb.atac.celltype.txt' -m 4 -s $gene_atac'gencode_pc_sort.TSSup1000kb.atac.celltype.matched.txt' -n 4 -o $gene_atac'gencode_pc_sort.TSSdown1000kb.atac.celltype.matched.txt' -k F
	time python $script_folder'rna_matrix/vlookup_uniq.py' -t $gene_atac'gencode_pc_sort.TSSup100kb.atac.celltype.txt' -m 4 -s $gene_atac'gencode_pc_sort.TSSup1000kb.atac.celltype.matched.txt' -n 4 -o $gene_atac'gencode_pc_sort.TSSup100kb.atac.celltype.matched.txt' -k F
	time python $script_folder'rna_matrix/vlookup_uniq.py' -t $gene_atac'gencode_pc_sort.TSSdown100kb.atac.celltype.txt' -m 4 -s $gene_atac'gencode_pc_sort.TSSup1000kb.atac.celltype.matched.txt' -n 4 -o $gene_atac'gencode_pc_sort.TSSdown100kb.atac.celltype.matched.txt' -k F
	time python $script_folder'rna_matrix/vlookup_uniq.py' -t $gene_atac'gencode_pc_sort.TSSexp1kb.atac.celltype.txt' -m 4 -s $gene_atac'gencode_pc_sort.TSSup1000kb.atac.celltype.matched.txt' -n 4 -o $gene_atac'gencode_pc_sort.TSSexp1kb.atac.celltype.matched.txt' -k F
	time python $script_folder'rna_matrix/vlookup_uniq.py' -t $gene_atac'gencode_pc_sort.TSSexp5kb.atac.celltype.txt' -m 4 -s $gene_atac'gencode_pc_sort.TSSup1000kb.atac.celltype.matched.txt' -n 4 -o $gene_atac'gencode_pc_sort.TSSexp5kb.atac.celltype.matched.txt' -k F
	time python $script_folder'rna_matrix/vlookup_uniq.py' -t $gene_atac'gencode_pc_sort.TSSexp10kb.atac.celltype.txt' -m 4 -s $gene_atac'gencode_pc_sort.TSSup1000kb.atac.celltype.matched.txt' -n 4 -o $gene_atac'gencode_pc_sort.TSSexp10kb.atac.celltype.matched.txt' -k F

	time python $script_folder'rna_matrix/vlookup_uniq.py' -t $rsem_matrix_folder'rsem_matrix.norm.rld_matrix.celltype.notmatched.txt' -m 4 -s $gene_atac'gencode_pc_sort.TSSup100kb.atac.celltype.matched.txt' -n 4 -o $rsem_matrix_folder'rsem_matrix.norm.rld_matrix.celltype.matched.txt' -k F
	time python $script_folder'rna_matrix/vlookup_uniq.py' -t $gene_atac'gencode_pc_sort.atac.celltype.txt' -m 4 -s $gene_atac'gencode_pc_sort.TSSup100kb.atac.celltype.matched.txt' -n 4 -o $gene_atac'gencode_pc_sort.atac.celltype.matched.txt' -k F


###### atac NCIS norm
	### get ncis table
	while read LINE
	do
		output_name=$(echo "$LINE" | awk '{print $3}')
		sig1=$(echo "$LINE" | awk '{print $4}')
		sig2=$(echo "$LINE" | awk '{print $5}')
		echo $output_name
		time python $script_folder'ncis_norm/get_ncis_t_a_b.py' -i $atac_reads_table$sig1 -j $atac_reads_table$sig2 -o $atac_ncis_table$output_name
	done < $input_folder'ncis_table_list.txt'

	### get the NCIS T-R model normed matrix
	time Rscript $script_folder'gene_atac/ncis_t_norm.R' $input_folder'ncis_table_t_norm.txt' $gene_atac'gencode_pc_sort.TSSexp10kb.atac.celltype.matched.txt' $gene_atac'gencode_pc_sort.TSSexp10kb.atac.celltype.matched.TRnormed.txt' $atac_ncis_table

	### get the NCIS R VS T pattern variance change-point
	time Rscript $script_folder'gene_atac/ncis_change_point.R' $input_folder'ncis_table_list.txt' $scale_factor_matrix'ncis_table_list.t_thresh.txt' $atac_ncis_table $atac_ncis_table_plot BinSeg

	### get scale factor table
	time Rscript $script_folder'gene_atac/ncis_scale_factor_matrix.R' $scale_factor_matrix'ncis_table_list.t_thresh.txt' $scale_factor_matrix'ncis_table_list.sf.txt' $atac_reads_table $atac_ncis_table_plot 500000 T

	### compare dif norm results
	time Rscript $script_folder'variance_check_atac_rna/check_variance_cor.R'

	### get ncis norm cRE signals
	time Rscript $script_folder'variance_check_atac_rna/norm_matrix.R' '/Volumes/MAC_Data/data/labs/hardison_lab/vision/atacseq/index_set_sig_matrix/celltype.index.5end.sorted.txt' '/Volumes/MAC_Data/data/labs/hardison_lab/vision/gene_rnaseq_atac/scale_factor_matrix/ncis_table_list.sf.txt' '/Volumes/MAC_Data/data/labs/hardison_lab/vision/atacseq/index_set_sig_matrix/celltype.index.5end.ncis_normed.sorted.txt'
	time Rscript $script_folder'variance_check_atac_rna/norm_matrix.R' '/Volumes/MAC_Data/data/labs/hardison_lab/vision/atacseq/index_set_sig_matrix/celltype.index_set.5end.sorted.txt' '/Volumes/MAC_Data/data/labs/hardison_lab/vision/gene_rnaseq_atac/scale_factor_matrix/ncis_table_list.sf.txt' '/Volumes/MAC_Data/data/labs/hardison_lab/vision/atacseq/index_set_sig_matrix/celltype.index_set.5end.ncis_normed.sorted.txt'





