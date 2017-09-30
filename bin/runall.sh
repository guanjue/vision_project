##################################
script_folder='/Volumes/MAC_Data/data/labs/hardison_lab/vision/bin/'
index_caller_script_folder='/Volumes/MAC_Data/data/labs/zhang_lab/01projects/index_caller/index_caller/bin/'

##################################
	###### initiate folders
	input_folder='input_data/'
	index_set_dir='index_set_matrix/'
	index_set_sig_dir='index_set_sig_matrix/'
	index_set_ideas_RE_dir='index_set_ideas_RE/'
	index_set_figure_dir='index_set_figure/'

	### mkdir index_set module folder
	if [ -d "$index_set_dir" ]; then  
    	rm -r $index_set_dir
	fi
	mkdir $index_set_dir

	### mkdir index_set module folder
	if [ -d "$index_set_sig_dir" ]; then  
    	rm -r $index_set_sig_dir
	fi
	mkdir $index_set_sig_dir

	### mkdir ideas state output folder
	if [ -d "$index_set_ideas_RE_dir" ]; then  
    	rm -r $index_set_ideas_RE_dir
	fi
	mkdir $index_set_ideas_RE_dir

	### mkdir figure folder
	if [ -d "$index_set_figure_dir" ]; then  
    	rm -r $index_set_figure_dir
	fi
	mkdir $index_set_figure_dir

##################################
	######## index set module
	### get binary matrix & signal matrix ###
	echo 'binary matrix and signal matrix'
	time python $script_folder'index_set/split_signal_binary_matrix.py' -i $input_folder'homerTable3.peaks.filtered.txt' -n 4 -a $index_set_dir'homerTable3.peaks.filtered.interval.txt' -x 4 -b $index_set_dir'homerTable3.peaks.filtered.signal.txt' -y 32 -c $index_set_dir'homerTable3.peaks.filtered.binary_pattern.txt' -z 60

	### get cell type average matrix
	echo 'get cell type average matrix'
	time python $script_folder'index_set/merge_cell_type_data.py' -i $index_set_dir'homerTable3.peaks.filtered.binary_pattern.txt' -m $input_folder'sample2celltype.txt' -n 2 -o $index_set_dir'celltype.binary_pattern.txt'
	time python $script_folder'index_set/merge_cell_type_data.py' -i $index_set_dir'homerTable3.peaks.filtered.signal.txt' -m $input_folder'sample2celltype.txt' -n 2 -o $index_set_dir'celltype.signal.txt'

	### get index sets (CORE!!!)
	echo 'get index sets (CORE!!!)'
	time python $script_folder'index_set/get_index_set.py' -i $index_set_dir'celltype.binary_pattern.txt' -r $input_folder'celltype.order.txt' -l $input_folder'signal_level_range.txt' -f $index_set_dir'celltype.index.sorted.txt' -s $index_set_dir'celltype.index_set.sorted.txt'

	### plot index set
	echo 'plot index set'
	time Rscript $script_folder'figures/plot_index_set_region_hist.R' $index_set_dir'celltype.index_set.sorted.txt' $index_set_figure_dir'index_hist.pdf' $index_set_figure_dir'index_hist_noxlim.pdf'
	time Rscript $script_folder'figures/plot_index_set_module.R' $index_set_dir'celltype.index_set.sorted.txt' $index_set_dir'celltype.index.sorted.txt' $index_set_dir'celltype.index_set_filtered.sorted.txt' $index_set_figure_dir'index_set_all.pdf' $index_set_figure_dir'index_set_thresh.pdf' $index_set_figure_dir'index.png' 200 black

##################################
### get input signal (tpm & rpkm) matrix
	echo 'get input signal (tpm & rpkm) matrix'
	#time bash run_get_input_matrices.sh

##################################
	######## homer signal
	### signal matrix sort
	echo 'signal matrix sort'
	time python $script_folder'sort_matrix/get_index_signal_matrix.py' -t $index_set_dir'celltype.signal.txt' -a 1 -s $index_set_dir'celltype.index.sorted.txt' -b 1 -r $input_folder'celltype.order.txt' -q 75 -o $index_set_sig_dir'celltype.index.signal.sorted.txt' -p $index_set_sig_dir'celltype.index_set.signal.sorted.txt'
	### filter the original matrix to 210k matrix
	echo 'filter the original matrix to 210k matrix'
	time python $script_folder'sort_matrix/vlookup.py' -t $input_folder'homerTable3.peaks.filtered.txt' -m 4 -s $index_set_dir'celltype.index.sorted.txt' -n 1 -o $input_folder'homerTable3.peaks.filtered.210k.txt' -k F
	### plot index set signal
	echo 'plot index set signal'
	time Rscript $script_folder'figures/plot_index_set_signal_module.R' $index_set_sig_dir'celltype.index_set.signal.sorted.txt' $index_set_sig_dir'celltype.index.signal.sorted.txt' $index_set_dir'celltype.index_set.sorted.txt' $index_set_figure_dir'index_set_signal_all.pdf' $index_set_figure_dir'index_set_signal_thresh.pdf' $index_set_figure_dir'index_signal.png' red 200 0.99 log2

	######## prepare signal matrix
	echo 'get IDEAS states matrix (210k)'
	#time bash run_get_input_matrices.sh

	######## DNA region TPM
	echo 'get cell type average matrix'
	time python $script_folder'index_set/merge_cell_type_data.py' -i $input_folder'homerTable3.peaks.filtered.tpm.txt' -o $index_set_dir'celltype.tpm.txt'
	### tpm matrix sort
	echo 'tpm matrix sort'
	time python $script_folder'sort_matrix/get_index_signal_matrix.py' -t $index_set_dir'celltype.tpm.txt' -a 1 -s $index_set_dir'celltype.index.sorted.txt' -b 1 -r $input_folder'celltype.order.txt' -q 75 -o $index_set_sig_dir'celltype.index.tpm.sorted.txt' -p $index_set_sig_dir'celltype.index_set.tpm.sorted.txt'
	### filter the original matrix to 210k matrix
	echo 'filter the original matrix to 210k matrix'
	time python $script_folder'sort_matrix/vlookup.py' -t $input_folder'homerTable3.peaks.filtered.tpm.txt' -m 1 -s $index_set_dir'celltype.index.sorted.txt' -n 1 -o $input_folder'homerTable3.peaks.filtered.tpm.210k.txt' -k F
	### plot index set TPM
	echo 'plot index set TPM'
	time Rscript $script_folder'figures/plot_index_set_signal_module.R' $index_set_sig_dir'celltype.index_set.tpm.sorted.txt' $index_set_sig_dir'celltype.index.tpm.sorted.txt' $index_set_dir'celltype.index_set.sorted.txt' $index_set_figure_dir'index_set_tpm_all.pdf' $index_set_figure_dir'index_set_tpm_thresh.pdf' $index_set_figure_dir'index_tpm.png' red 200 0.99 log2

	######## DNA region RPKM
	echo 'get cell type average matrix'
	time python $script_folder'index_set/merge_cell_type_data.py' -i $input_folder'homerTable3.peaks.filtered.rpkm.txt' -o $index_set_dir'celltype.rpkm.txt'
	### RPKM matrix sort
	echo 'RPKM matrix sort'
	time python $script_folder'sort_matrix/get_index_signal_matrix.py' -t $index_set_dir'celltype.rpkm.txt' -a 1 -s $index_set_dir'celltype.index.sorted.txt' -b 1 -r $input_folder'celltype.order.txt' -q 75 -o $index_set_sig_dir'celltype.index.rpkm.sorted.txt' -p $index_set_sig_dir'celltype.index_set.rpkm.sorted.txt'
	### filter the original matrix to 210k matrix
	echo 'filter the original matrix to 210k matrix'
	time python $script_folder'sort_matrix/vlookup.py' -t $input_folder'homerTable3.peaks.filtered.rpkm.txt' -m 1 -s $index_set_dir'celltype.index.sorted.txt' -n 1 -o $input_folder'homerTable3.peaks.filtered.rpkm.210k.txt' -k F
	### plot index set RPKM
	echo 'plot index set RPKM'
	time Rscript $script_folder'figures/plot_index_set_signal_module.R' $index_set_sig_dir'celltype.index_set.rpkm.sorted.txt' $index_set_sig_dir'celltype.index.rpkm.sorted.txt' $index_set_dir'celltype.index_set.sorted.txt' $index_set_figure_dir'index_set_rpkm_all.pdf' $index_set_figure_dir'index_set_rpkm_thresh.pdf' $index_set_figure_dir'index_rpkm.png' red 200 0.99 log2

##################################
### get IDEAS states matrix
	echo 'get IDEAS states matrix (210k)'
	#time bash run_get_RE_associate_matrices.sh

	### get index set most Frequent ideas state matrix (the output with $index_set_dir'celltype.index.sorted.txt' row order)
	time python $script_folder'index_set_RE_freq_matrix/most_freq_matrix.py' -t 'ideas_bb/DNA_regin_210k.celltype_sorted.txt' -a 1 -s $index_set_dir'celltype.index.sorted.txt' -b 1 -r $input_folder'celltype.order.txt' -o $index_set_ideas_RE_dir'celltype.index.ideas_RE.sorted.txt' -p $index_set_ideas_RE_dir'celltype.index_set.ideas_RE.sorted' -n 17

	### plot index set most Frequent ideas state matrix
	echo 'plot index set most Frequent ideas state matrix'
	time Rscript $script_folder'figures/plot_index_set_ideas_state_module.R' $index_set_ideas_RE_dir'celltype.index_set.ideas_RE.sorted.freq_state.txt' $index_set_ideas_RE_dir'celltype.index_set.ideas_RE.sorted.state.se.txt' $index_set_ideas_RE_dir'celltype.index.ideas_RE.sorted.txt' $index_set_dir'celltype.index_set.sorted.txt' $index_set_figure_dir'index_set_ideas_RE_all.pdf' $index_set_figure_dir'index_set_ideas_RE_thresh.pdf' $index_set_figure_dir'index_set_ideas_RE_SE_all.pdf' $index_set_figure_dir'index_set_ideas_RE_SE_thresh.pdf' $index_set_figure_dir'index_ideas_RE.png' $input_folder'state_color.txt' 200 0.99 log2

##################################
### specific analysis
	### run PCA for each sample based on 210k cRES
	time Rscript $script_folder'figures/sample_cor_and_PCA.R' $input_folder'homerTable3.peaks.filtered.210k.txt' $index_set_figure_dir'homerTable3.peaks.filtered.210k' 5 spearman
	time Rscript $script_folder'figures/sample_cor_and_PCA.R' $input_folder'homerTable3.peaks.filtered.tpm.210k.txt' $index_set_figure_dir'homerTable3.peaks.filtered.tpm.210k' 2 spearman
	time Rscript $script_folder'figures/sample_cor_and_PCA.R' $input_folder'homerTable3.peaks.filtered.rpkm.210k.txt' $index_set_figure_dir'homerTable3.peaks.filtered.rpkm.210k' 2 spearman

	### select cREs with at least 1 active ideas state in 16 cell types
	time python $script_folder'index_set_RE_freq_matrix/select_cREs_by_ideas_state.py' -t $index_set_ideas_RE_dir'celltype.index.ideas_RE.sorted.txt' -a $input_folder'ideas_active_state_label.txt' -o $index_set_ideas_RE_dir'celltype.index.sorted.active_state_filtered.txt'

	### select cREs with at atac-seq only ideas state in all 16 cell types
	time python $script_folder'index_set_RE_freq_matrix/select_cREs_by_all_one_ideas_state.py' -t $index_set_ideas_RE_dir'celltype.index.ideas_RE.sorted.txt' -a $input_folder'atac_only_ideas_state_label.txt' -o $index_set_ideas_RE_dir'celltype.index.sorted.atac_only_state_filtered.txt' -s 17
	### select cREs with at ctcf only ideas state in all 16 cell types
	time python $script_folder'index_set_RE_freq_matrix/select_cREs_by_all_one_ideas_state.py' -t $index_set_ideas_RE_dir'celltype.index.ideas_RE.sorted.txt' -a $input_folder'ctcf_only_ideas_state_label.txt' -o $index_set_ideas_RE_dir'celltype.index.sorted.ctcf_only_state_filtered.txt' -s 17

	### plot cRE number for each cell type 
	time Rscript $script_folder'figures/plot_celltype_cRE_hist.R' $index_set_dir'celltype.index.sorted.txt' $index_set_figure_dir'cRE_celltype_number.pdf'

	### plot the proportion of each ideas state in cRE for each cell type
	time Rscript $script_folder'figures/plot_celltype_cRE_IDEASpro_bar.R' $index_set_ideas_RE_dir'celltype.index.ideas_RE.sorted.txt' $input_folder'state_color.txt' $index_set_figure_dir'celltype_cRE_IDEASpro_bar.pdf'







