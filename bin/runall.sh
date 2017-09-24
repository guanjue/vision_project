##################################

##################################
### index set module
script_folder='/Volumes/MAC_Data/data/labs/hardison_lab/vision/bin/'
### get binary matrix & signal matrix ###
echo binary matrix and signal matrix
time python $script_folder'index_set/split_signal_binary_matrix.py' -i homerTable3.peaks.filtered.txt -n 4 -a homerTable3.peaks.filtered.interval.txt -x 4 -b homerTable3.peaks.filtered.signal.txt -y 32 -c homerTable3.peaks.filtered.binary_pattern.txt -z 60

### get cell type average matrix
time python $script_folder'index_set/merge_cell_type_data.py' -i homerTable3.peaks.filtered.binary_pattern.txt -o celltype.binary_pattern.txt
time python $script_folder'index_set/merge_cell_type_data.py' -i homerTable3.peaks.filtered.signal.txt -o celltype.signal.txt

### get index sets
time python $script_folder'index_set/get_index_set.py' -i celltype.binary_pattern.txt -r celltype.order.txt -l signal_level_range.txt -f celltype.index.sorted.txt -s celltype.index_set.sorted.txt

### plot index set
time Rscript $script_folder'figures/plot_index_set_region_hist.R' celltype.index_set.sorted.txt
time Rscript $script_folder'figures/plot_index_set_module.R' celltype.index_set.sorted.txt celltype.index.sorted.txt black 200 index_set_all.pdf index_set_thresh.pdf index.png black

##################################
### signal matrix sort
time python $script_folder'sort_matrix/get_index_signal_matrix.py' -t celltype.signal.txt -a 1 -s celltype.index.sorted.txt -b 1 -r celltype.order.txt -q 75 -o celltype.index.signal.sorted.txt -p celltype.index_set.signal.sorted.txt

### filter the original matrix to 210k matrix
time python $script_folder'sort_matrix/vlookup.py' -t homerTable3.peaks.filtered.txt -m 4 -s celltype.index.sorted.txt -n 1 -o homerTable3.peaks.filtered.210k.txt -k F

### plot index set signal
time Rscript $script_folder'figures/plot_index_set_signal_module.R' celltype.index_set.signal.sorted.txt celltype.index.signal.sorted.txt black 200 index_set_signal_all.pdf index_set_signal_thresh.pdf index_signal.png red celltype.index_set.sorted.txt

##################################
### enriched IDEAS states




