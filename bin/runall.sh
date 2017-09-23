
script_folder='/Volumes/MAC_Data/data/labs/hardison_lab/vision/bin/'
### get binary matrix & signal matrix ###
echo binary matrix and signal matrix
time python $script_folder'split_signal_binary_matrix.py' -i homerTable3.peaks.filtered.txt -n 4 -a homerTable3.peaks.filtered.interval.txt -x 4 -b homerTable3.peaks.filtered.signal.txt -y 32 -c homerTable3.peaks.filtered.binary_pattern.txt -z 60

### get cell type average matrix
time python $script_folder'merge_cell_type_data.py' -i homerTable3.peaks.filtered.binary_pattern.txt -o celltype.binary_pattern.txt
time python $script_folder'merge_cell_type_data.py' -i homerTable3.peaks.filtered.signal.txt -o celltype.signal.txt

### get index sets
time python $script_folder'get_index_set.py' 
