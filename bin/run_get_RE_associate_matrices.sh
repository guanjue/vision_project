##################################

##################################
### index set module
script_folder='/Volumes/MAC_Data/data/labs/hardison_lab/vision/bin/'

cd /Volumes/MAC_Data/data/labs/hardison_lab/vision/data_test
### get 210k DNA region bed file (interval sorted)
tail -n+2 'input_data/homerTable3.peaks.filtered.interval.210k.txt' | awk -F '\t' -v OFS='\t' '{print $2, $3, $4, $5}' | sort -k1,1 -k2,2n > 'input_data/DNA_intervals_210k.bed'
### get 210k DNA region bed file (index sorted)
tail -n+2 'input_data/homerTable3.peaks.filtered.interval.210k.txt' > 'ideas_bb/DNA_regin_210k_indexsort.bed'
### get 210k DNA region bed file (index sorted) only intervals
tail -n+2 'input_data/homerTable3.peaks.filtered.interval.210k.txt' > 'ideas_bb/DNA_regin_210k_indexsort_onlyinterval.bed'


while read LINE
do
	celltype_id_num=$(echo "$LINE" | awk '{print $1}')
	celltype_name=$(echo "$LINE" | awk '{print $2}')
  	echo $celltype_id_num
	echo $celltype_name

	### convert bigbed to bed file
	#$script_folder'index_set_RE_freq_matrix/bigBedToBed' 'ideas_bb/ideasVision'$celltype_id_num'.bb' 'ideas_bb/ideasVision'$celltype_id_num'.bed'
	#cat 'ideas_bb/ideasVision'$celltype_id_num'.bed' | sort -k1,1 -k2,2n > 'ideas_bb/ideasVision'$celltype_id_num'.sort.bed'
	#rm 'ideasVision'$celltype_id_num'.bed'
	
	### intersect ideas state
	bedtools window -a 'input_data/DNA_intervals_210k.bed' -b 'ideas_bb/ideasVision'$celltype_id_num'.sort.bed' -w 0 > 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.tmp.bed'
	
	### get input
	cat 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.tmp.bed' | awk -F '\t' -v OFS='\t' '{
		if ($6-$2>=0 && $7-$3<=0) print $1,$2,$3,$4,$7-$6,$8,$13,($7+$6-$3-$2)/2,$7-$6; 
		else if ($6-$2<0 && $7-$3>0) print $1,$2,$3,$4,$3-$2,$8,$13,($7+$6-$3-$2)/2,$7-$6;
		else if ($6-$2<0 && $7-$3<=0) print $1,$2,$3,$4,$7-$2,$8,$13,($7+$6-$3-$2)/2,$7-$6;
		else if ($6-$2>=0 && $7-$3>0) print $1,$2,$3,$4,$3-$6,$8,$13,($7+$6-$3-$2)/2,$7-$6
		}' > 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.bed'

	### get bed file with color
	python $script_folder'index_set_RE_freq_matrix/get_cRE_ideas_state.py' -a 'input_data/homerTable3.peaks.filtered.interval.210k.txt' -i 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.bed' -o 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.colored.bed'

	### add label names
	time python $script_folder'index_set_RE_freq_matrix/add_label_name.py' -i 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.colored.bed' -l 'input_data/ideas_state_label.txt' -o 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.colored.named.bed'

	### get ideas state matrix
	cat 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.colored.named.bed' | awk -F '\t' -v OFS='\t' '{print $4}' > 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.colored.named.txt'
	paste 'ideas_bb/DNA_regin_210k_indexsort.bed' 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.colored.named.txt' > 'ideas_bb/DNA_regin_210k.tmp.bed' && mv 'ideas_bb/DNA_regin_210k.tmp.bed' 'ideas_bb/DNA_regin_210k_indexsort.bed'

	### rm tmp files
	rm 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.colored.bed'
	rm 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.bed'
	rm 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.tmp.bed'
	rm 'ideas_bb/DNA_intervals_210k.'$celltype_id_num'.colored.named.txt'

done < input_data/ideas_cell_order.txt

### sort ideas state matrix cell type order (column order)
time python $script_folder'index_set_RE_freq_matrix/ideas_state_cell_type_sort.py' -i 'ideas_bb/DNA_regin_210k_indexsort.bed' -b 'input_data/ideas_cell_order.txt' -r 'input_data/homerTable3.peaks.filtered.txt' -o 'ideas_bb/DNA_regin_210k.celltype_sorted.txt' -n 17


