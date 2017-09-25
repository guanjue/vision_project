##################################

##################################
### index set module
script_folder='/Volumes/MAC_Data/data/labs/hardison_lab/vision/bin/'

### get 5p end reads count per DNA region
echo 'get 5p end reads count per DNA region'
# run in server (biostar)
# time bash get_reads_count_matrix_0913.sh
# bam file location: bam_file.txt
# peak file: homerTable3.peaks.filtered.bed
# output read count matrix: reads_count_matrix_5end.txt

######## homer TPM
### prepare TPM matrix
echo 'calculate tpm'
time python $script_folder'get_5p_matrix/reads_count_tpm.py' -i reads_count_matrix_5end.txt -o reads_count_matrix_5end_tpm.txt

### cell type sorting (column sorting)
echo 'sort tpm column order (come from bam file name sort order) by homerTable3.peaks.filtered.txt column order'
time python $script_folder'get_5p_matrix/cell_type_sort.py' -i reads_count_matrix_5end_tpm.txt -b bam_file.txt -r homerTable3.peaks.filtered.txt -o homerTable3.peaks.filtered.tpm.txt

######## homer RPKM
### prepare RPKM matrix
echo 'calculate tpm'
time python $script_folder'get_5p_matrix/reads_count_rpkm.py' -i reads_count_matrix_5end.txt -t total_reads_all.txt -o reads_count_matrix_5end_rpkm.txt

### cell type sorting (column sorting)
echo 'sort tpm column order (come from bam file name sort order) by homerTable3.peaks.filtered.txt column order'
time python $script_folder'get_5p_matrix/cell_type_sort.py' -i reads_count_matrix_5end_rpkm.txt -b bam_file.txt -r homerTable3.peaks.filtered.txt -o homerTable3.peaks.filtered.rpkm.txt
