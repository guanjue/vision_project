##################################

##################################
### index set module
script_folder='/Volumes/MAC_Data/data/labs/hardison_lab/vision/bin/'

### get 5p end reads count per DNA region
echo 'get 5p end reads count per DNA region'
# run in server (biostar)
# time bash /gpfs/home/gzx103/scratch/vision_clustering/get_reads_count_matrix_0913.sh
# time bash get_5p_matrix/get_reads_count_matrix_0913.sh

# bam file location: /gpfs/home/gzx103/scratch/vision_clustering/bam_file.txt
# bam file location: list_file_example/bam_file.txt

# peak file: homerTable3.peaks.filtered.bed
# output read count matrix: reads_count_matrix_5end.txt

######## homer TPM
### prepare TPM matrix
echo 'calculate tpm'
time python $script_folder'get_5p_matrix/reads_count_tpm.py' -i 'input_data/reads_count_matrix_5end.txt' -o 'input_data/reads_count_matrix_5end_tpm.txt'
time python $script_folder'get_5p_matrix/reads_count_tpm.py' -i 'input_data/reads_count_matrix_5end_whole_genome.txt' -o 'input_data/reads_count_matrix_5end_whole_genome_tpm.txt'

### cell type sorting (column sorting)
echo 'sort tpm column order (come from bam file name sort order) by homerTable3.peaks.filtered.txt column order'
time python $script_folder'get_5p_matrix/cell_type_sort.py' -i 'input_data/reads_count_matrix_5end_tpm.txt' -b 'input_data/bam_file.txt' -r 'input_data/homerTable3.peaks.filtered.txt' -o 'input_data/homerTable3.peaks.filtered.tpm.txt'

######## homer RPKM
### prepare RPKM matrix
echo 'calculate tpm'
time python $script_folder'get_5p_matrix/reads_count_rpkm.py' -i 'input_data/reads_count_matrix_5end.txt' -t 'input_data/total_reads_all.txt' -o 'input_data/reads_count_matrix_5end_rpkm.txt'

### cell type sorting (column sorting)
echo 'sort tpm column order (come from bam file name sort order) by homerTable3.peaks.filtered.txt column order'
time python $script_folder'get_5p_matrix/cell_type_sort.py' -i 'input_data/reads_count_matrix_5end_rpkm.txt' -b 'input_data/bam_file.txt' -r 'input_data/homerTable3.peaks.filtered.txt' -o 'input_data/homerTable3.peaks.filtered.rpkm.txt'


