#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -A open
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:[add paths to libraries if not found];

module load gcc
module load bedtools
module load samtools

work_folder='/storage/home/gzx103/scratch/vision_clustering/histones/h3k4me3/'
input_bam_folder='/storage/home/gzx103/group/projects/vision/histone/h3k4me3/bam/unsorted_bam/'
tmp_output_folder='/storage/home/gzx103/scratch/vision_clustering/histones/h3k4me3/tmp_files/'
matrix_folder='/storage/home/gzx103/scratch/vision_clustering/histones/h3k4me3/matrix/'

cd $work_folder

rm $matrix_folder'ongoing.txt'
rm $matrix_folder'total_reads_q30.txt'

tail -n+2 $matrix_folder'200_noblack.windows' | sort -k1,1 -k2,2n > $matrix_folder'target.sorted.bed'
cp $matrix_folder'target.sorted.bed' $matrix_folder'h3k4me3_reads_count_matrix_5end.txt'

for bam_file in $(cat $matrix_folder'h3k4me3_bamfile_list.txt')
do
	echo $bam_file
	echo $bam_file >> $tmp_output_folder'ongoing.txt'
	### sort bam
	samtools sort $input_bam_folder$bam_file $tmp_output_folder$bam_file'.sorted'

	### index bam file
	samtools index $tmp_output_folder$bam_file'.sorted.bam'

	### only select mapping quality greater than 30
	samtools view -q 30 -b $tmp_output_folder$bam_file'.sorted.bam' > $tmp_output_folder$bam_file'.sorted.q30.bam'

	### index q30 bam file
	samtools index $tmp_output_folder$bam_file'.sorted.q30.bam'

	### get total reads in the q30 bam file
	samtools view -F 0x4 $tmp_output_folder$bam_file'.sorted.q30.bam' | cut -f 1 | sort | uniq | wc -l >> $matrix_folder'total_reads_q30.txt'

	### bam to bed file
	bedtools bamtobed -i $tmp_output_folder$bam_file'.sorted.q30.bam' > $tmp_output_folder$bam_file'.sorted.q30.bam.bed'

	### get 5' end and shift (+strand: +0; -strand: -0)
	cat $tmp_output_folder$bam_file'.sorted.q30.bam.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1, $2+0, $2+0+1, $4, $5, $6; else print $1, $3-0-1, $3-0, $4, $5, $6}' > $tmp_output_folder$bam_file'.sorted.q30.5end.shifted.bam.bed'

	### intersected reads number
	bedtools intersect -a $matrix_folder'target.sorted.bed' -b $tmp_output_folder$bam_file'.sorted.q30.5end.shifted.bam.bed' -wa -c > $tmp_output_folder$bam_file'.sorted.q30.5end.shifted.bam.readscount.bed'
	
	### get the signal column
	cat $tmp_output_folder$bam_file'.sorted.q30.5end.shifted.bam.readscount.bed' | awk -F '\t' -v OFS='\t' '{print $4}' > $tmp_output_folder$bam_file'.sorted.q30.5end.shifted.bam.readscount.txt'
	
	### cbind to reads count matrix
	paste $matrix_folder'h3k4me3_reads_count_matrix_5end.txt' $tmp_output_folder$bam_file'.sorted.q30.5end.shifted.bam.readscount.txt' > $matrix_folder'h3k4me3_reads_count_matrix_5end_tmp.txt'
	
	### change names
	mv $matrix_folder'h3k4me3_reads_count_matrix_5end_tmp.txt' $matrix_folder'h3k4me3_reads_count_matrix_5end.txt'

	#rm tmp.bam 
	#rm tmp_bam.bam
	#rm tmp_bam.sort.bam
	#rm tmp_readcount.bed
	#rm tmp_readcount.txt
	#rm tmp_bam.sort.bam.bai
	#rm test.bed
	#rm test.5end.shifted.bed
	#rm tmp_readcount.bed
done



