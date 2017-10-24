#PBS -l walltime=40:00:00
#PBS -l nodes=2:ppn=8
#PBS -l mem=100gb

cd /gpfs/home/gzx103/scratch/vision_clustering

#tail -n+2 target.whole_genome.bed | sort -k1,1 -k2,2n > target.sorted_whole_genome.bed
cp target.sorted_whole_genome.bed reads_count_matrix_wg_5end.txt

rm ongoing.txt
rm total_reads.txt
for bam_file in $(cat bam_file.txt)
do
	echo $bam_file
	echo $bam_file >> ongoing.txt
	### intersected reads number
	bedtools intersect -a target.sorted_whole_genome.bed -b $bam_file'test.5end.shifted.bed' -wa -c > tmp_readcount_wg.bed

	cat tmp_readcount_wg.bed | awk -F '\t' -v OFS='\t' '{print $6}' > tmp_readcount_wg.txt
	paste reads_count_matrix_wg_5end.txt tmp_readcount_wg.txt > reads_count_matrix_wg_tmp.txt
	mv reads_count_matrix_wg_tmp.txt reads_count_matrix_wg_5end.txt
done


shuf -n 5000000 reads_count_matrix_wg_5end.txt > reads_count_matrix_wg_5end_sample.txt

