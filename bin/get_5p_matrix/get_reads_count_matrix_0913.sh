#PBS -l walltime=40:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=100gb

cd /gpfs/home/gzx103/scratch/vision_clustering

tail -n+2 homerTable3.peaks.filtered.bed | sort -k1,1 -k2,2n > target.sorted.bed
cp target.sorted.bed reads_count_matrix_5end.txt

rm ongoing.txt
rm total_reads.txt
for bam_file in $(cat bam_file.txt)
do
	echo $bam_file
	echo $bam_file >> ongoing.txt
	cp $bam_file tmp.bam
	### unique mapped reads
	samtools view -q 10 -b tmp.bam > tmp_bam.bam
	samtools sort tmp_bam.bam -o tmp_bam.sort.bam
	samtools index tmp_bam.sort.bam
	cp tmp_bam.sort.bam tmp_bam.sort.cp.bam
	cp tmp_bam.sort.bam.bai tmp_bam.sort.cp.bam.bai
	samtools view -F 0x4 tmp_bam.sort.cp.bam | cut -f 1 | sort | uniq | wc -l >> total_reads.txt
	#samtools view -F 0x40 bam_files/B_100072_result3.bam | cut -f1 | sort | uniq | wc -l >> total_reads.txt
	### bam to bed file
	bedtools bamtobed -i tmp_bam.sort.cp.bam > test.bed
	### get 5' end and shift (+strand: +4; -strand: -5
	cat test.bed | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1, $2+4, $2+4+1, $4, $5, $6; else print $1, $3-5-1, $3-5, $4, $5, $6}' > $bam_file'test.5end.shifted.bed'
	### intersected reads number
	bedtools intersect -a target.sorted.bed -b $bam_file'test.5end.shifted.bed' -wa -c > tmp_readcount.bed
	#bedtools multicov -bams tmp_bam.sort.bam -bed target.sorted.bed > tmp_readcount.bed
	cat tmp_readcount.bed | awk -F '\t' -v OFS='\t' '{print $6}' > tmp_readcount.txt
	paste reads_count_matrix_5end.txt tmp_readcount.txt > reads_count_matrix_tmp.txt
	mv reads_count_matrix_tmp.txt reads_count_matrix_5end.txt

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
