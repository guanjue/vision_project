#PBS -l walltime=80:00:00
#PBS -l nodes=2:ppn=8
#PBS -l mem=100gb

cd /gpfs/home/gzx103/scratch/vision_clustering/gene_atac_sig

cp gencode_pc_sort.bed gencode_pc_sort.atac.txt
cp gencode_pc_sort.exp10kb.bed gencode_pc_sort.exp10kb.atac.txt
cp gencode_pc_sort.TSSexp10kb.bed gencode_pc_sort.TSSexp10kb.atac.txt

rm ongoing.txt
rm total_reads_atac.txt
for bam_file in $(cat bam_file_gene_atac.txt)
do
	echo $bam_file
	echo $bam_file >> ongoing.txt
	cp $bam_file tmp.bam
	### unique mapped reads
	samtools view -q 10 -b tmp.bam > tmp_bam.bam
	samtools sort tmp_bam.bam -o tmp_bam.sort.bam
	samtools index tmp_bam.sort.bam
	samtools view -F 0x4 tmp_bam.sort.bam | cut -f 1 | sort | uniq | wc -l >> total_reads_atac.txt
	#samtools view -F 0x40 bam_files/B_100072_result3.bam | cut -f1 | sort | uniq | wc -l >> total_reads.txt
	### bam to bed file
	bedtools bamtobed -i tmp_bam.sort.bam > test.bed
	### get 5' end and shift (+strand: +4; -strand: -5
	cat test.bed | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1, $2+4, $2+4+1, $4, $5, $6; else print $1, $3-5-1, $3-5, $4, $5, $6}' > test.5end.shifted.bed

	### intersected reads number
	bedtools intersect -a gencode_pc_sort.bed -b test.5end.shifted.bed -wa -c > gencode_pc_sort.bed.readcount.bed
	### extract signal column
	cat gencode_pc_sort.bed.readcount.bed | awk -F '\t' -v OFS='\t' '{print $7}' > gencode_pc_sort.bed.readcount.txt
	### merge to matrix
	paste gencode_pc_sort.atac.txt gencode_pc_sort.bed.readcount.txt > gencode_pc_sort.atac.tmp.txt 
	mv gencode_pc_sort.atac.tmp.txt gencode_pc_sort.atac.txt

	### intersected reads number
	bedtools intersect -a gencode_pc_sort.exp10kb.bed -b test.5end.shifted.bed -wa -c > gencode_pc_sort.exp10kb.bed.readcount.bed
	### extract signal column
	cat gencode_pc_sort.exp10kb.bed.readcount.bed | awk -F '\t' -v OFS='\t' '{print $7}' > gencode_pc_sort.exp10kb.bed.readcount.txt
	### merge to matrix
	paste gencode_pc_sort.exp10kb.atac.txt gencode_pc_sort.exp10kb.bed.readcount.txt > gencode_pc_sort.exp10kb.atac.tmp.txt
	mv gencode_pc_sort.exp10kb.atac.tmp.txt gencode_pc_sort.exp10kb.atac.txt

	### intersected reads number
	bedtools intersect -a gencode_pc_sort.TSSexp10kb.bed -b test.5end.shifted.bed -wa -c > gencode_pc_sort.TSSexp10kb.bed.readcount.bed
	### extract signal column
	cat gencode_pc_sort.TSSexp10kb.bed.readcount.bed | awk -F '\t' -v OFS='\t' '{print $7}' > gencode_pc_sort.TSSexp10kb.bed.readcount.txt
	### merge to matrix
	paste gencode_pc_sort.TSSexp10kb.atac.txt gencode_pc_sort.TSSexp10kb.bed.readcount.txt > gencode_pc_sort.TSSexp10kb.atac.tmp.txt
	mv gencode_pc_sort.TSSexp10kb.atac.tmp.txt gencode_pc_sort.TSSexp10kb.atac.txt

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
