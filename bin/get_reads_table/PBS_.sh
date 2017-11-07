#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=40gb

cd /gpfs/home/gzx103/scratch/index_caller/test_data
### get reads
#time python reads_count_reads.py -i reads_count_matrix_5end_whole_genome.add_id.txt -o reads_count_matrix_5end_whole_genome_reads.txt
time python reads_count_reads_count.py -i reads_count_matrix_5end_whole_genome.add_id.txt -o reads_count_matrix_5end_whole_genome_reads.txt
#time python cell_type_sort.py -i reads_count_matrix_5end_whole_genome_reads.txt -b bam_file.txt -r homerTable3.peaks.filtered.txt -o homerTable3.peaks.filtered.reads.txt
time python cell_type_sort.py -i reads_count_matrix_5end_whole_genome_reads.txt -b bam_file.txt -r homerTable3.peaks.filtered.txt -o homerTable3.peaks.filtered.reads.txt
#time python merge_cell_type_data.py -i homerTable3.peaks.filtered.reads.txt -m sample2celltype.txt -n 2 -o celltype.reads.txt
time python merge_cell_type_data.py -i homerTable3.peaks.filtered.reads.txt -m sample2celltype.txt -n 2 -o celltype.reads.txt


cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$2}' > lsk.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$3}' > cmp.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$4}' > gmp.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$5}' > mep.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$6}' > cfue.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$7}' > ery.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$8}' > cfumk.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$9}' > meg.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$10}' > mono.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$11}' > neu.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$12}' > b.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$13}' > nk.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$14}' > tcd4.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$15}' > tcd8.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$16}' > g1e.wg.reads.txt
cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$17}' > er4.wg.reads.txt


var=1
while read LINE
do
	output_name=$(echo "$LINE" | awk '{print $1}')
	echo $output_name
	var=$((var+1))
	cat homerTable3.peaks.filtered.reads.txt | awk -v var="$var" '{print $1, $var}' > $output_name'.wg.reads.txt'
done < homerTable3.peaks.filtered.txt_column_name.txt

time tmp.sh

#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$2}' > lsk.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$3}' > cmp.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$4}' > gmp.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$5}' > mep.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$6}' > cfue.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$7}' > ery.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$8}' > cfumk.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$9}' > meg.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$10}' > mono.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$11}' > neu.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$12}' > b.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$13}' > nk.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$14}' > tcd4.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$15}' > tcd8.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$16}' > g1e.wg.reads.txt
#cat celltype.reads.txt | awk -F '\t' -v OFS='\t' '{print $1,$17}' > er4.wg.reads.txt



