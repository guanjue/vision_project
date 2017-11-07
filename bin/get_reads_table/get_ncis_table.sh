#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=40gb

cd /gpfs/home/gzx103/scratch/index_caller/test_data

time python get_atac_t_ncis_fixwin.py -i lsk.wg.reads.txt -j cmp.wg.reads.txt -o lsk_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i cmp.wg.reads.txt -j cmp.wg.reads.txt -o cmp_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i gmp.wg.reads.txt -j cmp.wg.reads.txt -o gmp_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i mep.wg.reads.txt -j cmp.wg.reads.txt -o mep_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i cfumk.wg.reads.txt -j cmp.wg.reads.txt -o cfumk_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i meg.wg.reads.txt -j cmp.wg.reads.txt -o meg_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i mono.wg.reads.txt -j cmp.wg.reads.txt -o mono_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i cfue.wg.reads.txt -j cmp.wg.reads.txt -o cfue_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i ery.wg.reads.txt -j cmp.wg.reads.txt -o ery_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i neu.wg.reads.txt -j cmp.wg.reads.txt -o neu_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i g1e.wg.reads.txt -j cmp.wg.reads.txt -o g1e_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i er4.wg.reads.txt -j cmp.wg.reads.txt -o er4_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i b.wg.reads.txt -j cmp.wg.reads.txt -o b_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i nk.wg.reads.txt -j cmp.wg.reads.txt -o nk_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i tcd4.wg.reads.txt -j cmp.wg.reads.txt -o tcd4_vs_cmp.200.txt
time python get_atac_t_ncis_fixwin.py -i tcd8.wg.reads.txt -j cmp.wg.reads.txt -o tcd8_vs_cmp.200.txt

