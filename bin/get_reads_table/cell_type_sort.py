import os
import numpy as np
################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return data0
################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()

################################################################################################

def cell_type_sort(input_matrix, bam_file_order, given_col_order, tpm_celltype_sorted_matrix):
	### read matrix
	data = read2d_array(input_matrix, str)

	### split row info & row data
	data_info = np.transpose([data[:,3]])
	data_sig = np.array(data[:,4:], dtype = float)
	#print(data_sig.shape)

	### read bam files (when generate the input matrix) orders
	bam_files = open(bam_file_order, 'r')

	bam_order = []
	for records in bam_files:
		header_tmp = records.split('/')[1].split('_')
		bam_order.append(header_tmp[0] + '_' + header_tmp[1])
	bam_order = np.array(bam_order)
	#print(bam_order)

	### data matrix to data dict
	data_sig_dict = {}
	for i in range(0, data_sig.shape[1]):
		#print(str(i)+': '+bam_order[i])
		data_sig_dict[bam_order[i]] = data_sig[:,i]

	### get given column order
	given_order = open(given_col_order, 'r')	
	given_order = given_order.readline().split('\t')[4:32]

	### sort column by given col order in dict
	data_col_sorted = data_info
	for j in range(0, len(given_order)):
		celltype_col = np.transpose([data_sig_dict[given_order[j]]])
		data_col_sorted = np.concatenate((data_col_sorted, celltype_col), axis = 1)

	### write header
	result = open(tpm_celltype_sorted_matrix, 'w')
	#result.write('chrom'+'\t'+'start'+'\t'+'end'+'\t'+'name'+'\t')
	result.write('name'+'\t')
	for j in range(0, len(given_order)-1):
		result.write(given_order[j]+'\t')
	result.write(given_order[len(given_order)-1]+'\n')

	### write info
	for records in data_col_sorted:
		for i in range(0, len(records)-1):
			result.write(records[i]+'\t')
		result.write(records[len(records)-1]+'\n')

	result.close()

############################################################################
#time python cell_type_sort.py -i reads_count_matrix_5end_tpm.txt -b bam_file.txt -r homerTable3.peaks.filtered.txt -o homerTable3.peaks.filtered.tpm.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:b:r:o:")
	except getopt.GetoptError:
		print 'time python cell_type_sort.py -i input_tpm_matrix -b bam_file_list -r target_col_order -o input_tpm_colsorted_matrix'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python cell_type_sort.py -i input_tpm_matrix -b bam_file_list -r target_col_order -o input_tpm_colsorted_matrix'
			sys.exit()
		elif opt=="-i":
			input_matrix=str(arg.strip())
		elif opt=="-b":
			bam_file_order=str(arg.strip())
		elif opt=="-r":
			given_col_order=str(arg.strip())
		elif opt=="-o":
			tpm_celltype_sorted_matrix=str(arg.strip())

	cell_type_sort(input_matrix, bam_file_order, given_col_order, tpm_celltype_sorted_matrix)

if __name__=="__main__":
	main(sys.argv[1:])