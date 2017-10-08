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
### read 2d list
def read2d_list(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		data0.append(tmp)
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
### sample2matrix
def sample2matrix(signal_matrix, col_header, sample2celltype_list):
	### initiate header vector & signal matrix
	cellytpe_merge_matrix_header = []
	cellytpe_merge_matrix_data = np.array([], dtype=np.float).reshape(0, signal_matrix.shape[0])

	### merge sample based on sample2celltype_list
	for i in range(0, len(sample2celltype_list)):
		### column header
		tmp_header = col_header[i]

		### get sample number
		sample_size = len(sample2celltype_list[i]) - 1

		### average column signal (ith cell type, first sample)
		tmp_sig = signal_matrix[ :,int(sample2celltype_list[i][1]) ]
		### (ith cell type, rest samples)
		for j in range(1,sample_size):
			tmp_sig = tmp_sig + signal_matrix[:,int(sample2celltype_list[i][1+j])]
		tmp_sig_average = tmp_sig / sample_size

		### append header vector
		cellytpe_merge_matrix_header = cellytpe_merge_matrix_header + [ tmp_header ]
		### append signal data matrix
		cellytpe_merge_matrix_data = np.concatenate((cellytpe_merge_matrix_data, [tmp_sig_average]), axis=0)

	### transpose data matrix
	cellytpe_merge_matrix_data = np.transpose(cellytpe_merge_matrix_data)
	### reshape header 
	cellytpe_merge_matrix_header = np.array([cellytpe_merge_matrix_header])

	### add column header to data matrix
	cellytpe_merge_matrix_info = np.concatenate((cellytpe_merge_matrix_header, cellytpe_merge_matrix_data), axis=0)
	
	return cellytpe_merge_matrix_info

################################################################################################
def merge_cell_type_data(inputfile, sig_col, sample2celltype, outputfile):
	data0 = read2d_array(inputfile, str)

	### get row info & add header
	data0_info = np.array([range(0, sig_col-1)])
	print(data0_info)
	data0_info_tmp = data0[:,0:sig_col-1]
	data0_info = np.concatenate((data0_info, data0_info_tmp), axis = 0)
	print(data0_info)
	### extract sample signal matrix
	data0_signal = np.array(data0[:,sig_col-1:], dtype=float)
	print(data0_signal.shape)
	### read sample2celltype_list
	sample2celltype_list = read2d_list(sample2celltype, str)
	print(sample2celltype_list)
	### get cell type sample names (column header)
	header = []
	for i in range(0, len(sample2celltype_list)):
		header.append(sample2celltype_list[i][0])

	### get sample2matrix
	cellytpe_merge_matrix_info = sample2matrix(data0_signal, header, sample2celltype_list)
	print('cellytpe_merge_matrix_info.shape')
	print(cellytpe_merge_matrix_info.shape)
	print('data0_info.shape')
	print(data0_info.shape)	
	### add DNA region ID
	cellytpe_merge_matrix = np.concatenate((data0_info, cellytpe_merge_matrix_info), axis=1)
	
	### write table
	write2d_array(cellytpe_merge_matrix, outputfile)


############################################################################
#time python merge_cell_type_data_rsem.py -i homerTable3.peaks.filtered.binary_pattern.txt -n 6 -m sample2celltype.txt -o celltype.binary_pattern.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:n:m:o:")
	except getopt.GetoptError:
		print 'time python merge_cell_type_data.py -i inputfile -n sig_col -m sample2celltype -o outputfile'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python merge_cell_type_data.py -i inputfile -n sig_col -m sample2celltype -o outputfile'
			sys.exit()
		elif opt=="-i":
			inputfile=str(arg.strip())
		elif opt=="-n":
			sig_col=int(arg.strip())
		elif opt=="-m":
			sample2celltype=str(arg.strip())
		elif opt=="-o":
			outputfile=str(arg.strip())

	merge_cell_type_data(inputfile, sig_col, sample2celltype, outputfile)

if __name__=="__main__":
	main(sys.argv[1:])
