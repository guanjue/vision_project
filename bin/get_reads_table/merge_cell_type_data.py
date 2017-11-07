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
### sample2matrix
def sample2matrix(signal_matrix, col_header, sample2celltype_matrix, sample_size):
	### initiate header vector & signal matrix
	cellytpe_merge_matrix_header = []
	cellytpe_merge_matrix_data = np.array([], dtype=np.float).reshape(0, signal_matrix.shape[0])

	### merge sample based on sample2celltype_matrix
	for i in range(0, sample2celltype_matrix.shape[0]):
		### column header
		tmp_header = col_header[int(sample2celltype_matrix[i, 1])].split('_')[0]

		### average column signal
		tmp_sig = signal_matrix[:,int(sample2celltype_matrix[i, 1])]
		for j in range(1,sample_size):
			tmp_sig = tmp_sig + signal_matrix[:,int(sample2celltype_matrix[i, 1+j])]
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
def merge_cell_type_data(inputfile, sample2celltype, sample_size, outputfile):
	data0 = read2d_array(inputfile,str)
	### get cell type sample names (column header)
	header = data0[0,1:]
	data0_id = np.transpose([data0[:,0]])
	data0_signal = np.array(data0[1:,1:], dtype=float)

	### read sample2celltype_matrix
	sample2celltype_matrix = read2d_array(sample2celltype,str)

	### get sample2matrix
	cellytpe_merge_matrix_info = sample2matrix(data0_signal, header, sample2celltype_matrix, sample_size)
	
	### add DNA region ID
	cellytpe_merge_matrix = np.concatenate((data0_id, cellytpe_merge_matrix_info), axis=1)
	
	### write table
	write2d_array(cellytpe_merge_matrix, outputfile)


############################################################################
#time python merge_cell_type_data.py -i homerTable3.peaks.filtered.binary_pattern.txt -m sample2celltype.txt -n 2 -o celltype.binary_pattern.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:m:n:o:")
	except getopt.GetoptError:
		print 'time python merge_cell_type_data.py -i inputfile -m sample2celltype -n sample_size -o outputfile'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python merge_cell_type_data.py -i inputfile -m sample2celltype -n sample_size -o outputfile'
			sys.exit()
		elif opt=="-i":
			inputfile=str(arg.strip())
		elif opt=="-m":
			sample2celltype=str(arg.strip())
		elif opt=="-n":
			sample_size=int(arg.strip())
		elif opt=="-o":
			outputfile=str(arg.strip())

	merge_cell_type_data(inputfile, sample2celltype, sample_size, outputfile)

if __name__=="__main__":
	main(sys.argv[1:])
