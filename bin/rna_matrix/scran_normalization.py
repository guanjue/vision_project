import os
import numpy as np

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used, sep):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split(sep)]
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
def scran_normalization(data_matrix, signal_col, scale_factor, outputname):
	###### read inputs
	### read data matrix
	data_matrix = read2d_array(data_matrix, str, '\t')

	### split header; info matrix; signal matrix
	header = [data_matrix[0,:]]
	data_matrix_info = data_matrix[1:, 0:signal_col-1]
	data_matrix_sig = np.array(data_matrix[1:, signal_col-1:], dtype=float)

	print('data_matrix_sig.shape:')
	print(data_matrix_sig.shape)
	######
	data_sf = read2d_array(scale_factor, float, '\t')

	###### normalize based on scran sf
	data_scran_matrix = data_matrix_sig * data_sf[0]

	### add back data info matrix
	data_scran_matrix = np.concatenate((data_matrix_info, data_scran_matrix),  axis = 1)
	### add back header
	data_scran_matrix = np.concatenate((header, data_scran_matrix),  axis = 0)
	
	###### write scran normed data
	write2d_array(data_scran_matrix, outputname)


############################################################################
############################################################################
#time python scran_normalization.py -i all_rna_sampe_target_class.txt -n 6 -s scran_scale_factor_file.txt -o all_rna_sampe_target_class_scran.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:n:s:o:")
	except getopt.GetoptError:
		print 'time python scran_normalization.py -i data_matrix -n signal_col -s scale_factor -o outputname'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python scran_normalization.py -i data_matrix -n signal_col -s scale_factor -o outputname'
			sys.exit()
		elif opt=="-i":
			data_matrix=str(arg.strip())
		elif opt=="-n":
			signal_col=int(arg.strip())	
		elif opt=="-s":
			scale_factor=str(arg.strip())	
		elif opt=="-o":
			outputname=str(arg.strip())		

	scran_normalization(data_matrix, signal_col, scale_factor, outputname)

if __name__=="__main__":
	main(sys.argv[1:])