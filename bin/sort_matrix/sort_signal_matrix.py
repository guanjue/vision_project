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
### write index sets 
def get_index_set(inputmatrix, given_column_order, given_signal_level, sorted_index_output, sorted_index_set_output):
	### read given cell type order
	column_oder = read2d_array(given_column_order,str)
	order = np.array(np.transpose(column_oder[:,0]),dtype=int)
	header = np.array([column_oder[:,1]])
	print('cell type order:')
	print(header)
