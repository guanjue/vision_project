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
### read 2d array to dict
def read2d_array_dict(filename, colsort, id_col):
	data=open(filename,'r')
	data0={}
	for records in data:
		tmp = np.array([x.strip() for x in records.split('\t')])
		### col sort
		tmp = tmp[colsort]
		data0[tmp[id_col]]=tmp
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
def get_index_signal_matrix(target_matrix, t_id_col, source_matrix, s_id_col, given_column_order, target_matrix_sorted):
	### read given cell type order
	column_oder = read2d_array(given_column_order,str)
	order = np.array(np.transpose(column_oder[:,0]),dtype=int)
	header = np.array([column_oder[:,1]])
	print('cell type order:')
	print(header)

	### read target matrix
	target_mark0 = read2d_array_dict(inputmatrix, order, t_id_col-1)

	### read source matrix
	source_mark0 = read2d_array(inputmatrix,str)

	target_mark_sorted = []
	for records in source_mark0:
		target_mark_sorted.append( target_mark0[records[s_id_col-1]] )


	write2d_array(target_mark_sorted, target_matrix_sorted)




############################################################################
#time python get_index_signal_matrix.py -t celltype.signal.txt -a 1 -s celltype.index.sorted.txt -b 1 -r celltype.order.txt -o celltype.index.signal.sorted.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"ht:a:s:b:r:o:")
	except getopt.GetoptError:
		print 'time python get_index_set.py -i inputmatrix -r celltype_order -l signal_level_range -f index_output -s index_set_output'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python get_index_set.py -i inputmatrix -r celltype_order -l signal_level_range -f index_output -s index_set_output'
			sys.exit()
		elif opt=="-t":
			target_matrix=str(arg.strip())
		elif opt=="-a":
			t_id_col=int(arg.strip())
		elif opt=="-s":
			source_matrix=str(arg.strip())		
		elif opt=="-b":
			s_id_col=int(arg.strip())
		elif opt=="-r":
			given_column_order=str(arg.strip())
		elif opt=="-o":
			target_matrix_sorted=str(arg.strip())

	get_index_signal_matrix(target_matrix, t_id_col, source_matrix, s_id_col, given_column_order, target_matrix_sorted)

if __name__=="__main__":
	main(sys.argv[1:])
