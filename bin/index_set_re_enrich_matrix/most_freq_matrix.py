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
	data.readline() # skip header
	data0={}
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		### remove labels
		for i in range(1, len(tmp)):
			tmp[i] = tmp[i].split(';')[1]
		### split id & data
		id_info = tmp[0]
		data_info = np.array(tmp)
		### col sort
		data_info = data_info[colsort]
		data0[id_info] = data_info
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
def get_index_signal_matrix(target_matrix, t_id_col, source_matrix, s_id_col, given_column_order, quantile, target_matrix_sorted_output, sorted_index_set_signal_output):
	### read given cell type order
	column_oder = read2d_array(given_column_order,str)
	order = np.array(np.transpose(column_oder[:,0]),dtype=int)
	header = np.array([column_oder[:,1]])
	print('cell type order:')
	print(header)

	### read target matrix
	target_mark0 = read2d_array_dict(target_matrix, order, t_id_col-1)
	### read source matrix
	source_mark0 = read2d_array(source_matrix,str)

	### sort target matrix
	target_mark_sorted = []
	for records in source_mark0:
		id_label = records[s_id_col-1]
		### re-array the target matrix by the source matrix id order
		if id_label in target_mark0: ### re-array the target matrix by the source matrix id order
			target_mark_sorted.append( target_mark0[id_label] )

	### add header
	source_mark_id = np.transpose([source_mark0[1:,0]])
	source_mark_id_sorted = np.concatenate((np.array([['name']]), source_mark_id),axis = 0)
	target_mark_sorted_info = np.concatenate((header, target_mark_sorted), axis = 0)
	### add row names
	target_mark_info_sorted_matrix = np.concatenate((source_mark_id_sorted, target_mark_sorted_info), axis = 1)

	### write sorted signal matrix
	write2d_array(target_mark_info_sorted_matrix, target_matrix_sorted_output)


	### get index sets median signal
	index_set_sig_dict = {}
	index_set = []
	for index_info, signal_info in zip( source_mark0[1:,:], target_mark_sorted):
		index_merge = ''
		for i in index_info[1:]:
			index_merge = index_merge+'_'+str(i)
		### get the number of DNA regions in index set
		if 'NA' in signal_info:
			print(signal_info)
		if not (index_merge in index_set_sig_dict):
			index_set_sig_dict[index_merge] = [ signal_info ]
			index_set.append(index_merge)
		else:
			index_set_sig_dict[index_merge].append(signal_info)

	### write index set median signal matrix
	def write_index_set_signal_matrix(output_name, index_set_sig_dict, index_set_array, header, quantile):
		result = open(output_name,'w')
		### write header
		result.write('name'+'\t')
		for i in range(0,header.shape[1]-1):
			result.write(header[0,i]+'\t')
		result.write(str(header[0,header.shape[1]-1])+'\n')

		### write index set signal info
		for pattern in index_set_array:
			### get index set signal median for each cell type
			index_set_ideas_matrix = np.array(index_set_sig_dict[pattern], dtype = int)
			signal_vector = []
			for i in range(0,index_set_ideas_matrix.shape[1]):
				index_set_ideas_matrix_celltype_tmp = index_set_ideas_matrix[:,i]
				most_frequent_state = np.argmax(np.bincount(index_set_ideas_matrix_celltype_tmp))
				signal_vector.append(most_frequent_state)

			### signal vector info
			result.write(pattern[1:]+'\t')
			for i in range(0, len(signal_vector)-1):
				result.write(str(signal_vector[i])+'\t') 
			result.write(str(signal_vector[len(signal_vector)-1])+'\n')
		result.close()

	write_index_set_signal_matrix(sorted_index_set_signal_output, index_set_sig_dict, index_set, header, quantile)

############################################################################
############################################################################
#time python get_index_signal_matrix.py -t celltype.signal.txt -a 1 -s celltype.index.sorted.txt -b 1 -r celltype.order.txt -q 75 -o celltype.index.signal.sorted.txt -p celltype.index_set.signal.sorted.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"ht:a:s:b:r:q:o:p:")
	except getopt.GetoptError:
		print 'time python get_index_signal_matrix.py -t target_matrix -a target_id_col -s source_matrix -b source_id_col -r celltype.order -q quantile -o index.signal.output -p index_set.signal.output'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python get_index_signal_matrix.py -t target_matrix -a target_id_col -s source_matrix -b source_id_col -r celltype.order -q quantile -o index.signal.output -p index_set.signal.output'
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
		elif opt=="-q":
			quantile=float(arg.strip())
		elif opt=="-o":
			target_matrix_sorted=str(arg.strip())
		elif opt=="-p":
			sorted_index_set_signal_output=str(arg.strip())

	get_index_signal_matrix(target_matrix, t_id_col, source_matrix, s_id_col, given_column_order, quantile, target_matrix_sorted, sorted_index_set_signal_output)

if __name__=="__main__":
	main(sys.argv[1:])