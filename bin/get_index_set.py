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


def get_index_set(inputmatrix, given_column_order, given_signalrange, sorted_index_output, sorted_index_set_output):
	### read given cell type order
	column_oder = read2d_array('celltype.order.txt',str)
	order = np.array(np.transpose(column_oder[:,0]),dtype=int)
	header = np.array([column_oder[:,1]])
	print('cell type order:')
	print(header)

	### read given sorting signal level
	signal_level_range = read2d_array('signal_level_range.txt',float)
	print('signal level range:')
	print(signal_level_range)

	### read primary matrix
	prime_mark0 = read2d_array('celltype.binary_pattern.txt',str)
	prime_mark_id = np.transpose([prime_mark0[1:,0]])
	### sort cell types by given cell type order
	prime_mark_info_sorted = np.array(prime_mark0[1:,order], dtype = float)

	### change signal to user defined signal levels
	for levels in signal_level_range:
		prime_mark_info_sorted[((prime_mark_info_sorted>=levels[1]) * (prime_mark_info_sorted<levels[2]))] = levels[0]


	print(prime_mark_id.shape)
	### sort rows based on index pattern
	for i in range(0,16)[::-1]:
		row_order = prime_mark_info_sorted[:,i].argsort(kind='mergesort')
		prime_mark_info_sorted = prime_mark_info_sorted[row_order,:]
		prime_mark_id = prime_mark_id[row_order,:]
	prime_mark_id_sorted = np.concatenate((np.array([['DNA_region_id']]), prime_mark_id),axis = 0)
	print(header.shape)
	print(prime_mark_info_sorted.shape)
	print(prime_mark_id.shape)
	### add header
	prime_mark_info_sorted_info = np.concatenate((header, prime_mark_info_sorted), axis = 0)
	### add row names
	prime_mark_info_sorted_matrix = np.concatenate((prime_mark_id_sorted, prime_mark_info_sorted_info), axis = 1)
	### write the all DNA region binary pattern
	write2d_array(prime_mark_info_sorted_matrix, 'celltype.binary_pattern.sorted.txt')

	prime_mark_info_sorted_matrix_reliable = [prime_mark_info_sorted_matrix[0]]
	for i in range(1, prime_mark_info_sorted_matrix.shape[0]):
		pattern = np.array(prime_mark_info_sorted_matrix[i,1:], dtype = float)
		info = np.sum(pattern)
		if info > 0:
			prime_mark_info_sorted_matrix_reliable.append(prime_mark_info_sorted_matrix[i,:])

	write2d_array(prime_mark_info_sorted_matrix_reliable, 'celltype.binary_pattern.sorted.filtered.txt')


	index_set_dict = {}
	index_set = []
	for index in prime_mark_info_sorted:
		index_merge = ''
		for i in index:
			index_merge = index_merge+'_'+str(i)
		if not (index_merge in index_set_dict):
			index_set_dict[index_merge] = 1
			index_set.append(index_merge)
		else:
			index_set_dict[index_merge] = index_set_dict[index_merge]+1


	result = open('celltype.index_set.sorted.txt','w')
	result.write('index'+'\t')
	for i in range(0,header.shape[1]-1):
		result.write(header[0,i]+'\t')
	result.write(str(header[0,header.shape[1]-1])+'\n')

	for pattern in index_set:
		pattern_num = index_set_dict[pattern]
		binary_label = pattern.split('_')
		if pattern_num >= 200:
			result.write(pattern[1:]+'\t')
			for i in range(1, len(binary_label)-1):
				result.write(str(float(binary_label[i]) * pattern_num)+'\t')
			result.write(str(float(binary_label[len(binary_label)-1]) * pattern_num)+'\n')
	result.close()


############################################################################
#time python split_signal_binary_matrix.py -i homerTable3.peaks.filtered.txt -a homerTable3.peaks.filtered.interval.txt -x 4 -b homerTable3.peaks.filtered.signal.txt -y 32 -c homerTable3.peaks.filtered.binary_pattern.txt -z 60

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:a:b:c:x:y:z:")
	except getopt.GetoptError:
		print 'time python split_signal_binary_matrix.py -i inputfile a output1 -x output1_col -b output2 -y output2_col -c output3 -z output3_col'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python split_signal_binary_matrix.py -i inputfile a output1 -x output1_col -b output2 -y output2_col -c output3 -z output3_col'
			sys.exit()
		elif opt=="-i":
			inputfile=str(arg.strip())
		elif opt=="-a":
			output1=str(arg.strip())
		elif opt=="-b":
			output2=str(arg.strip())		
		elif opt=="-c":
			output3=str(arg.strip())
		elif opt=="-x":
			output1_col=int(arg.strip())
		elif opt=="-y":
			output2_col=int(arg.strip())		
		elif opt=="-z":
			output3_col=int(arg.strip())

	get_index_set(inputfile, output1, output1_col, output2, output2_col, output3, output3_col)

if __name__=="__main__":
	main(sys.argv[1:])
