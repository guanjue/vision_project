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


### read given cell type order
colmun_oder = read2d_array('celltype.order.txt',str)
order = np.array(np.transpose(colmun_oder[:,0]),dtype=int)
header = np.transpose(colmun_oder[:,1])
print('cell type order:')
print(header)

### read given sorting signal level
signal_level_range = read2d_array('signal_level_range.txt',float)
print('signal level range:')
print(signal_level_range)

### read primary matrix
prime_mark0 = read2d_array('celltype.binary_pattern.txt',str)
prime_mark_id = prime_mark0[:,0]
### sort cell types by given cell type order
prime_mark_info_sorted = np.array(prime_mark0[1:,order], dtype = float)

### change signal to user defined signal levels
for levels in signal_level_range:
	prime_mark_info_sorted[((prime_mark_info_sorted>=levels[1]) * (prime_mark_info_sorted<levels[2]))] = levels[0]

### sort rows based on index pattern
for i in range(0,16)[::-1]:
	prime_mark_info_sorted = prime_mark_info_sorted[-prime_mark_info_sorted[:,i].argsort(kind='mergesort'),:]

write2d_array(prime_mark_info_sorted, 'celltype.binary_pattern.sorted.txt')

index_set_dict = {}
index_set = []
for index in prime_mark_info_sorted:
	if not (index in index_set_dict):
		index_set_dict[index] = 1
		index_set.append(index)
	else:
		index_set_dict[index] = index_set_dict[index]+1

result = open('celltype.index_set.sorted.txt','w')

for pattern in index_set:
	pattern_num = index_set_dict[pattern]
	for i in range(0, len(pattern)-1):
		result.write(pattern[i] * pattern_num)

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
