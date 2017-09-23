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
def split_signal_binary_matrix(inputfile, output1, output1_col, output2, output2_col, output3, output3_col):
	data_all = read2d_array(inputfile,str)

	### split signal matrix & binary matrix
	interval_info = data_all[:,0:output1_col]
	signal = data_all[:,output1_col:output2_col]
	binary = data_all[:,output2_col:output3_col]
	### write results
	write2d_array(interval_info, output1)
	write2d_array(signal, output2)
	write2d_array(binary, output3)


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

	split_signal_binary_matrix(inputfile, output1, output1_col, output2, output2_col, output3, output3_col)

if __name__=="__main__":
	main(sys.argv[1:])
