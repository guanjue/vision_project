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
### remove columns
def remove_columns(input_matrix, remove_col_list, output_matrix):
	### read reads count matrix
	data0 = read2d_array(input_matrix, str)

	### get remove col list
	#print(remove_col_list.split(','))
	rm_col = [int(x.strip())-1 for x in remove_col_list.split(',')]

	### remove columns
	data0 = np.delete(data0, rm_col, axis=1)

	### write output
	write2d_array(data0, output_matrix)

############################################################################
#time python remove_columns.py -i input_matrix -t rm_cols:5,6,7 -o output_matrix

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:r:o:")
	except getopt.GetoptError:
		print 'time python remove_columns.py -i input_matrix -t rm_cols:5,6,7 -o output_matrix'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python remove_columns.py -i input_matrix -t rm_cols:5,6,7 -o output_matrix'
			sys.exit()
		elif opt=="-i":
			input_matrix=str(arg.strip())
		elif opt=="-r":
			remove_col_list=str(arg.strip())
		elif opt=="-o":
			output_matrix=str(arg.strip())

	remove_columns(input_matrix, remove_col_list, output_matrix)

if __name__=="__main__":
	main(sys.argv[1:])
