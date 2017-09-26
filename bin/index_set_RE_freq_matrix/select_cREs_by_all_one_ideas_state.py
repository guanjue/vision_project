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
def select_cREs_by_all_one_ideas_state(target_matrix, active_state_list, output_cRE, skip_state):
	###### read inputs
	### read index_matrix
	index_matrix = read2d_array(target_matrix, str)
	### read active_state
	active_state = int(read2d_array(active_state_list, str)[0,0])
	print(active_state)

	###
	cRE_active_mark = []
	cRE_active_mark.append(index_matrix[0,:])
	for records in index_matrix[1:,:]:
		check_num = 0 ### initialize non-active state num 
		for s in records[1:]: ### check each state in 16 cell types
			if int(s) != skip_state: ### skip info cell type info
				if int(s) != active_state: ### check if non-active state appears
					check_num = check_num +1
		if check_num == 0+1: ### check_num == 0 means all states are in active state list
			cRE_active_mark.append(records)
	### write output
	write2d_array(cRE_active_mark, output_cRE)
############################################################################
############################################################################
#time python select_cREs_by_all_one_ideas_state.py -t celltype.index.sorted.txt -a ideas_active_state_label.txt -o celltype.index.sorted.active_state_filtered.txt -s 17

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"ht:a:o:s:")
	except getopt.GetoptError:
		print 'time python select_cREs_by_all_one_ideas_state.py -t target_matrix -a active_state_list -o output_cRE -s skip_state'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python select_cREs_by_all_one_ideas_state.py -t target_matrix -a active_state_list -o output_cRE -s skip_state'
			sys.exit()
		elif opt=="-t":
			target_matrix=str(arg.strip())
		elif opt=="-a":
			active_state_list=str(arg.strip())
		elif opt=="-o":
			output_cRE=str(arg.strip())		
		elif opt=="-s":
			skip_state=str(arg.strip())		

	select_cREs_by_all_one_ideas_state(target_matrix, active_state_list, output_cRE, skip_state)

if __name__=="__main__":
	main(sys.argv[1:])