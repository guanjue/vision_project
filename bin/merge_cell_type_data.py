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
def merge_cell_type_data(inputfile, outputfile):
	data0 = read2d_array(inputfile,str)
	### get cell type sample names (column header)
	header = data0[0,1:]
	data0_id = np.transpose([data0[:,0]])
	data0_signal = np.array(data0[1:,1:], dtype=float)

	### get cell header split by '_' (remove sample id)
	lsk_header = [header[0].split('_')[0]]
	cmp_header = [header[3].split('_')[0]]
	gmp_header = [header[5].split('_')[0]]
	mep_header = [header[7].split('_')[0]]
	cfue_header = [header[9].split('_')[0]]
	ery_header = [header[11].split('_')[0]]
	cfumk_header = [header[13].split('_')[0]]
	meg_header = [header[15].split('_')[0]]
	mono_header = [header[17].split('_')[0]]
	neu_header = [header[19].split('_')[0]]
	b_header = [header[20].split('_')[0]]
	nk_header = [header[21].split('_')[0]]
	tcd4_header = [header[22].split('_')[0]]
	tcd8_header = [header[23].split('_')[0]]
	g1e_header = [header[24].split('_')[0]]
	er4_header = [header[26].split('_')[0]]

	cellytpe_merge_matrix_header = np.array( ([lsk_header + cmp_header + gmp_header + mep_header + cfue_header + ery_header + cfumk_header + meg_header + mono_header + neu_header + b_header + nk_header + tcd4_header + tcd8_header + g1e_header + er4_header]) )
	
	### get cell average
	lsk = np.array([(data0_signal[:,0] + data0_signal[:,1])/2])
	cmp = np.array([(data0_signal[:,2] + data0_signal[:,3])/2])
	gmp = np.array([(data0_signal[:,4] + data0_signal[:,5])/2])
	mep = np.array([(data0_signal[:,6] + data0_signal[:,7])/2])
	cfue = np.array([(data0_signal[:,8] + data0_signal[:,9])/2])
	ery = np.array([(data0_signal[:,10] + data0_signal[:,11])/2])
	cfumk = np.array([(data0_signal[:,12] + data0_signal[:,13])/2])
	meg = np.array([(data0_signal[:,14] + data0_signal[:,15])/2])
	mono = np.array([(data0_signal[:,16] + data0_signal[:,17])/2])
	neu = np.array([(data0_signal[:,18] + data0_signal[:,19])/2])
	b = np.array([data0_signal[:,20]])
	nk = np.array([data0_signal[:,21]])
	tcd4 = np.array([data0_signal[:,22]])
	tcd8 = np.array([data0_signal[:,23]])
	g1e = np.array([(data0_signal[:,24] + data0_signal[:,25])/2])
	er4 = np.array([(data0_signal[:,26] + data0_signal[:,27])/2])

	### merge columns
	cellytpe_merge_matrix_data = np.transpose(np.concatenate((lsk, cmp, gmp, mep, cfue, ery, cfumk, meg, mono, neu, b, nk, tcd4, tcd8, g1e, er4), axis=0))
	
	### add column header
	cellytpe_merge_matrix_info = np.concatenate((cellytpe_merge_matrix_header, cellytpe_merge_matrix_data), axis=0)
	
	### add DNA region ID
	cellytpe_merge_matrix = np.concatenate((data0_id, cellytpe_merge_matrix_info), axis=1)
	
	### write table
	write2d_array(cellytpe_merge_matrix, outputfile)


############################################################################
#time python merge_cell_type_data.py -i homerTable3.peaks.filtered.binary_pattern.txt -o celltype.binary_pattern.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:o:")
	except getopt.GetoptError:
		print 'time python merge_cell_type_data.py -i inputfile -a output1 -o outputfile'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python merge_cell_type_data.py -i inputfile -a output1 -o outputfile'
			sys.exit()
		elif opt=="-i":
			inputfile=str(arg.strip())
		elif opt=="-o":
			outputfile=str(arg.strip())

	merge_cell_type_data(inputfile, outputfile)

if __name__=="__main__":
	main(sys.argv[1:])
