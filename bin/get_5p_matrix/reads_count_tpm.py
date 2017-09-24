import numpy as np

def reads_count_tpm(input_readscount_matrix, output_tpm_matrix):
	data = open(input_readscount_matrix,'r')
	data0 = []
	### read reads count matrix
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		data0.append(tmp)
	data0 = np.array(data0)
	data.close()

	### split the row info & row data
	data0_sig = np.array(data0[:,5-1:], dtype = float)

	### get reads count per bp 
	r_l = []
	for sigs in data0_sig:
		l = sigs[0]
		r_l_tmp = sigs[1:] / l
		r_l.append(r_l_tmp)
	r_l = np.array(r_l)

	### get sum of (reads count per bp) 
	r_l_sum = np.sum(r_l, axis = 0)

	### get tpm & write output
	result = open(output_tpm_matrix, 'w')
	for records, r_l_sig in zip(data0, r_l):
		### write row infos (chr start end id)
		result.write(records[0]+'\t'+records[1]+'\t'+records[2]+'\t'+records[3]+'\t')
		for i in range(0,len(r_l_sig)-1):
			### calculate tpm
			tpm = r_l_sig[i] / r_l_sum[i] * 10**6
			result.write(str(tpm)+'\t')
		tpm = r_l_sig[len(r_l_sig)-1] / r_l_sum[len(r_l_sig)-1] * 10**6
		result.write(str(tpm)+'\n')

	result.close()

############################################################################
#time python reads_count_tpm.py -i reads_count_matrix_5end.txt -o reads_count_matrix_5end_tpm.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:o:")
	except getopt.GetoptError:
		print 'time python reads_count_tpm.py -i input_readscount_matrix -o output_tpm_matrix'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python reads_count_tpm.py -i input_readscount_matrix -o output_tpm_matrix'
			sys.exit()
		elif opt=="-i":
			input_readscount_matrix=str(arg.strip())
		elif opt=="-o":
			output_tpm_matrix=str(arg.strip())

	reads_count_tpm(input_readscount_matrix, output_tpm_matrix)

if __name__=="__main__":
	main(sys.argv[1:])