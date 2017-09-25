import numpy as np

def reads_count_rpkm(input_readscount_matrix, total_reads_count_table, output_rpkm_matrix):
	data = open(input_readscount_matrix,'r')
	data0 = []
	### read reads count matrix
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		data0.append(tmp)
	data0 = np.array(data0)
	data.close()

	### read total reads
	data_t = open(total_reads_count_table,'r')
	data0_t = []
	### read matrix
	for records in data_t:
		tmp = [x.strip() for x in records.split('\t')]
		data0_t.append(tmp)
	data0_t = np.array(data0_t, dtype = float)
	data_t.close()

	### get rpkm & write output
	result = open(output_rpkm_matrix, 'w')

	for records in data0:
		### write chr start end name
		result.write(records[0]+'\t'+records[1]+'\t'+records[2]+'\t'+records[3]+'\t')
		### write rpkm
		for i in range(5,len(records)-1):
			C = float(records[i])
			N = float(data0_t[i-5])
			L = float(records[4])
			### calculate rpkm
			rpkm = round((10**9 * C)/(N * L), 2)
			result.write(str(rpkm)+'\t')
		C = float(records[len(records)-1])
		N = float(data0_t[len(records)-1-5])
		L = float(records[4])
		rpkm = round((10**9 * C)/(N * L), 2)

		result.write(str(rpkm)+'\n')
	result.close()


############################################################################
#time python reads_count_rpkm.py -i reads_count_matrix_5end.txt -t total_reads_all.txt -o reads_count_matrix_5end_rpkm.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:t:o:")
	except getopt.GetoptError:
		print 'time python reads_count_rpkm.py -i input_readscount_matrix -t total_reads_count_table -o output_rpkm_matrix'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python reads_count_rpkm.py -i input_readscount_matrix -t total_reads_count_table -o output_rpkm_matrix'
			sys.exit()
		elif opt=="-i":
			input_readscount_matrix=str(arg.strip())
		elif opt=="-t":
			total_reads_count_table=str(arg.strip())
		elif opt=="-o":
			output_rpkm_matrix=str(arg.strip())

	reads_count_rpkm(input_readscount_matrix, total_reads_count_table, output_rpkm_matrix)

if __name__=="__main__":
	main(sys.argv[1:])