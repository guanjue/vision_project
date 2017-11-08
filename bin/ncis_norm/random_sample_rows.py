import numpy as np

def random_sample_rows(inputfile, sample_number, random_seed, use_header, outputfile):
	data = open(inputfile,'r')
	result = open(outputfile,'w')

	### write header
	if use_header == 'T':
		result.write(data.readline())

	### read matrix
	data0 = []
	### read reads count matrix
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		data0.append(tmp)
	data0 = np.array(data0)
	data.close()

	### random sampling
	np.random.seed(random_seed)
	idx = np.random.randint(data0.shape[0], size=sample_number)

	### sampling matrix
	data0_sample = data0[idx, :]

	### write output
	for records in data0_sample:
		for i in range(0,len(records)-1):
			result.write(records[i]+'\t')
		result.write(records[len(records)-1]+'\n')
	result.close()

###########################################
# time python random_sample_rows.py -i homerTable3.peaks.filtered.reads.txt -n 1000000 -s 2017 -d T -o whole_genome_sampling_regions.reads.txt
import getopt
import sys

def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:n:s:d:o:")
	except getopt.GetoptError:
		print 'python random_sample_rows.py -i <input matrix> -n <random sampling number> -s <random seed> -d <use header> -o <output file>'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python random_sample_rows.py -i <input matrix> -n <random sampling number> -s <random seed> -d <use header> -o <output file>'
			sys.exit()
		elif opt=="-i":
			inputfile=str(arg.strip())
		elif opt=="-n":
			sample_number=int(arg.strip())
		elif opt=="-s":
			random_seed=int(arg.strip())
		elif opt=="-d":
			use_header=str(arg.strip())
		elif opt=="-o":
			outputfile=str(arg.strip())
	random_sample_rows(inputfile, sample_number, random_seed, use_header, outputfile)

if __name__=="__main__":
	main(sys.argv[1:])

