import numpy as np

def get_ncis_t_a_b(matrix1, matrix2, output):
	data1 = open(matrix1,'r')
	data2 = open(matrix2,'r')

	data1.readline()
	data2.readline()

	data1_line = data1.readlines()
	data2_line = data2.readlines()

	for w in range(1,2):
		d1_tri1 = {}
		d1_tri2 = {}
		c = 0
		r1_vec = []
		r2_vec = []
		t_vec = []

		for i in range(1,10000):
			d1_tri1[int(i)] =[ 0.0 ]
			d1_tri2[int(i)] =[ 0.0 ]

		l_0 = 13489996#len(data1)
		l = l_0 - l_0%w
		index_all = range(0, l-w, w)
		#index_sample = np.random.choice(index_all, 13000000)
		for i in index_all: ### loop matrix at x axis
			if c%100000 == 0:
				print(c)
			r1_vec = []
			r2_vec = []

			for j in range(0,w):
				d1 = [ x.strip() for x in data1_line[i+j].split('\t') ]
				d2 = [ x.strip() for x in data2_line[i+j].split('\t') ]
				r1 = float(d1[1])*2
				r2 = float(d2[1])*2
				r1_vec.append(r1)
				r2_vec.append(r2)
				t = r1+r2
				#t_vec.append(t)
				c = c+1
			r1_vec_sum = np.sum(r1_vec)
			r2_vec_sum = np.sum(r2_vec)
			t_vec_sum = r1_vec_sum+r2_vec_sum
			if t_vec_sum in d1_tri1:
				d1_tri1[int(t_vec_sum)].append( r1_vec_sum )
				d1_tri2[int(t_vec_sum)].append( r2_vec_sum )
			else:
				d1_tri1[int(t_vec_sum)] = [ r1_vec_sum ]
				d1_tri2[int(t_vec_sum)] = [ r2_vec_sum ]
		data1.close()
		data2.close()

		print('read table DONE')
		#d1_tri_sum = d1_tri1.keys() ### get t for NCIS
	
		#t_list = np.unique(d1_tri_sum) ### get unique t 
		#print(t_list)
		#t_order = np.argsort(t_list) ### sort t id
		#t_list_sort = t_list[t_order] ### sort t

		d1_tri_tsum =[]
		for t in range(1,10000):
			d1_tmp = np.sum(d1_tri1[t]) ### get matrix1 sum, when matrix1+matrix2 = ti
			d2_tmp = np.sum(d1_tri2[t]) ### get matrix1 sum, when matrix1+matrix2 = ti
			d1_tri_tsum.append( [ t, d1_tmp, d2_tmp ] )
		d1_tri_tsum = np.array(d1_tri_tsum)
		print('ncis table DONE')
		#print(d1_tri_tsum)
		result = open(output+'.'+str(w*200)+'.txt','w')
		for records in d1_tri_tsum:
			result.write(str(records[0])+'\t'+str(records[1])+'\t'+str(records[2])+'\n')
		result.close()


###########################################
# time python get_atac_t_ncis.py -i LSK_BM.h3k4me3rep.100056.tab -j input_tab/LSK_BM.h3k4me3rep.100056.input.tab -o LSK_BM.h3k4me3rep.100056.ncsi.matrix.txt
import getopt
import sys

def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:j:o:")
	except getopt.GetoptError:
		print 'python extract_statepair_position.py -i <idea state at each position for two conditions> -t <target cell data col num> -t <source cell data col num> -o <output file>'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python extract_statepair_position.py -i <idea state at each position for two conditions> -t <target cell data col num> -t <source cell data col num> -o <output file>'
			sys.exit()
		elif opt=="-i":
			matrix1=str(arg.strip())
		elif opt=="-j":
			matrix2=str(arg.strip())
		elif opt=="-o":
			output=str(arg.strip())
	get_ncis_t_a_b(matrix1, matrix2, output)

if __name__=="__main__":
	main(sys.argv[1:])



