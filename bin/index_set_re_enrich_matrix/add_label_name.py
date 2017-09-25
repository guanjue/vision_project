#enhancer_like_medium_signal=[8]
#enhancer_like_low_signal=[13]
#elongating_transcription=[9]
#enhancer_like_active_k27ac=[12]
#ctcf_and_or_enhancer_like=[14]
#promoter_like=[16]
#repressed=[0,1,3,4,5,6,7]
#atac_only=[15]
#low_signal=[2,11]
#quiescent=[10]
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

def add_label_name(inputfile, ideas_state_label, outputfile):
	### read ideas state lable table (label at the 2nd column)
	ideas_state_labels = read2d_array(ideas_state_label, str)[:,1]

	data0=open(inputfile,'r')
	result=open(outputfile,'w')
	### write label name to beside the state id
	for records in data0:
		tmp=[x.strip() for x in records.split('\t')]
		if tmp[3] !='NA':
			result.write(tmp[0]+'\t'+tmp[1]+'\t'+tmp[2]+'\t'+ideas_state_labels[int(tmp[3])]+';'+tmp[3]+'\t'+tmp[4]+'\t'+tmp[5]+'\t'+tmp[6]+'\t'+tmp[7]+'\t'+tmp[8]+'\n')
		else:
			result.write(tmp[0]+'\t'+tmp[1]+'\t'+tmp[2]+'\t'+tmp[3]+';'+tmp[3]+'\t'+tmp[4]+'\t'+tmp[5]+'\t'+tmp[6]+'\t'+tmp[7]+'\t'+tmp[8]+'\n')
	data0.close()
	result.close()


############################################################################
### python add_label_name.py -i used_DNA_intervals.neutrophil.8.colored.bed -l ideas_state_label.txt -o used_DNA_intervals.neutrophil.8.colored_named.bed
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"h:i:l:o:")
	except getopt.GetoptError:
		print 'python add_label_name.py -i used_DNA_intervals.neutrophil.8.colored.bed -l ideas_state_label -o used_DNA_intervals.neutrophil.8.colored_named.bed'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python add_label_name.py -i used_DNA_intervals.neutrophil.8.colored.bed -l ideas_state_label -o used_DNA_intervals.neutrophil.8.colored_named.bed'
			sys.exit()
		elif opt=="-i":
			inputfile=str(arg.strip())
		elif opt=="-l":
			ideas_state_label=str(arg.strip())
		elif opt=="-o":
			outputfile=str(arg.strip())

	add_label_name(inputfile, ideas_state_label, outputfile)
if __name__=="__main__":
	main(sys.argv[1:])
