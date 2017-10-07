def vlookup_bed(bed_file1,bed_file2,output_name, keep2nd):
	# read the first table into a list
	data01=open(bed_file1,'r')
	data011={}
	for records in data01:
		tmp=[x.strip() for x in records.split('\t')]
		bed_info = tmp[0]+'_'+tmp[1]+'_'+tmp[2]
		data011[bed_info]=tmp
	# read the second table into a dict
	data02=open(bed_file2,'r')
	data021=[]
	for records in data02:
		tmp2=[x.strip() for x in records.split('\t')]
		data021.append(tmp2)

	#print(data021)
	data03=[]
	for records in data021:
		# test if second table bed_info is in first table bed_info
		bed_info2=records[0]+'_'+records[1]+'_'+records[2]
		if bed_info2 in data011:
			match_data=[]
			# save the first table info
			for rec1 in data011[bed_info2]:
				match_data.append(rec1)
			if keep2nd == 'T':
				# save the second table info
				for rec2 in records:
					match_data.append(rec2)
			### append to final sorted matrix
			data03.append(match_data)

	data04=open(output_name,'w')
	for records in data03:
		for i in range(0,len(records)-1):
			data04.write(records[i]+'\t')
		data04.write(records[len(records)-1]+'\n')
	data04.close()
	data01.close()
	data02.close()


import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"ht:s:o:k:")
	except getopt.GetoptError:
		print 'python vlookup_bed.py -t first_table_file -s second_table_file -o output_name -k keep2nd_matrix'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python vlookup_bed.py -t first_table_file -s second_table_file -o output_name -k keep2nd_matrix'
			sys.exit()
		elif opt=="-t":
			bed_file1=str(arg.strip())
		elif opt=="-s":
			bed_file2=str(arg.strip())
		elif opt=="-o":
			output_name=str(arg.strip())
		elif opt=="-k":
			keep2nd=str(arg.strip())
			
	vlookup_bed(bed_file1,bed_file2,output_name, keep2nd)

if __name__=="__main__":
	main(sys.argv[1:])