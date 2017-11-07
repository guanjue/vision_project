data = open('reads_count_matrix_5end_whole_genome.txt','r')
result = open('reads_count_matrix_5end_whole_genome.add_id.txt','w')

i=0
for records in data:
	tmp = records.split()
	result.write(tmp[0]+'\t'+tmp[1]+'\t'+tmp[2]+'\t'+str(i)+'\t'+str(int(tmp[2])-int(tmp[1]))+'\t' )
	for j in range(3,len(tmp)-1):
		result.write(tmp[j]+'\t')
	result.write(tmp[len(tmp)-1]+'\n')
	i = i+1
	if i%1000000 ==0:
		print(i)
result.close()
data.close()

