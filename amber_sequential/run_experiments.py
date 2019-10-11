import os
from subprocess import call
import shutil

runs = 10
hours = 24
secs = hours*3600
cores = 4

directories_temp = os.walk('.').next()[1]
# directories_temp = directories_temp[1:]
directories = []

for dir_temp in directories_temp:
	if(dir_temp  != 'histograms' and dir_temp != 'tools'):
		directories.append(dir_temp)

print directories
print "\nGen stride files... ",
for directory in directories:
	os.chdir(directory)
	# call(['pwd'])
	call(["stride",directory+'.pdb','-f'+directory+'.str'])
	os.chdir('..')
print "Done\n"

proteins = []
for directory in directories:
	os.chdir(directory)
	f = open(directory+'.str','r')
	amino_seq = ''
	amino_stride = ''
	amino_id = ''
	for line in f:
		if line[0:3] == 'SEQ':
			amino_seq += line.strip().split()[2]
			amino_id = line.strip().split()[4]
		if line[0:3] == 'ASG':
			amino_stride += line.strip().split()[5]
	print amino_id, len(amino_seq), amino_seq
	print amino_id, len(amino_stride), amino_stride,'\n'
	proteins.append([amino_id, amino_seq, amino_stride])
	f.close()
	os.chdir('..')

# f = open('set_proteins.txt','w')
for protein in proteins:
	for i in range(0,runs):
		newnamef = protein[0]+'-'+str(i+1)+'.nab'
		shutil.copy2('memetic.nab', protein[0]+'/'+newnamef)
		os.chdir(protein[0])
		fdata = None
		f = open(newnamef,'r')
		fdata = f.read()
		f.close()
		fdata = fdata.replace('#define RUN X','#define RUN\t\t\t\t\t'+str(i+1))
		fdata = fdata.replace('#define PROTEIN_ID "XXXX"','#define PROTEIN_ID\t\t\t"'+protein[0]+'"')
		fdata = fdata.replace('#define SEQ_SIZE XXXX','#define SEQ_SIZE\t\t\t'+str(len(protein[1])))
		fdata = fdata.replace('#define SEQ_STRING "XXXX"','#define SEQ_STRING\t\t\t"'+protein[1]+'"')
		fdata = fdata.replace('#define SEQ_SECSTRUCT "XXXX"','#define SEQ_SECSTRUCT\t\t"'+protein[2]+'"')
		fdata = fdata.replace('#define ARRAY_SIZE XXXX','#define ARRAY_SIZE\t\t\t'+str(len(protein[1])*6*78))
		f = open(newnamef,'w')
		f.write(fdata)
		f.close()
		call(['nab','-o',newnamef[:-4],newnamef,'../ran.c','../time.c'])
		# call(['./'+newnamef[:-4]])
		if os.path.exists('run_'+str(i+1)):
		    shutil.rmtree('run_'+str(i+1))
		os.makedirs('run_'+str(i+1))
		os.chdir('..')
	#f.write(protein[0] + " " + str(len(protein[1])) + " " + protein[1] + " " +protein[2]+ " \n")
# f.close()





