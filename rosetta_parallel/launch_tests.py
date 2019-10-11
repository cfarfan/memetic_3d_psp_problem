import os
import subprocess
import shutil
import urllib2
import logging

#protein_ids = ['1AB1','1AIL','1ENH','1L2Y','1ROP','2F4K','2MQ8','2MTW','2P5K','2P81','3P7K','1ACW','1DFN','1K43','1Q2K','1WQC','2JUC','2MR9','2MW1','2P6J','2PMR','3V1A']
#protein_ids = ['1ACW','1DFN','1K43','1L2Y','1Q2K','1WQC','2MTW','2P81']
#protein_ids = ['1ACW']
#protein_ids = ['1Q2K','1WQC','2MTW','2P81']
protein_ids = ['1ACW','1DFN','1K43','1L2Y']

test_per_protein = 1
tests_concurrent = 2

config = ['13 10 3 True','13 10 3 False','10 5 2 True','10 5 2 False']
#config = ['13 10 3 True','13 10 3 False']

def download_pdbs(pdb_list):
	rcsb_restful_url = 	'http://www.rcsb.org/pdb/files/'
	download_path 	 =	'protein_pdbs/'
	
	if not os.path.exists(download_path):
		os.makedirs(download_path)

	logging.info('\033[32mDownloading proteins: %s\033[39m' % (pdb_list) )

	for protein_id in pdb_list:
		if(os.path.isfile(os.path.join(download_path,protein_id+'.pdb'))):
			logging.info('\033[33mProtein %s already exists in %s\033[39m' % (protein_id, download_path) )
		else:
			req = rcsb_restful_url+protein_id+'.pdb'
			response = urllib2.urlopen(req)
			pdb = response.read()
			f = open(os.path.join(download_path,protein_id+'.pdb'),'w')
			f.write(pdb)
			f.close()
			logging.info( '\033[32mProtein %s has been downloaded from %s to %s\033[39m' % (protein_id, req, os.path.abspath(download_path)))

def run():

	tests = []
	for protein in protein_ids:
		for conf in config:
			for i in range(test_per_protein):
				tests.append(protein+' '+conf)

	test_run = 0
	process = []
	while not (len(tests) == 0):
		if test_run < tests_concurrent:
			test_protein = tests.pop()
			test_run += 1
			cmd = """python memetic_paralell.py %s >/dev/null 2>&1""" % (test_protein)
			logging.info('\033[32m'+cmd+'\033[39m')
			process.append(subprocess.Popen(cmd, shell=True))
		else:
			for p in process:
				p.wait()
			test_run = 0
			process = []

def main():
	logging.basicConfig(format='%(message)s', level=logging.INFO)
	# download_pdbs(protein_ids)
	run()

if __name__ == '__main__':
	main()
