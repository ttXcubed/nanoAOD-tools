import os
import subprocess
import sys
import argparse
import yaml
from collections import OrderedDict

parser = argparse.ArgumentParser()

parser.add_argument('--sgn', dest='sgn', action='store_true', help='Signal', default=False)

args = parser.parse_args()

subprocess.call("voms-proxy-init --rfc --voms cms -valid 192:00", shell = True)
subprocess.call("echo $X509_USER_PROXY", shell = True)
subprocess.call("export X509_USER_PROXY=/afs/desy.de/user/g/gmilella/.globus/x509", shell = True)

listOfDataSets = []
listOutputDir = []
yaml_file_dict = {}

with open('file_list_NN.yaml') as yaml_f:
	try:
	    	yaml_file_dict = yaml.safe_load(yaml_f)
		#yam_file_list = yaml_file.split("\n")
	except yaml.YAMLError as exc:
		print(exc)

for key in yaml_file_dict.keys():
	listOfDataSets.append(yaml_file_dict[key][0])
	listOutputDir.append(yaml_file_dict[key][1])

#print(listOfDataSets)
#print(listOutputDir)

dasGoCommand = 'dasgoclient -query="file dataset=file dataset={DATASET} {INSTANCE}"'

queueStr = ""


with open('condor_submission', 'r') as condor_f:
	condor_sub_file = condor_f.read()
	condor_sub_file = condor_sub_file.replace('EXE','executable_jobs_tmp_topNN.sh')

with open('condor_submission_new', 'w+') as condor_f_new: 
	condor_f_new.write(condor_sub_file)	

	for index, datasetName in enumerate(listOfDataSets):
		if 'USER' in datasetName:
			listOfFiles = subprocess.check_output(dasGoCommand.format(DATASET=datasetName,INSTANCE="instance=prod/phys03"), shell = True).split("\n")[:-1]
		else:	
			#dasGoCommand += "
			listOfFiles = subprocess.check_output(dasGoCommand.format(DATASET=datasetName,INSTANCE=" "), shell = True).split("\n")[:-1]
		
		
		for i, inFile in enumerate(listOfFiles):
		
			if os.path.exists('/pnfs/desy.de/cms/tier2/'+inFile):
				#print("capocchione al sugo ", inFile)
				condor_f_new.write('arguments = "-i /pnfs/desy.de/cms/tier2/'+inFile+' /nfs/dust/cms/user/gmilella/topNN'+listOutputDir[index]+'"\n')
				
			else:
				#print("salsizz arraganat ", inFile)
				condor_f_new.write('arguments = "-i root://cms-xrd-global.cern.ch/'+inFile+' /nfs/dust/cms/user/gmilella/topNN'+listOutputDir[index]+'"\n')
				
			condor_f_new.write('Output       = /nfs/dust/cms/user/gmilella/topNN'+listOutputDir[index]+'/log/log.$(Process).out\nError       = /nfs/dust/cms/user/gmilella/topNN'+listOutputDir[index]+'/log/log.$(Process).err\nLog       = /nfs/dust/cms/user/gmilella/topNN'+listOutputDir[index]+'/log/log.$(Process).log\nqueue\n')

			
			if not os.path.isdir("/nfs/dust/cms/user/gmilella/topNN"+listOutputDir[index]):
				os.makedirs("/nfs/dust/cms/user/gmilella/topNN"+listOutputDir[index])
				os.makedirs("/nfs/dust/cms/user/gmilella/topNN"+listOutputDir[index]+"/log")
			else: continue
			
				
			

		

     
