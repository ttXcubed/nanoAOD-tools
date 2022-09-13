import os
import subprocess
import sys
import argparse
import yaml
from collections import OrderedDict

parser = argparse.ArgumentParser()

parser.add_argument('--sgn', dest='sgn', action='store_true', help='Signal', default=False)
parser.add_argument('--year', dest='year', type=str, help='Year', default='2017')

args = parser.parse_args()

#subprocess.call("voms-proxy-init --rfc --voms cms -valid 192:00", shell = True)
#subprocess.call("echo $X509_USER_PROXY", shell = True)
#subprocess.call("export X509_USER_PROXY=/afs/desy.de/user/g/gmilella/.globus/x509", shell = True)

listOfDataSets = []
listOutputDir = []
yaml_file_dict = {}

with open('file_list.yaml') as yaml_f:
	try:
	    	yaml_file_dict = yaml.safe_load(yaml_f)
		#yam_file_list = yaml_file.split("\n")
	except yaml.YAMLError as exc:
		print(exc)

#print(yaml_file_dict)
#print(type(yaml_file_dict))

if args.sgn:
	for key in yaml_file_dict.keys():
		listOfDataSets.append(yaml_file_dict[key][0])
		listOutputDir.append(yaml_file_dict[key][1])	

else:
	for year in yaml_file_dict.keys():
		if year != args.year: continue 
		else:
			for key in yaml_file_dict[year].keys():
				listOfDataSets.append(yaml_file_dict[year][key][0])
				listOutputDir.append(yaml_file_dict[year][key][1])
#print(listOfDataSets)
#print(listOutputDir)


dasGoCommand = 'dasgoclient -query="file dataset=file dataset={DATASET} {INSTANCE}"'

queueStr = ""


with open('condor_submission', 'r') as condor_f:
	condor_sub_file = condor_f.read()
	
	if args.sgn:
		condor_sub_file = condor_sub_file.replace('EXE','executable_jobs_sgn_'+args.year+'.sh')
	else: condor_sub_file = condor_sub_file.replace('EXE','executable_jobs_bkg_'+args.year+'.sh')

with open('condor_submission_new', 'w+') as condor_f_new: 
	condor_f_new.write(condor_sub_file)	

	for index, datasetName in enumerate(listOfDataSets):
		if 'USER' in datasetName:
			listOfFiles = subprocess.check_output(dasGoCommand.format(DATASET=datasetName,INSTANCE="instance=prod/phys03"), shell = True).split("\n")[:-1]
			if listOfFiles == []:
				print("WRONG DATASET: ", datasetName)
		else:	
			#dasGoCommand += "
			listOfFiles = subprocess.check_output(dasGoCommand.format(DATASET=datasetName,INSTANCE=" "), shell = True).split("\n")[:-1]
			if listOfFiles == []:
				print("WRONG DATASET: ", datasetName)

		
		for i, inFile in enumerate(listOfFiles):
			#print(inFile)
			if args.sgn: 
			
				if os.path.exists('/pnfs/desy.de/cms/tier2/'+inFile):
					condor_f_new.write('arguments = "-i /pnfs/desy.de/cms/tier2/'+inFile+' /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'"\n')	
				else: 
					condor_f_new.write('arguments = "-i root://cms-xrd-global.cern.ch/'+inFile+' /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'"\n')
			
				condor_f_new.write('Output       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'/log/log.$(Process).out\nError       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'/log/log.$(Process).err\nLog       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'/log/log.$(Process).log\nqueue\n')
				
				if not os.path.isdir('/nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]):
					os.makedirs('/nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index])
					os.makedirs('/nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+"/log")
				else: continue
				
			else: 
			
				if os.path.exists('/pnfs/desy.de/cms/tier2/'+inFile):
					condor_f_new.write('arguments = "-i /pnfs/desy.de/cms/tier2/'+inFile+' /nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'"\n')	
				else: 
					condor_f_new.write('arguments = "-i root://cms-xrd-global.cern.ch/'+inFile+' /nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'"\n')

			
				condor_f_new.write('Output       = /nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'/log/log.$(Process).out\nError       = /nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'/log/log.$(Process).err\nLog       = /nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'/log/log.$(Process).log\nqueue\n')
				
				if not os.path.isdir('/nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]):
					os.makedirs('/nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index])
					os.makedirs('/nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+"/log")
				else: continue
			
if args.sgn:
	with open('executable_jobs_sgn_tmp.sh', 'r') as exe_f_tmp:
		exe_file = exe_f_tmp.read()
else: 
	with open('executable_jobs_bkg_tmp.sh', 'r') as exe_f_tmp:
		exe_file = exe_f_tmp.read()
exe_file=exe_file.replace('YEAR', args.year)
					
if args.sgn:
	with open('executable_jobs_sgn_'+args.year+'.sh', 'w') as exe_f:
		exe_f.write(exe_file)
else:
	with open('executable_jobs_bkg_'+args.year+'.sh', 'w') as exe_f:
		exe_f.write(exe_file)		

		

     
