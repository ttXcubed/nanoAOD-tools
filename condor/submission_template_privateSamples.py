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

with open('file_list_privateSamples.yaml') as yaml_f:
	try:
	    	yaml_file_dict = yaml.safe_load(yaml_f)
		#yam_file_list = yaml_file.split("\n")
	except yaml.YAMLError as exc:
		print(exc)

if args.sgn:
	listOfDataSets.append(yaml_file_dict['sgn'][0])
	listOutputDir.append(yaml_file_dict['sgn'][1])
		
else:
	for year in yaml_file_dict['bkg'].keys():
		if year != args.year: continue 
		else:
			for key in yaml_file_dict['bkg'][year].keys():
				listOfDataSets.append(yaml_file_dict['bkg'][year][key][0])
				listOutputDir.append(yaml_file_dict['bkg'][year][key][1])
#print(listOfDataSets)
#print(listOutputDir)


mass_values = ['400','600','800','1000','1250','1500','1750','2000','2250','2500','3000']
relwidth_values = ['4','10','20','50']

with open('condor_submission', 'r') as condor_f:
	condor_sub_file = condor_f.read()
	
	if args.sgn:
		condor_sub_file = condor_sub_file.replace('EXE','executable_jobs_sgn_'+args.year+'.sh')
	else: condor_sub_file = condor_sub_file.replace('EXE','executable_jobs_bkg_'+args.year+'.sh')

command_str = ''
command_dict = dict()

#command_list = condor_sub_file
#print(command_list)

batch_number = 1 
batch_count = 0

for index, datasetName in enumerate(listOfDataSets):
	
	if args.sgn:
		for width in relwidth_values:
			for mass in mass_values:
				if args.year=='2017':
				
					if os.path.exists('/pnfs/desy.de/cms/tier2/store/user/gmilella/'+datasetName+width+'/NANOAOD/mass_'+mass):
					
						for file in os.listdir('/pnfs/desy.de/cms/tier2/store/user/gmilella/'+datasetName+width+'/NANOAOD/mass_'+mass):
								
							batch_count += 1
							if batch_count>4999*batch_number:
								command_dict[batch_number] = command_str
								batch_number += 1
								command_str = ''
							command_str += 'arguments = "-i /pnfs/desy.de/cms/tier2/store/user/gmilella/'+datasetName+width+'/NANOAOD/mass_'+mass+'/'+file+' /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer"\nOutput       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer/log/log.$(Process).out\nError       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer/log/log.$(Process).err\nLog       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer/log/log.$(Process).log\nqueue\n '
							command_dict[batch_number] = command_str
							
								
							#condor_f_new.write('arguments = "-i /pnfs/desy.de/cms/tier2/store/user/gmilella/'+datasetName+width+'/NANOAOD/mass_'+mass+'/'+file+' /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer"\n')     
						    
		    					#condor_f_new.write('Output       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer/log/log.$(Process).out\nError       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer/log/log.$(Process).err\nLog       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer/log/log.$(Process).log\nqueue\n')
			
					if not os.path.isdir('/nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer'):
						os.makedirs('/nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer')
						os.makedirs('/nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer/log')
					else: continue
				else: 
					print('/pnfs/desy.de/cms/tier2/store/user/gmilella/'+datasetName+'_'+args.year+'/NANOAOD/mass_'+mass)
					if os.path.exists('/pnfs/desy.de/cms/tier2/store/user/gmilella/'+datasetName+'_'+args.year+'/NANOAOD/mass_'+mass):
						print('si')
						for file in os.listdir('/pnfs/desy.de/cms/tier2/store/user/gmilella/'+datasetName+width+'_'+args.year+'/NANOAOD/mass_'+mass):
							print(file)
							batch_count += 1
							if batch_count>4999*batch_number:
								command_dict[batch_number] = command_str
								batch_number += 1
								command_str = ''
							command_str += 'arguments = "-i /pnfs/desy.de/cms/tier2/store/user/gmilella/'+datasetName+width+'/NANOAOD/mass_'+mass+'/'+file+' /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer"\nOutput       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer/log/log.$(Process).out\nError       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer/log/log.$(Process).err\nLog       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer/log/log.$(Process).log\nqueue\n '
							command_dict[batch_number] = command_str
							
							#condor_f_new.write('arguments = "-i /pnfs/desy.de/cms/tier2/store/user/gmilella/'+datasetName+width+'_'+args.year+'/NANOAOD/mass_'+mass+'/'+file+' /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer"\n')     
						    
		    					#condor_f_new.write('Output       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer/log/log.$(Process).out\nError       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer/log/log.$(Process).err\nLog       = /nfs/dust/cms/user/gmilella/sgn_'+args.year+'/'+listOutputDir[index]+'_mass'+mass+'_width'+width+'_ntuplizer/log/log.$(Process).log\nqueue\n')
			

					else: 
						print('no')
						continue	
					
	else:

		if os.path.exists('/pnfs/desy.de/cms/tier2/store/user/gmilella/'+datasetName):

			for subdir, dirs, files in os.walk('/pnfs/desy.de/cms/tier2/store/user/gmilella/'+datasetName):
			    for file in files:
				#print os.path.join(subdir, file)
				filepath = subdir + os.sep + file

				batch_count += 1
				if batch_count>4999*batch_number:
					command_dict[batch_number] = command_str
					batch_number += 1
					command_str = ''
				command_str += 'arguments = "-i '+filepath+' /nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'"\nOutput       = /nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'/log/log.$(Process).out\nError       = /nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'/log/log.$(Process).err\nLog       = /nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'/log/log.$(Process).log\nqueue\n'
				command_dict[batch_number] = command_str
				
				#condor_f_new.write('arguments = "-i '+filepath+' /nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'"\n')
				#condor_f_new.write('Output       = /nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'/log/log.$(Process).out\nError       = /nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'/log/log.$(Process).err\nLog       = /nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'/log/log.$(Process).log\nqueue\n')
				
		if not os.path.isdir('/nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]):
				os.makedirs('/nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index])
				os.makedirs('/nfs/dust/cms/user/gmilella/bkg_'+args.year+'/'+listOutputDir[index]+'/log')
		else: continue	

#print(batch_number)
#print(batch_count)
#print(command_dict.keys())
#print(command_dict[1])


for batch in command_dict.keys():
	with open('condor_submission_new_'+str(batch), 'w+') as condor_f_new: 
		condor_f_new.write(condor_sub_file)
		condor_f_new.write(command_dict[batch])	
	
			
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


     
