import os,sys
import argparse
import json
import subprocess
import shutil

parser = argparse.ArgumentParser(description='prepare crab code')
parser.add_argument('-f', dest='file', default='', help='json file input')
parser.add_argument('-m', dest='mode', default='', help='work mode')
args = parser.parse_args()

def get_abbre(name,sample_type,year):
    if sample_type == 'MC':
        return name.split('/')[1] + '_' + year
    elif sample_type == 'data':
        return name.split('/')[1] + '_' + name.split('/')[2].split('-')[0]

def prepare_crab(name,sample_type,year,era):

    abbre_name = get_abbre(name,sample_type,year) 
    if not os.path.exists('crabcode_' + year):
        os.mkdir("crabcode_" + year)

    print ("------> preparing submit code for",abbre_name)
    with open('crabcode_' + year + '/' + abbre_name + '_cfg.py', 'w+') as f:
        f.write('from WMCore.Configuration import Configuration \n\n')

        f.write('config = Configuration()\n')
        f.write('config.section_("General")\n')
        f.write('config.General.requestName = "' + abbre_name + '"\n')
        f.write('config.General.transferLogs = False \n')
        f.write('config.General.workArea = "crab' + year + '"\n\n')

        f.write('config.section_("JobType")\n')
        f.write('config.JobType.pluginName = "Analysis"\n')
        f.write('config.JobType.psetName = "PSet.py"\n')
        f.write('config.JobType.scriptExe = "./WWG_crab_script.sh" \n')
        f.write('config.JobType.inputFiles = ["../../../../scripts/haddnano.py","../WWG_fakelepton/WWG_postproc.py","../WWG_fakelepton/WWGfakelepton_Module.py","../WWG_fakelepton/WWG_keep_and_drop.txt","../WWG_fakelepton/WWG_output_branch.txt","../WWG_fakelepton/DAS_filesearch.py"] #hadd nano will not be needed once nano tools are in cmssw \n')
#	f.write('config.JobType.scriptArgs = ["isdata=' + sample_type + '","year=' + year + '","era=' + era + '"] \n')
        f.write('config.JobType.scriptArgs = ["isdata=' + sample_type + '","year=' + year + '","era=' + era + '"] \n')
        f.write('config.JobType.sendPythonFolder  = True\n')
        f.write('config.JobType.allowUndistributedCMSSW = True \n\n')

        f.write('config.section_("Data")\n')
        f.write('config.Data.inputDataset = "' + name + '" \n')
        f.write('#config.Data.inputDBS = "phys03"\n')
        f.write('config.Data.inputDBS = "global"\n')
        f.write('#config.Data.splitting = "FileBased"\n')
        f.write('#config.Data.unitsPerJob = 1\n')

        if sample_type == 'MC':
           f.write('config.Data.splitting = "FileBased"\n')
           f.write('config.Data.unitsPerJob = 1\n')
        elif year == '2018':
            f.write('config.Data.splitting = "LumiBased"\n')
            f.write('config.Data.unitsPerJob = 50\n')
            f.write('config.Data.lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt" \n\n')
        elif year == '2017':
            f.write('config.Data.splitting = "LumiBased"\n')
            f.write('config.Data.unitsPerJob = 50\n')
            f.write('config.Data.lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt" \n\n')
        elif year == '2016':
            f.write('config.Data.splitting = "LumiBased"\n')
            f.write('config.Data.unitsPerJob = 50\n')
            f.write('config.Data.lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt" \n\n')

#        f.write('config.Data.outLFNDirBase ="/store/user/sdeng/WWG_analysis/' + sample_type + '/' + year + '"\n')
        f.write('config.Data.publication = False\n')
        f.write('config.Data.ignoreLocality = True\n')
        f.write('config.Data.allowNonValidInputDataset = True\n')
        f.write('config.Data.outputDatasetTag = "' + abbre_name + '" \n\n')

        f.write('config.section_("Site")\n')
        f.write('config.Site.storageSite = "T2_CH_CERN"\n')
        f.write('config.Site.whitelist = ["T2_US_MIT","T2_US_Wisconsin","T2_US_Purdue","T2_US_UCSD","T2_US_Caltech","T2_US_Nebraska"] \n')
        f.close()

def submit(name,sample_type,year):

    abbre_name = get_abbre(name,sample_type,year)

    if not os.path.exists(f'crabcode_{year}/{abbre_name}_cfg.py'):
        print ("crabcode for ",abbre_name," not existed, skipping")
        return True

    r=subprocess.run(args=f"crab submit -c crabcode_{year}/{abbre_name}_cfg.py",shell=True,stdout=subprocess.PIPE,encoding='utf-8')
    if 'Success' in r.stdout:
        print ("--------> submit info:","submit crab jobs for",abbre_name)
    else:
        print ("--------> submit info:","\033[31mfail\033[0m to submit for",abbre_name)

def kill(name,sample_type,year):

    abbre_name = get_abbre(name,sample_type,year)

    if not os.path.exists(f'crab{year}/crab_{abbre_name}'):
        print ("crab log for ",abbre_name," not existed, skipping \n")
        return True

    r=subprocess.run(args=f"crab kill -d crab{year}/crab_{abbre_name}" ,shell=True,stdout=subprocess.PIPE,encoding='utf-8')
    print (r.stdout,'\n')

    shutil.rmtree(f'crab{year}/crab_{abbre_name}')
    print (f'crab{year}/crab_{abbre_name} has been removed')


def status(name,sample_type,year):

    abbre_name = get_abbre(name,sample_type,year)

    if not os.path.exists(f'crab{year}/crab_{abbre_name}'):
        print ("crab log for ",abbre_name," not existed, skipping \n")
        return True

    r=subprocess.run(args=f"crab status -d crab{year}/crab_{abbre_name}" ,shell=True,stdout=subprocess.PIPE,encoding='utf-8')
    print (r.stdout,'\n')


def hadd_help(name,sample_type,year):

    abbre_name = get_abbre(name,sample_type,year)
    store_path = '/eos/user/s/sdeng/WWG_analysis'
    first_name = name.split('/')[1]

    if os.path.exists(f'{abbre_name}.root'):
        print (f'{abbre_name} already existed, skipping')
        return True

    if not (os.path.exists(f'{store_path}/{sample_type}/{year}/{first_name}/{abbre_name}')):
        print (f'results for {abbre_name} not existed in {store_path}/{sample_type}/{year}/{first_name}/{abbre_name}, skipping\n')
        return True
    
    if not (len(os.listdir(f'{store_path}/{sample_type}/{year}/{first_name}/{abbre_name}')) == 1 ):
        print (f'more than 1 result for {abbre_name}, Please check {store_path}/{sample_type}/{year}/{first_name}/{abbre_name}\n')
        return True

    run_number = os.listdir(f'{store_path}/{sample_type}/{year}/{first_name}/{abbre_name}')[0]
    path = f'{store_path}/{sample_type}/{year}/{first_name}/{abbre_name}/{run_number}/0000/'
    print (f'hadding root files in {path}')
    r=subprocess.run(args=f"python $CMSSW_BASE/src/PhysicsTools/NanoAODTools/scripts/haddnano.py {abbre_name}.root {path}/*.root ", shell=True,stdout=subprocess.PIPE,encoding='utf-8')
    
    if os.path.exists(f'{abbre_name}.root'):
        print (f'hadd complete, please check {abbre_name}.root\n')
    else:
        print (f'hadd \033[31mfail\033[0m!!')

def report_lumi(name,sample_type,year):

    abbre_name = get_abbre(name,sample_type,year)
    if not os.path.exists(f'crab{year}/crab_{abbre_name}'):
        print ("crab log for ",abbre_name," not existed, skipping \n")
        return True

    r=subprocess.run(args=f"crab report -d crab{year}/crab_{abbre_name}" ,shell=True,stdout=subprocess.PIPE,encoding='utf-8')
    print (r.stdout,'\n')

    if not os.path.exists(f'lumi_{year}'):
        os.mkdir(f'lumi_{year}')
    
    shutil.copy(f'crab{year}/crab_{abbre_name}/results/notFinishedLumis.json', f'lumi_{year}/{abbre_name}.json')

def resubmit(name,sample_type,year):

    abbre_name = get_abbre(name,sample_type,year)

    if not os.path.exists(f'crab{year}/crab_{abbre_name}'):
        print ("crab log for ",abbre_name," not existed, skipping \n")
        return True

    print (f"resubmitting {abbre_name}\n")
    r = subprocess.run(args=f"crab resubmit -d crab{year}/crab_{abbre_name}" ,shell=True,stdout=subprocess.PIPE,encoding='utf-8')
    print (r.stdout,"\n")


if __name__=='__main__':
    
    with open(args.file, "r") as f:
        jsons = json.load(f)
        f.close()

    if args.mode == 'prepare':
        for dataset in jsons:
            prepare_crab(dataset['name'], dataset['type'], str(dataset['year']), dataset['era'])
    
    if args.mode == 'submit':
        for dataset in jsons:
            submit(dataset['name'], dataset['type'], str(dataset['year']))

    if args.mode == 'kill':
        for dataset in jsons:
            kill(dataset['name'], dataset['type'], str(dataset['year']))
    
    if args.mode == 'status':
        for dataset in jsons:
            status(dataset['name'], dataset['type'], str(dataset['year']))

    if args.mode == 'hadd':
        for dataset in jsons:
            hadd_help(dataset['name'], dataset['type'], str(dataset['year']))

    if args.mode == 'report':
        for dataset in jsons:
            if dataset['type'] == 'data':
                report_lumi(dataset['name'], dataset['type'], str(dataset['year']))
    
    if args.mode == 'resubmit':
        for dataset in jsons:
            resubmit(dataset['name'], dataset['type'], str(dataset['year']))
