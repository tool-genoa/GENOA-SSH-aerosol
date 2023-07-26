# -*- coding: utf-8 -*-
#================================================================================
#
#     GENOA v1.0: the GENerator of reduced Organic Aerosol mechanism
#
#     Copyright (C) 2023 CEREA (ENPC) - INERIS.
#     GENOA is distributed under GPL v3.
#
#================================================================================

import os
import json
import time
import shutil
import numpy as np
from copy import deepcopy

from Parameters import RunSets, Roptions, prefix, pathSSH_rdc, locs, \
                       pathNewChem, tail, pathSSH_rdc, \
                       Ttot, DeltaT, Tnow, initfile, pathInitFiles, \
                       pathNewRes, ipvocs
from AutoTrainingSeries import auto_training_srs
from AutoTrainingParallel import auto_training_para
from AutoTesting import auto_testing
from ReductionStrategy import trim_scheme, reaction_merge
from DataStream import read_chem_sets,to_SSH_sets
from SSHSetUp import run_SSH, set_SSH
from Functions import create_folder

def auto_reduction(APSetups = RunSets[3], ATSetups = RunSets[4], locs = locs):
    """
        Objectives:
            generate a full set of reduction from given training dataset.

        Inputs:
            parameters are read from Parameter.py            
        Outputs:
            a record of reduction with all valiad reduced schemes
    """
    # record total reduction time
    t0 = time.perf_counter()
    print('\n\nStart reduction.')

    # check restart: update values for AP
    if APSetups['restart']:
        refile = None
        path_sav_rec = '{0:s}/{1:s}_recs'.format(pathNewChem,tail)
        if APSetups['restart'] != 'optional' and isinstance(APSetups['restart'], str):
            refile = APSetups['restart']
        else:
            if os.path.exists(path_sav_rec):
                print('Try to search the restart file from ', path_sav_rec,' ...')
                refiles = [i for i in os.listdir(path_sav_rec) if 'restart.' in i]
                if refiles != []:
                    n = max([int(i.split('.')[1]) for i in refiles])
                    refile = '{:s}/restart.{:d}'.format(path_sav_rec,n)

        if refile is None or not os.path.exists(refile):
            print('Restart file is not found in ',path_sav_rec, refile)
            if APSetups['restart'] == 'optional':
                print('start reduction with no restart.')
                APSetups['restart'] = False
            else:
                raise FileNotFoundError('Stop program cuz restart file does not find.',APSetups['restart'])
        else:
            print('Find restart file: ',refile)
            APSetups['restart'] = refile
            with open(refile, 'r') as f: re_data = json.loads(f.read())
            for i in re_data:
                if i in list(APSetups.keys()) and re_data[i] != APSetups[i]:
                    # find changes
                    print('--- Change APSetups[',i,'] from ',APSetups[i],' to ',re_data[i])
                    APSetups[i] = re_data[i]
        tag_reuse = True
    else:
        # reuse generated files
        if APSetups['Reuse_chem']: tag_reuse = True
        else: tag_reuse = False
        
    # update parameters
    if APSetups['IDchemRef'] is None or APSetups['refChemPath'] is None:
        raise ValueError('Information is missing for reduction. check IDchemRef & refChemPath.')

    if ATSetups['IDchemRef'] is None:
        ATSetups['IDchemRef'] = APSetups['IDchemRef']
    if ATSetups['refChemPath'] is None:
        ATSetups['refChemPath'] = APSetups['refChemPath']

    # prepare SSH
    set_SSH(None,None,pathSSH_rdc,mode='fast')
    # clean up old files in SSH path
    cmd = 'rm -rf {:s}/Results_*;'.format(pathSSH_rdc)
    for j in os.listdir('{:s}/src/include/CHEMISTRY/'.format(pathSSH_rdc)):
        if prefix in j and j != 'BCARY':
            cmd += 'rm -rf {:s}/src/include/CHEMISTRY/{:s};'.format(pathSSH_rdc,j)
    os.system(cmd)

    # prepare ref concs.
    if APSetups['IDchemRef'] and APSetups['refChemPath']:
        ichem = APSetups['IDchemRef']
        ipath = APSetups['refChemPath']
        print('Read IDchemRef and refChemPath: ',ichem,ipath)
        
        # check format
        if os.path.exists('{0:s}/{1:s}/{1:s}.mol'.format(ipath,ichem)):
            ipath = '{0:s}/{1:s}/{1:s}'.format(ipath,ichem)
        elif os.path.exists('{0:s}/{1:s}.mol'.format(ipath,ichem)):
            ipath = '{0:s}/{1:s}'.format(ipath,ichem)
            print('IDchemRef format: ipath/ichem.reactions ',ipath)
        else:
            raise FileNotFoundError('IDchemRef .mol file is not found.')
    else:
        raise ValueError('IDchemRef/refChemPath can not be None. Please define in configuration file.', APSetups['IDchemRef'],APSetups['refChemPath'])
        
    path_sav_chem = '{0:s}/{1:s}_chems'.format(pathNewChem,tail)
    path_sav_res = '{0:s}/{1:s}'.format(pathNewRes,tail)

    if tag_reuse:
        # update chem names
        if APSetups['IDchemFake'] == None: APSetups['IDchemFake'] = ichem+'FA'
        if APSetups['IDchemPre'] == None: APSetups['IDchemPre'] = ichem+'P'
        
        for i in APSetups['IDchemRef'], APSetups['IDchemFake'],  APSetups['IDchemPre']:
            # check chems exist
            if not os.path.exists('{0:s}/{1:s}/{1:s}.mol'.format(path_sav_chem,i)):
                raise FileNotFoundError('Can not use Resue_chem option. Chems file is missing: ',path_sav_chem,i)
            # check results exist
            if not os.path.exists('{:s}/Results_{:s}'.format(path_sav_res,i)):
                raise FileNotFoundError('Can not use Resue_chem option. Results file is missing: ',path_sav_res,i)

    else:
        # generate folders if not exists
        for i in path_sav_chem, path_sav_res: create_folder(i)

        # read chems
        rc,sp=read_chem_sets(speciesfile='{:s}.mol'.format(ipath),
                             reactionfile='{:s}.reactions'.format(ipath),
                             speciesType='SSH',
                             reactionType='SSH')
        # save ref chem
        to_SSH_sets(path_sav_chem,ichem,rc,sp,'10')

    # change refChemPath
    APSetups['refChemPath'] = path_sav_chem


    # generate Fake chem
    if APSetups['IDchemFake'] == None: 
        ichemfk = ichem+'FA'
        if tag_reuse:
            print('Reuse IDchemFake: ',ichemfk)
        else:
            # save fake chem
            rc0,sp0 = deepcopy(rc), deepcopy(sp)
            to_SSH_sets(path_sav_chem,ichemfk,rc0,sp0,'10',tag_fake=True)
            print('Generate IDchemFake from IDchemRef: ',ichemfk)
        APSetups['IDchemFake'] = ichemfk
    else: 
        ichemfk = APSetups['IDchemFake']
        print('Read IDchemFake: ',ichemfk)

    # generate/read pre chem
    if APSetups['IDchemPre'] == None:
        ichempr = ichem+'P'
        # generate pre chem
        rc,sp,trim_sp=trim_scheme(rc,sp)
        if trim_sp != []: print('Trimmed useless species from the ref scheme.',trim_sp)
        rc = reaction_merge(rc,sp)
        to_SSH_sets(path_sav_chem,ichempr,rc,sp,'10')
        print('Generate IDchemPre from IDchemRef: ',ichempr)
        APSetups['IDchemPre'] = ichempr
        
    else: # check prechem and prepath
        ichempr = APSetups['IDchemPre']
        ipath = APSetups['preChemPath']
        if ipath == None: 
            raise ValueError('preChemPath should not be None with given IDchemPre.',ichempr, ipath)
        else:
            print('Read IDchemPre: ',ichempr)
            # check format: ichem.reactions or ichem/ichem.reactions
            if os.path.exists('{0:s}/{1:s}/{1:s}.mol'.format(ipath,ichempr)):
                # clean old if need
                if ipath != path_sav_chem:
                    shutil.rmtree('{:s}/{:s}'.format(path_sav_chem,ichempr), ignore_errors=True)
                    shutil.copytree('{:s}/{:s}'.format(ipath,ichempr),
                                '{:s}/{:s}'.format(path_sav_chem,ichempr))
            elif os.path.exists('{0:s}/{1:s}.mol'.format(ipath,ichempr)):
                ipath = '{0:s}/{1:s}'.format(ipath,ichempr)
                print('IDchemPre format: ipath/ichem.reactions ',ipath)
                # read pre chems
                rc,sp=read_chem_sets(speciesfile='{:s}.mol'.format(ipath),
                         reactionfile='{:s}.reactions'.format(ipath),
                         speciesType='SSH',
                         reactionType='SSH')
                # sav
                to_SSH_sets(path_sav_chem,ichempr,rc,sp,'10')
            else:
                raise FileNotFoundError('IDchemPre .mol file is not found.')
    # change preChemPath
    APSetups['preChemPath'] = path_sav_chem

    # get concs
    chems = []
    if APSetups['fakePath'] == None: 
        APSetups['fakePath'] = '{:s}/Results_{:s}'.format(path_sav_res,ichemfk)
        if ichemfk not in chems: chems.append(ichemfk)

    if APSetups['refPath'] == None: 
        APSetups['refPath'] = '{:s}/Results_{:s}'.format(path_sav_res,ichem)
        if ichem not in chems: chems.append(ichem)

    if APSetups['prePath'] == None: 
        APSetups['prePath'] = '{:s}/Results_{:s}'.format(path_sav_res,ichempr)
        if ichempr not in chems: chems.append(ichempr)

    if not tag_reuse:
        # get soa concentrations
        for i in chems: #chems = [ichemfk, ichem, ichempr]
            # check if it exists
            
            if i == ichemfk:
                icompile = 'complete_compile'
                j = 0
            else: 
                icompile = 'fast_compile'
                j = ipvocs
            print('Run simulations for chem: ',i,' with mode: ',icompile)
            run_SSH(i,path_sav_chem,pathSSH_rdc,path_sav_res,
                     i.replace(prefix,''),'Results_{:s}'.format(i),
                     Ttot,DeltaT,[[],locs],Tnow,
                     icompile,initfile,pathInitFiles,
                     ipvoc=j)

    # reduction style
    if  APSetups['tag_para']: #  parallel
        print('Reduction in parallel.')
        
        # generate concs for pre-teasting conditions
        if ichem in chems and not tag_reuse:
            print('Run reference and/or previous valid mechanisms.')

            # get pre-testing conditions
            nPreTest = APSetups['nPreTest']
            if isinstance(nPreTest, str):
                # read from given file
                if not os.path.exists(nPreTest):
                    print('Not find nPreTest file.', nPreTest)
                    raise FileNotFoundError(nPreTest)
                else:
                    print('read nPreTest from file: ', nPreTest)
                    with open (nPreTest,'r') as f:
                        prelocs = json.loads(f.read())
                        
            elif isinstance(nPreTest, int):
                # read from testing file
                Test_file = ATSetups['Test_file']
                if Test_file is not None:
                    if not os.path.exists(Test_file):
                        print('Not find Test_file.', Test_file)
                        raise FileNotFoundError(Test_file)
                    else:
                        print('read nPreTest: ', nPreTest, 'from file ', Test_file)
                        with open (Test_file,'r') as f:
                            prelocs = json.loads(f.read())[0:nPreTest]
                else:
                    raise OSError('Not provied Test_file.',nPreTest)
                    
            run_SSH(ichem,path_sav_chem,pathSSH_rdc,path_sav_res,
                     ichem.replace(prefix,''),'Results_{:s}'.format(ichem),
                     Ttot,DeltaT,[[],prelocs],Tnow,
                     'fast_compile',initfile,pathInitFiles,
                     ipvoc=ipvocs)

        # prepare a version of SSH
        print('Setting default SSH-aerosol in ',pathSSH_rdc+'/ssh')
        set_SSH(APSetups['IDchemRef'],APSetups['refChemPath'],pathSSH_rdc+'/ssh')
        # update compile option to speed up the reduction
        with open (pathSSH_rdc+'/ssh/compile','r') as f: info=f.read().splitlines()
        with open (pathSSH_rdc+'/ssh/compile','w+') as f:
            for i in info:
                if 'scons' in i and i[0] == "s": f.write('scons -j8 --implicit-cache\n') #-Q
                else:  f.write(i+'\n')
        #print('Set new settings: scons -j8 --implicit-cache in compile.')
        
        # auto_training
        print('\nStart training.')
        IDchem_now, chem_path_now = auto_training_para(Setups = APSetups,
                                    locs = locs)
    else: # training in series
        print('Reduction in series.')
        # auto_training
        print('Start training.')
        IDchem_now, chem_path_now = auto_training_srs(Setups = APSetups,
                                    locs = locs)
    t1 = time.perf_counter()
    print('GENOA training is finished. Time used: ', t1-t0)

    if IDchem_now:
        print('The name of the final reduced mechanism is {:s}, path: {:s}.'.format(IDchem_now,chem_path_now))

    ## testing
    if IDchem_now:
        print('Start testing, save the results by IDname ', IDchem_now)
        #ATSetups['loc_num'] = False # try all in the list
        #ATSetups['savname'] = IDchem_now # save rdc concs
        emax_ind, emax, eave, eave_max = auto_testing(Setups = ATSetups, IDchem = IDchem_now, 
                                                  chempath = chem_path_now, savpath = chem_path_now,
                                                  trd_max = 1, out_file = True, sav_soapath = True,
                                                  ipvoc = ipvocs) # testing is without ipvocs
        t2 = time.perf_counter()
        print('GENOA testing is finished. Time used: ', t2-t1, emax_ind, emax, eave, eave_max)    
        print('The average (maximum) errors of testing conditions is {:6.4f} ({:6.4f}).'.format(eave,emax))
    else:
        print('No training as no IDchem is found.')

    # clean up
    #print('Clean up tmp files.')
    # remove Results_ folder, chem, 
    
    # END
    t3 = time.perf_counter()
    print('Reduction is finished. Total time: ',t3-t0)

    return None

if __name__ == '__main__':
    auto_reduction()
