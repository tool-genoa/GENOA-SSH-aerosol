# -*- coding: utf-8 -*-
#================================================================================
#
#     GENOA v1.0: the GENerator of reduced Organic Aerosol mechanism
#
#     Copyright (C) 2022 CEREA (ENPC) - INERIS.
#     GENOA is distributed under GPL v3.
#
#================================================================================

import os
import json
import time
import numpy as np
from copy import deepcopy

from Parameters import RunSets, Roptions, prefix, pathSSH_rdc, locs, \
                       pathNewChem, tail, pathSSH_rdc, \
                       Ttot, DeltaT, Tnow, initfile, pathInitFiles, \
                       pathNewRes
from AutoTrainingSeries import auto_training_srs
from AutoTesting import auto_testing
from ReductionStrategy import trim_scheme
from DataStream import read_chem_sets,to_SSH_sets,reaction_merge
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
    if  APSetups['from_ref']:
        print('Prepare initialisation files from the reference chemical mechanism.')
        ichem = APSetups['IDchemRef']
        ipath = APSetups['refChemPath']
        ichemfk = ichem+'FA'
        ichempr = ichem+'P'
        chems = [ichem, ichemfk, ichempr]
        # generate folders if not exists
        path_sav_chem = '{0:s}/{1:s}_chem'.format(pathNewChem,tail)
        path_sav_res = '{0:s}/{1:s}'.format(pathNewRes,tail)
        for i in path_sav_chem, path_sav_res: create_folder(i)
        # get chems
        rc,sp=read_chem_sets(speciesfile='{0:s}/{1:s}/{1:s}.mol'.format(ipath,ichem),
                             reactionfile='{0:s}/{1:s}/{1:s}.reactions'.format(ipath,ichem),
                             speciesType='SSH',
                             reactionType='SSH')
        # save ref chem
        to_SSH_sets(path_sav_chem,ichem,rc,sp,0)
        # save fake chem
        rc0,sp0 = deepcopy(rc), deepcopy(sp)
        to_SSH_sets(path_sav_chem,ichemfk,rc0,sp0,0,tag_fake=True)
        # save pre chem
        rc,sp,trim_sp=trim_scheme(rc,sp)
        if trim_sp != []: print('Trimmed useless species from the ref scheme.',trim_sp)
        rc = reaction_merge(rc,sp)
        to_SSH_sets(path_sav_chem,ichempr,rc,sp,0)
        # get soa concentrations
        for i in chems:
            if i == ichemfk:
                icompile = 'complete_compile'
            else: icompile = 'fast_compile'
            run_SSH(i,path_sav_chem,pathSSH_rdc,path_sav_res,
                     i.replace(prefix,''),'Results_{:s}'.format(i),
                     Ttot,DeltaT,[[],locs],Tnow,
                     icompile,initfile,pathInitFiles)
        # record
        APSetups['IDchemFake'] = ichemfk
        APSetups['IDchemPre'] = ichempr
        APSetups['preChemPath'] = path_sav_chem
        APSetups['refPath'] = '{:s}/Results_{:s}'.format(path_sav_res,ichem)
        APSetups['fakePath'] = '{:s}/Results_{:s}'.format(path_sav_res,ichemfk)
        APSetups['prePath'] = '{:s}/Results_{:s}'.format(path_sav_res,ichempr)
    else:
        ichem = None
        print('Read starting information for reduction.')
        if APSetups['prePath'] == None or APSetups['IDchemFake'] == None:
            raise ValueError('No input prePath or IDchemFake. Please use option I of [training].')

    # training in series
    print('Running GENOA training... This process may take few days.')
    IDchem_now, chem_path_now = auto_training_srs(Setups = APSetups,
                                NoLumps = Roptions['FreezeSpecies'],
                                locs = locs)
    t1 = time.perf_counter()
    print('GENOA training is finished. Time used: ', t1-t0)
    print('The name of the final reduced mechanism is {:s}, path: {:s}.'.format(IDchem_now,chem_path_now))

    ## testing
    ATSetups['Trytimes'] = False # try all in the list
    emax_ind, emax, eave, eave_max = auto_testing(Setups = ATSetups, IDchem = IDchem_now, 
                                                  chempath = chem_path_now, savpath = chem_path_now,
                                                  trd_max = 1, out_file = True, sav_soapath = False)
    t2 = time.perf_counter()
    print('GENOA testing is finished. Time used: ', t2-t1, emax_ind, emax, eave, eave_max)    
    print('The average (maximum) errors of testing conditions is {:6.4f} ({:6.4f}).'.format(eave,emax))

    # END
    print('Reduction is finished. Total time: ',t2-t0)

    return None

if __name__ == '__main__':
    auto_reduction()
