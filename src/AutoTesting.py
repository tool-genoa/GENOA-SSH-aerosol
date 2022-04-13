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
import sys
import json
import numpy as np
from datetime import datetime

from Functions import isint, get_locations
from SSHSetUp import set_SSH,get_namelist,run_SSHsimulation
from Parameters import RunSets, pathNewChem, pathNewRes, \
                       pathSSH_rdc, pathSSH_ref,\
                       IDchem, Ttot, DeltaT, Tnow, prefix, \
                       pathInitFiles, initfile, namelist_pre

def auto_testing(Setups = RunSets[4], IDchem = IDchem, chempath = pathNewChem,
                    ind = [], trd_max = 1, out_file = True, savpath = pathNewChem, 
                    pathTest = pathInitFiles, sav_soapath = False):
    """
        Objectives:
            Test the input reduced chemical scheme on the initial conditions extracted from CHIMERE model

        Inputs:
            
        Outputs:
    """
    # Part A: Initialisation
    # Assign parameters
    Test_file = Setups['Test_file']

    # error criterias
    err_max = Setups['err_max']
    err_out = Setups['err_out']

    # current case
    tail = IDchem.replace(prefix,'')
    if prefix not in IDchem: IDchem = prefix + IDchem
    ResultFolder = 'Results_{:s}'.format(IDchem)

    # ref case
    IDchem_ref = Setups['IDchemRef']
    refPath = Setups['refChemPath']

    orgPath = Setups['orgpath'] # the folder name if save ref case for multiple reduction
    if not orgPath and sav_soapath:
        orgPath = 'ref'

    tail_ref = IDchem_ref.replace(prefix,'')

    ResultFolder_ref = 'Results_{:s}'.format(IDchem_ref)

    if ind == []:
        ## Obtain all locs for test cases
        if Test_file is not None:
            if not os.path.exists(Test_file):
                print('Not find Test_file.', Test_file)
                raise FileNotFoundError(Test_file)
            with open (Test_file,'r') as f:
                ind = json.loads(f.read())
        else:
            raise OSError('not provied Test_file.')

        if isint(Setups['Trytimes']): # number of testing. use for general random tests
            Trytimes = int(Setups['Trytimes'])
            if Trytimes < len(ind) and Trytimes > 0: 
                ind = ind[0:Trytimes]

    # record
    if out_file:
        ## initialize record files
        os.makedirs(savpath, exist_ok=True)
        recordfile = '{:s}/Testing_{:s}'.format(savpath,IDchem)
        print('Save Testing file: ', recordfile)
        # reopen for the output
        fall = open (recordfile+'_all','w+') # record all info
        # output info in record files
        tmp = 'Testing. Err max {:f} Try inds: {:d}\n'.format(err_max,len(ind))

        if isinstance(orgPath, str):
            tmp += 'PathTest:{:s}, IDchem_ref:{:s}, IDchem_ref:{:s}, orgPath: {:s}, Namelist:{:s}\n'.format(pathTest, IDchem_ref, IDchem, orgPath, namelist_pre)
        else:
            tmp += 'PathTest:{:s}, IDchem_ref:{:s}, IDchem_ref:{:s}, orgPath:None, Namelist:{:s}\n'.format(pathTest, IDchem_ref, IDchem, namelist_pre)
        fall.write(tmp)

    else: fall = ''

    # compile SSH-aerosol for reduction case
    set_SSH(IDchem,chempath,pathSSH_rdc)
    set_SSH(IDchem_ref,refPath,pathSSH_ref)

    #Namelist = '{:s}/namelist_test.ssh'.format(pathSSH_rdc)
    ssh_namelist, ssh_inds = get_namelist()

    ## run test cases
    # init variables
    nall = 0 # record the total number of test cases
    nout = 0 # record the number of cases with err > err_max
    emax = 0.0 # record current max error
    eall =0.0 # sum up all error to compute the average error
    ealls = {} # all error pre Tnow
    for inow in Tnow: ealls['{:d}h'.format(inow)] = 0.0

    emax_case = '' # record the case with current max error
    out_ind = [] # save the loc with max err
    
    itrd = 0
    nind = len(ind)

    pathSSH_ref_abs = os.path.abspath(pathSSH_ref)

    while itrd < nind:

        if trd_max == 1:
            trd_val = nind - itrd
        elif trd_max == 2: # two times
            if itrd == 0: trd_val = 12
            else: trd_val = nind - itrd
        else: #average
            if itrd == 0: trd_val = int(nind/trd_max)
            
        ## 1. initialization of init files, locs, locations
        if trd_val + itrd <= nind:
            locs = ind[itrd:]
        else:
            locs = ind[itrd:trd_val + itrd]

        itrd += trd_val
        if initfile == 'storage': 
            ifile = 'storage'
            locations = get_locations(locs)

        ## 2. Run simulations in SSH-aerosol
        # reference case
        if orgPath: # check; if non, run reference case and save Organics_1.txt in the format: m/y/x/0h.txt
            tag_savorg = 0 # if need to saved ref. 0: default ref exists
            results_ref = [] # absolute path!
            for il in locs:
                for inow in Tnow:
                    tmp = '{:s}/{:s}/m{:d}/y{:d}/x{:d}/{:d}h.txt'.format(pathSSH_ref_abs,orgPath,il[2],il[0],il[1],inow)
                    if not os.path.exists(tmp): tag_savorg += 1
                    else: results_ref.append(tmp)
                    if tag_savorg: break

        else: tag_savorg = 1 # no previous saved ref

        if tag_savorg: # run ref cases
            res_ref_now, errs = run_SSHsimulation(tail_ref,ResultFolder_ref,pathSSH_ref,
                                  ssh_namelist, ssh_inds,
                                  Ttot,DeltaT,
                                  [locations,locs],
                                  Tnow, pathInit = pathTest,
                                  initfile = ifile,
                                  savorg = orgPath)
            results_ref = ['{:s}/{:s}'.format(pathSSH_ref_abs, i) for i in res_ref_now]
        #else:
            #print('use the ref files in {:s}/{:s}'.format(pathSSH_ref,orgPath)) # print to check

        # reduction case
        if sav_soapath: isavorg = 'rdc'
        else: isavorg = False
        results_rdc, errs = run_SSHsimulation(tail,ResultFolder,pathSSH_rdc,
                                              ssh_namelist, ssh_inds,
                                              Ttot,DeltaT,[locations,locs],
                                              Tnow,pathInit = pathTest,
                                              initfile = ifile, 
                                              ref_file = results_ref,
                                              savorg = isavorg,
                                              ioutType = 1)

        ## 3. Analysis results
        errs = [i[0] for i in errs]

        n = len(errs)
        err_now = np.max(errs)# current max one
        # save info if need
        if out_file:
            for i in range(n):
                if errs[i] > err_max: nout += 1
                fall.write('Run:{:d}\tOut:{:d}\t{:s}\t{:6.4f}\n'.format(i+nall+1,nout,results_rdc[i].replace(ResultFolder+'/',''), errs[i]))
            fall.write('Current Max: {:6.4f}\n'.format(err_now))
            fall.flush()

        # record loc and error
        eall += np.sum(errs)
        nall += len(errs)
        for i in range(n):
            inow = results_rdc[i].split('_')[-2] # shoule be Xh
            ealls[inow] += errs[i]

        if err_now > emax: # compare to the current recorded max err
            emax = err_now
            y,x,m = locs[int(np.where(errs == err_now)[0][0] / (len(Tnow)*len(Ttot)))] # gedit ix
            out_ind = [[y,x,m]]
            emax_case = 'm{:d}y{:d}x{:d}'.format(m,y,x)
            
        if err_now > err_out: # compare to the criteria
            break

    # sum
    # average error
    err_ave = eall/nall
    if len(ealls) != len(Tnow):
        print(ealls, Tnow, len(ealls), len(Tnow))
        raise ValueError('AT: number in ealls not fits number in Tnow.')
    nt = nall/len(ealls)
    err_ave_max = max(ealls.values())/nt
    
    if out_file:
        tmp = 'Ave error pre Tnow:'
        for inow in Tnow: 
            err_tmp = ealls['{:d}h'.format(inow)]/nt
            tmp += ' {:d} {:6.4f};'.format(inow, err_tmp)

        fall.write('-----finish-----total run {:d}-----out run {:d}-----\n Case with max err: {:s}. Ave err {:6.4f}. Ave err max {:6.4f}\n'.format(nall,nout,emax_case,err_ave, err_ave_max))
        fall.write(tmp)
        fall.close()
        #print('AT: ', out_ind,emax)

    # mv soa results
    if sav_soapath:
        path = '{:s}/{:s}_test'.format(sav_soapath, tail)
        create_folder(path,del_exist = True)
        # save rdc results
        os.system('mv {:s}/rdc {:s}/'.format(pathSSH_rdc,path))
        # save ref results
        os.system('mv {:s}/{:s} {:s}/'.format(pathSSH_ref,orgPath,path))

    if out_ind == []:
        print('Testing result []',out_ind, emax, err_ave, err_ave_max,tail_ref,pathSSH_ref,tail,pathSSH_rdc)
    return out_ind, emax, err_ave, err_ave_max

if __name__ == '__main__':
    auto_testing()
