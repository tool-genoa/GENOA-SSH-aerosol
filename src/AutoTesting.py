# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  AutoTesting.py runs pre-testing/testing processes under given conditions. 
#
# ================================================================================

import os
import sys
import json
import shutil
import numpy as np

from datetime import datetime

from Functions import isint, get_locations
from SSHSetUp import set_SSH,get_namelist,run_SSHsimulation
from Parameters import RunSets, pathNewChem, pathNewRes, \
                       pathSSH_rdc, pathSSH_ref,\
                       IDchem, Ttot, DeltaT, Tnow, prefix, \
                       pathInitFiles, initfile, namelist_pre, \
                       ipvocs, primaryVOCs
                       
def auto_testing(Setups = RunSets[4], IDchem = IDchem, chempath = pathNewChem,
                    ind = [], trd_max = 1, out_file = True, savpath = pathNewChem, 
                    pathTest = pathInitFiles, sav_soapath = False, ipvoc=ipvocs):
    """
        Objectives:
            Test the input reduced chemical scheme on the initial conditions extracted from CHIMERE model

    """
    # Part A: Initialisation
    # Assign parameters
    Test_file = Setups['Test_file']

    # error criterias
    err_max = Setups['err_max']
    err_out = Setups['err_out']

    # savname for output file
    if isinstance(Setups['savname'],str): savname = '_'+Setups['savname']
    else: savname = '_all'

    # current case
    tail = IDchem.replace(prefix,'')
    if prefix not in IDchem: IDchem = prefix + IDchem
    ResultFolder = 'Results_{:s}'.format(IDchem)

    # ref case
    IDchem_ref = Setups['IDchemRef']
    refPath = Setups['refChemPath']

    orgPath = Setups['orgpath'] # the folder name if save ref case for multiple reduction

    # output SOA yield or SOA concs.
    if '_wgas' in namelist_pre:
        tag_wgas = True
        savname += '_SOAyield'
    else:
        tag_wgas = False
            
    if not orgPath and sav_soapath:
        if savname == '_all': orgPath = 'ref'
        else: orgPath = 'ref'+savname

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

        if isint(Setups['loc_num']): # number of testing. use for general random tests
            loc_num = int(Setups['loc_num'])
            if loc_num < len(ind) and loc_num > 0: 
                ind = ind[0:loc_num]
                print('AT loc_num: Testing on ',loc_num,' conditions.')

    # record
    if out_file:
        ## initialize record files
        os.makedirs(savpath, exist_ok=True)
        recordfile = '{:s}/Testing_{:s}'.format(savpath,IDchem)
        print('Save Testing file: ', recordfile)
        # reopen for the output

        # output info in record files
        fall = open (recordfile+savname,'w+') # record all info

        tmp = 'Testing. Err max {:f} Try inds: {:d}\n'.format(err_max,len(ind))
        if ipvoc: tmp += 'ipvoc: {:d}\n'.format(ipvoc)

        if isinstance(orgPath, str):
            tmp += 'PathTest:{:s}, IDchem_ref:{:s}, IDchem_ref:{:s}, orgPath: {:s}, Namelist:{:s}\n'.format(pathTest, IDchem_ref, IDchem, orgPath, namelist_pre)
        else:
            tmp += 'PathTest:{:s}, IDchem_ref:{:s}, IDchem_ref:{:s}, orgPath:None, Namelist:{:s}\n'.format(pathTest, IDchem_ref, IDchem, namelist_pre)
        fall.write(tmp)

    else: fall = ''

    # compile SSH-aerosol for reduction case
    set_SSH(IDchem,chempath,pathSSH_rdc)    # rdc
    set_SSH(IDchem_ref,refPath,pathSSH_ref) # ref

    # generate the basic namelist for testing
    ssh_namelist, ssh_inds = get_namelist()

    ## run test cases
    # init variables
    nall = 0  # record the total number of test cases
    #nout = 0 # record the number of cases with err > err_max
    emax = 0.0 # record current max error
    nt = len(Tnow)
    eall = np.zeros((ipvoc+1,nt)) # sum up all error to compute the average error
    nsoa = 0 # number of simulations with ave soa <= 0.01

    emax_case = '' # record the case with current max error
    out_ind = [] # save the loc with max err
    
    itrd = 0
    nind = len(ind) # number of conditions

    # path to SSH-aerosol that simulates ref mechanism
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
            # check if need to do simulation
            for j in range(ipvoc+1):
                for il in locs:
                    for inow in Tnow:
                        if j: tmp = '{:s}/{:s}/m{:d}/y{:d}/x{:d}/SOAs_{:d}h.txt.{:d}'.format(pathSSH_ref_abs,orgPath,il[2],il[0],il[1],inow,j)
                        else: tmp = '{:s}/{:s}/m{:d}/y{:d}/x{:d}/SOAs_{:d}h.txt'.format(pathSSH_ref_abs,orgPath,il[2],il[0],il[1],inow)
                        if not os.path.exists(tmp): tag_savorg += 1
                        if j == 0: results_ref.append(tmp)
                        #if tag_savorg: break

        else: tag_savorg = 1 # no previous saved ref

        if tag_savorg: # run ref cases
            res_ref_now, errs = run_SSHsimulation(tail_ref,ResultFolder_ref,pathSSH_ref,
                                  ssh_namelist, ssh_inds,
                                  Ttot,DeltaT,
                                  [locations,locs],
                                  Tnow, pathInit = pathTest,
                                  initfile = ifile,
                                  savorg = orgPath,
                                  ipvoc = ipvoc)
            #results_ref = ['{:s}/{:s}'.format(pathSSH_ref_abs, i) for i in res_ref_now]
        #else:
            #print('use the ref files in {:s}/{:s}'.format(pathSSH_ref,orgPath)) # print to check

        # reduction case
        if sav_soapath: 
            if savname == '_all': isavorg = 'rdc'
            else: isavorg = 'rdc'+savname
        else: isavorg = False
        results_rdc, errs = run_SSHsimulation(tail,ResultFolder,pathSSH_rdc,
                                              ssh_namelist, ssh_inds,
                                              Ttot,DeltaT,[locations,locs],
                                              Tnow,pathInit = pathTest,
                                              initfile = ifile, 
                                              ref_file = results_ref,
                                              savorg = isavorg,
                                              ioutType = 1,
                                              ipvoc = ipvoc)

        ## 3. Analysis results
        n = len(locs) # number of conditions
        errs = [i[0] for i in errs] # only keep ref
        errs = np.array(errs).reshape(ipvoc+1,n,nt) # ipvoc, nloc, nnows
        
        err_now = np.max(errs)# current max one
        # save info if need
        if out_file:
            for i in range(n):
                for k in range(nt):
                    #if errs[i] > err_max: nout += 1
                    c = k+i*nt
                    il = results_rdc[c].replace(ResultFolder+'/','').split('/')[0]
                    info = 'Run:{:d}\t{:s}_{:d}h'.format(nall+c+1,il,Tnow[k])
                    # m6y10x10 -> m6/y10/x10
                    il = '{:s}/{:s}/{:s}'.format(pathSSH_ref_abs,orgPath,il.replace('y','/y').replace('x','/x'))
                    for j in range(ipvoc+1):
                    
                        # get SOA conc
                        if j == 0: sfile = '{:s}/SOAs_{:d}h.txt'.format(il,Tnow[k])
                        else: sfile = '{:s}/SOAs_{:d}h.txt.{:d}'.format(il,Tnow[k],j)

                        if tag_wgas:
                        
                            gconcs = []
                            # get gas concs
                            for x in primaryVOCs:
                                if j == 0: gfile = '{:s}/{:s}_{:d}h.txt'.format(il,x,Tnow[k])
                                else: gfile = '{:s}/{:s}_{:d}h.{:d}.txt'.format(il,x,Tnow[k],j)
                                with open (gfile,'r') as f: infile = f.read().splitlines()
                                gconc = np.zeros(len(infile))
                                for y in range(len(infile)):
                                    gconc[y] = float(infile[y])
                                gconcs.append(gconc)
                                
                            # get soa concs
                            with open (sfile,'r') as f: infile = f.read().splitlines()
                            soas = np.zeros(len(infile))
                            for y in range(len(infile)):
                                soas[y] = float(infile[y])
                            
                            # gas soa yields
                            isoa = []
                            for y in range(1,len(infile)):
                                # get yields at index y
                                dm_voc = 0.0
                                for x in range(len(primaryVOCs)):
                                    k1 = gconcs[x][0] - gconcs[x][y]
                                    if k1 > 0.: dm_voc += k1
                                if dm_voc > 0.: 
                                    isoa.append((soas[y] - soas[0])/dm_voc * 100.)
                            if isoa != []: isoa = sum(isoa)/len(isoa)
                            else: isoa = 0.
                        else:
                        
                            isoa = 0.0
                            # get average soa
                            with open (sfile,'r') as f: infile = f.read().splitlines()
                            for x in infile:
                                isoa += [float(y) for y in x.split(' ') if y !=''][0]
                            if infile != []: isoa /= len(infile)
                            # ave value too low
                            if isoa <= 0.01: 
                                fall.write('Low soa! {:s}\t{:6.4f} werr {:6.4f}\n'.format(sfile,isoa,errs[j][i][k]))
                                nsoa += 1
                        info += '\t{:6.4f} ({:6.4f})'.format(errs[j][i][k], isoa)
                    fall.write(info+'\n')
            info = str(nall)+' Current Max:'
            for j in range(ipvoc+1):  info += '\t{:6.4f}'.format(np.max(errs[j]))
            fall.write(info+'\n')
            fall.flush()

        # record loc and error
        nall += n
        for j in range(ipvoc+1):
            tmp = np.sum(errs[j],axis=0)
            for k in range(nt): eall[j][k] += tmp[k]

        # check
        if err_now > emax: # compare to the current recorded max err
            emax = err_now
            y,x,m = locs[int(np.where(errs == err_now)[1][0])] # gedit ix
            out_ind = [[y,x,m]]
            emax_case = 'm{:d}y{:d}x{:d}'.format(m,y,x)
            
        if err_now > err_out: # compare to the criteria
            break

    # sum
    # average error
    if nall > 0: eall /= nall
    else: raise ValueError('nall <= 0',nall)
    
    err_ave = np.average(eall)
    # max average error
    err_ave_max = 0.0
    tmp = np.average(eall,axis=0) # average ipvocs

    #print('tmp: ',eall, eall.shape,tmp, tmp.shape,ipvoc+1,nt)
    for k in range(nt):
        err_ave_max = max(err_ave_max,tmp[k])
        
    if out_file:
        info = 'Time: \n'
        for k in range(nt):
            info += '{:d}h'.format(Tnow[k])
            for j in range(ipvoc+1): info += '\t{:6.4f}'.format(eall[j][k])
            info += '\n'
            
        fall.write('-----finish-----total run {:d}-----\n Case with max err: {:s}. Ave err {:6.4f}. Ave err max {:6.4f}\n'.format(nall*(ipvoc+1)*2, emax_case, err_ave, err_ave_max)) # nout
        fall.write(info)
        print('Find simulations with low soa:', nsoa)
        fall.write('Find {:d} simulations with low soa concs: ave <= 0.01 ug/m3.\n'.format(nsoa))
        fall.close()
        #print('AT: ', out_ind,emax)

    # mv soa results
    if sav_soapath and isinstance(sav_soapath,str):
        path = '{:s}/{:s}_test'.format(sav_soapath, tail)
        print('Move all testing results to ', path)
        create_folder(path,del_exist = True)
        # save rdc results
        shutil.move(pathSSH_rdc,path)
        # save ref results
        shutil.move(pathSSH_ref+'/'+orgPath,path)

    # print info to check
    if out_ind == []:
        print('Testing result []',out_ind, emax, err_ave, err_ave_max,tail_ref,pathSSH_ref,tail,pathSSH_rdc)
        
    return out_ind, emax, err_ave, err_ave_max

if __name__ == '__main__':
    auto_testing()
