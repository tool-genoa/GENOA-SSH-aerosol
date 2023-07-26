# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  SSHSetUp.py sets up the aerosol box model SSH-aerosol that simulate aerosol
#
#   concentrations for reduction evaluation within the GENOA program.
#
# ================================================================================

import os
import time
import shutil
import numpy as np
from datetime import datetime
from multiprocessing import Pool

from Parameters import namelist_pre, outputs, \
                       pathSSH_sav, ncpu, primaryVOCs
from Functions import is_same_file,replace_in_file, \
                      get_locations, create_folder, \
                      move_results

# changed item in namelist !!! ssh_items and ssh_vals need to be corresponding to each other
ssh_items = ['initial_time =',
              'final_time =',
              'delta_t =',
              'init_gas_conc_file =',
              'init_aero_conc_mass_file =',
              'output_directory =',
              'cst_gas_file =',
              'ref_soa_conc_file =', # for error analysis
              'pre_soa_conc_file = ',
              'latitude =',
              'longitude =',
              'meteo_file =',
              'cst_aero_file =']

def get_ssh_index_and_file(namelist_text, items = ssh_items):
    n = len(items)
    inds = [None] * n # list of index in namelist
    for line in namelist_text:
        for i in range(n): # check if contain
            if items[i] == line[:len(items[i])]:
                if inds[i] == None:
                    inds[i] = namelist_text.index(line) # record ind
                else:
                    print(line, i, items[i], inds[i], namelist_text.index(line))
                    print('SSU: get_ssh_index_and_file: find multiple items')
                    raise ValueError('check namelist and ssh_items.')
                    
    if None in inds:
        print('SSU: get_ssh_index_and_file: namelist not contain all the required items',inds)
        raise ValueError('check namelist and ssh_items.')
    return inds


def run_one_sim(ipath,iloc,itimes, pathInit, rPath,
                sshf, ssh_inds,
                iref_file ='---', ipre_file='---',
                ioutType = 1, isavorg = None, ipvoc=0):

    if len(iloc) == 3: ilc = get_locations([iloc])[0]
    else: ilc = iloc

    iT,iD,inow,imonth = itimes # pass time
    # write new changes of simulation
    ssh_vals=[' {:},\n'.format(imonth + inow * 3600.),
              ' {:},\n'.format(imonth + inow * 3600.+iT),
              ' {:d},\n'.format(iD),
              '"{:s}",\n'.format('{:s}/{:s}/init_gas_{:d}h.dat'.format(pathInit,rPath,inow)),
              '"{:s}",\n'.format('{:s}/{:s}/init_aero_{:d}h.dat'.format(pathInit,rPath,inow)),
              '"{:s}",\n'.format(ipath),
              '"{:s}",\n'.format('{:s}/{:s}/gas.cst'.format(pathInit,rPath)),
              '"{:s}",\n'.format(iref_file),
              '"{:s}",\n'.format(ipre_file),
              '{:6.2f},\n'.format(ilc[0]),
              '{:6.2f},\n'.format(ilc[1]),
              '"{:s}",\n'.format('{:s}/{:s}/meteo.dat'.format(pathInit,rPath)),
              '"{:s}",\n'.format('{:s}/{:s}/aero.cst'.format(pathInit,rPath))]

    # clean & create directory to save results
    if ipvoc: ipath1 = ipath + '.{:d}'.format(ipvoc)
    else: ipath1 = ipath
    
    create_folder(ipath1,del_exist = True)

    # write new namelist
    inamelist = '{:s}/namelist.ssh'.format(ipath1)
    for i,j in enumerate(ssh_inds):
        sshf[j] = ssh_items[i] + ssh_vals[i]
    with open (inamelist, 'w+') as f:
        for i in sshf: f.write(i)

    # launch simulations
    toto = '{:s}/toto'.format(ipath1)
    if ipvoc: 
        os.system('./ssh-aerosol {:s} {:d} 1>{:s}'.format(inamelist, ipvoc, toto))
    else:
        os.system('./ssh-aerosol {:s} 1>{:s}'.format(inamelist, toto))

    # copy organic concs. if need
    if isavorg: 
        os.makedirs(isavorg, exist_ok=True)
        if ipvoc:
            shutil.copy('{:s}/aero/Organics_1.txt'.format(ipath1),
                        '{:s}/SOAs_{:d}h.txt.{:d}'.format(isavorg,inow,ipvoc))
        else:
            shutil.copy('{:s}/aero/Organics_1.txt'.format(ipath1),
                        '{:s}/SOAs_{:d}h.txt'.format(isavorg,inow))
        # wgas
        if '_wgas' in namelist_pre:
            for s in primaryVOCs:
                if ipvoc:
                    shutil.copy('{:s}/gas/{:s}.txt'.format(ipath1,s),
                                '{:s}/{:s}_{:d}h.{:d}.txt'.format(isavorg,s,inow,ipvoc))
                else:
                    shutil.copy('{:s}/gas/{:s}.txt'.format(ipath1,s),
                                '{:s}/{:s}_{:d}h.txt'.format(isavorg,s,inow))

    # analysis errors
    errs = [None, None] # init
    with open (toto) as f:
        for i in f.read().splitlines():
            if 'err_ref' in i and ':' in i: 
                errs[0] = float(i.split(':')[1])
            elif 'err_pre' in i and ':' in i:
                errs[1] = float(i.split(':')[1])
            elif 'sumXk is zero' in i:
                raise ValueError('sumXk is zero is found.',path_ssh,ind,namelist)
    return errs

def run_SSH(IDchem,pathChem,pathSSH,pathResult,tail,ResultFolder,
            Ttot,DeltaT,locs,Tnow,mode,initfile,pathInit,savorg = False,
            ref_file = None, pre_file = None, ioutType = 1, ipvoc=0):
    """build and run ssh-aerosol simulations with given chemistry folder"""

    # clean and compile SSH-aerosol with the new chem
    set_SSH(IDchem,pathChem,pathSSH,mode = mode)

    # write the namelist with the repository of the new chem
    Namelist = '{:s}/namelist_run.ssh'.format(pathSSH)
    ssh_namelist, ssh_inds = get_namelist()

    # Run simulations with settings provided in Parameters.py
    results,errs = run_SSHsimulation(tail,ResultFolder,pathSSH,
                        ssh_namelist, ssh_inds,
                        Ttot,DeltaT,locs,Tnow,pathInit,
                        initfile,savorg,
                        ref_file, pre_file,
                        ioutType,ipvoc)
    if pathResult: # save results in pathResult
        path_old = '{:s}/{:s}'.format(pathSSH,ResultFolder)
        path_new = '{:s}/{:s}'.format(pathResult,ResultFolder)
        move_results(path_old, path_new, items = [i.split('/')[1] for i in results])

    return results,errs

def set_SSH(IDchem,pathChem,pathSSH,pathSSH_sav = pathSSH_sav, mode = 'fast_compile'):
    """
        Objectives:
            modify ssh with the new chemistry
    """
    # check exists
    if not os.path.exists('{:s}/src/ssh-aerosol.f90'.format(pathSSH)):
        shutil.rmtree(pathSSH, ignore_errors=True) # clean
        shutil.copytree(pathSSH_sav[0],pathSSH) # default: fast

    # check mode
    if 'complete' in mode:
        if not is_same_file('{:s}/ssh-aerosol.com'.format(pathSSH_sav[1]), '{:s}/src/ssh-aerosol.f90'.format(pathSSH)):
            print('Running ssh-aerosol-genoa at the complete mode...')
            shutil.copy('{:s}/ssh-aerosol.com'.format(pathSSH_sav[1]),
                        '{:s}/src/ssh-aerosol.f90'.format(pathSSH))
    elif 'fast' in mode:
        if not is_same_file('{:s}/ssh-aerosol.sim'.format(pathSSH_sav[1]), '{:s}/src/ssh-aerosol.f90'.format(pathSSH)):
            print('Running ssh-aerosol-genoa at the fast mode...')
            shutil.copy('{:s}/ssh-aerosol.sim'.format(pathSSH_sav[1]),
                        '{:s}/src/ssh-aerosol.f90'.format(pathSSH))
    else:
        print(mode, 'SSU: mode not exist')
        raise TypeError('Check ssh mode.')

    if 'compile' in mode:
        # update chemistry if name is not BCARY
        if IDchem != 'BCARY':
            #tag = 0
            for i in ['reactions', 'species', 'aer.vec', 'mol', 'RO2']:
                tmp0 = '{0:s}/{1:s}/{1:s}.{2:s}'.format(pathChem,IDchem,i) #new
                tmp1 = '{0:s}/src/include/CHEMISTRY/BCARY/BCARY.{1:s}'.format(pathSSH,i) #old
                if not is_same_file(tmp0, tmp1): shutil.copy(tmp0, tmp1)
                    #if i in ['reactions', 'species']: tag += 1
            #if tag: # remove pervious executive files
            #    os.system('{0:s}/src/include/CHEMISTRY/BCARY/*.f90'.format(pathSSH))
        cmd = './compile >tmptoto  2>&1'
        if 'clean' in mode: 
            print('Complete clean SSH before compile...')
            cmd = './clean >tmptoto;' + cmd

        # go to SSH-aerosl path and compile
        wd = os.getcwd()
        os.chdir(pathSSH)
        os.system(cmd)

        # recompile if need
        with open ('tmptoto','r') as f: totofile = f.read()    
        if 'scons: done building targets.' in totofile: tag = 0
        else:
            tag = 1
            # if error, run second time
            if 'error' in totofile:
                print('Has error: run again...')
                os.system(cmd)

                # check message of the second time running
                with open ('tmptoto','r') as f: totofile = f.read()
                if 'scons: done building targets.' in totofile:
                    tag = 0
                else:
                    # check the third time
                    print('Run the third time...')
                    os.system('./clean >tmptoto  2>&1; ./compile >tmptoto 2>&1')
                    with open ('tmptoto','r') as f: totofile = f.read()
                    if 'scons: done building targets.' in totofile:
                        tag = 0
                    else:# need to check manually
                        print('Still has error in compile.')
                if tag:
                    print('SSH: find error in the compiling of ssh-aerosol.', pathSSH)
                    raise RuntimeError('Check the tmptoto file in ssh-aerosol.')
            else:
                print('SSH: stopped and not found error in the compile message im tmptoto', pathSSH)
                raise RuntimeError('Check ssh-aerosol compile.')
        # return to the pervious path
        os.chdir(wd)

def get_namelist(pathSSH_sav='../', outConc = True):
    """
        Objectives:
            rewrite namelist file for running all simulations
    """

    with open ('{:s}/src/{:s}'.format(pathSSH_sav,namelist_pre),'r') as f: ssh_namelist = f.readlines()

    # write setting for outputs
    if outConc: ioutType = 1
    else: ioutType = 0

    items = ['output_type = '] # output
    news =  ['{:d},\n'.format(ioutType)]

    if outputs['aero'] != []: #out_aero = '"---"' # no output
        out_aero = '' # write output
        for i in outputs['aero']:
            out_aero += '"{:s}",'.format(i)
        n_out_aero = len(outputs['aero'])
        items += ['n_output_aero =', 'output_aero =']
        news += [ '{:d},\n'.format(n_out_aero), '{:s}\n'.format(out_aero)]

    if outputs['gas'] != []: #out_gas = '"---"' # no output
        out_gas = '' # write output
        for i in outputs['gas']:
            out_gas += '"{:s}",'.format(i)
        n_out_gas = len(outputs['gas'])
        items += ['n_output_gas =', 'output_gas =']
        news += ['{:d},\n'.format(n_out_gas), '{:s}\n'.format(out_gas)]

    # change namelist
    if items != []:
        inds = get_ssh_index_and_file(ssh_namelist, items)
        for i,j in enumerate(inds):
            ssh_namelist[j] = items[i] + news[i]

    # get items that need to be changed in simulations
    ssh_inds = get_ssh_index_and_file(ssh_namelist)

    return ssh_namelist, ssh_inds

def run_SSHsimulation(tail,ResultFolder,pathSSH,
                        ssh_namelist, ssh_inds,
                        Ttot,DeltaT,locs,Tnow,pathInit,
                        initfile = 'storage',savorg = False,
                        ref_file = None, pre_file = None,
                        ioutType = 1, ipvoc=0):
    """run simulations
       Ttot: an array of the simulation duration
       DeltaT: simulation time step
       gas: path of input inital gas-phase distribution
       aer: path of input initial aerosol phase distribution
       Tnow: marker of <aer> and the simulation start hour at the start date
    """

    if savorg and not isinstance(savorg,str): 
        print('SSU: if provides, savorg should be a string.', savorg)
        raise ValueError('check savorg')

    # for testing with generating init files
    if '_' in initfile:
        initfile,inittail = initfile.split('_')

    # prepare simulations: build lists
    results = [] # save all paths
    sim_times = [] # in the format iT,iD,inow,imonth
    rPaths = [] # storage file reminder
    isavorgs = []
    ilocs = []

    if len(locs[0]) == len(locs[1]) and len(locs[0][0]) == 2:
        ilocations = locs[0]
    else:
        ilocations = locs[1] # get location later
          

    for il in range(len(locs[1])):

        # location index in the map
        y,x,m = locs[1][il]
        
        # if sav org, the path
        if savorg: orgpath = '{:s}/m{:d}/y{:d}/x{:d}/'.format(savorg,m,y,x)
        else: orgpath = False
        isavorgs.append(orgpath)
        
        # file names
        if initfile == 'storage':
            lc = 'm{:d}y{:d}x{:d}'.format(m,y,x)
            rPath = 'm{:d}/y{:d}/x{:d}'.format(m,y,x)
        else:
            if initfile == 'month': lc = 'm{:d}y{:d}x{:d}'.format(m,y,x)
            elif initfile == 'Test': lc = '{:s}_{:d}'.format(inittail,il)
            else: lc = 'y{:d}x{:d}'.format(y,x)
            rPath = lc
        rPaths.append(rPath)

        # time of the month
        imonth = (datetime(2015,m+1,1,0) - datetime(2015,1,1,0)).total_seconds()
        for inow in Tnow: # for different time
            for idx in range(len(Ttot)): # for different time settings

                # time settings
                iT=Ttot[idx] #total time
                iD=DeltaT[idx] #time step
                sim_times.append([iT,iD,inow,imonth])

                # save path
                ilocs.append(il)
                ipath = '{:s}/{:s}/c_{:d}_{:d}_{:d}h_{:s}'.format(ResultFolder,lc,iT,iD,inow,tail)
                results.append(ipath)

    # prepare error files
    n = len(results)
    if ref_file == None: iref_files = ['---']*n
    elif isinstance(ref_file, list):
        if len(ref_file) != n:
            print(len(ref_file), n, ref_file)
            raise ValueError('SSU: ref_file length not fit.')
        elif '.txt' in ref_file[0]:
            iref_files = ref_file
        else:
            iref_files = ['{:s}/aero/Organics_1.txt'.format(ref_file[i]) for i in range(n)]
    else:
        print(ref_file)
        raise TypeError('SSU: not recognize ref_file type.')

    if pre_file == None: ipre_files = ['---']*n
    elif isinstance(pre_file, list):
        if len(pre_file) != n:
            print(len(pre_file), n, pre_file)
            raise ValueError('SSU: pre_file length not fit.')
        elif '.txt' in pre_file[0]:
            ipre_files = pre_file
        else:
            ipre_files = ['{:s}/aero/Organics_1.txt'.format(pre_file[i]) for i in range(n)]
    else:
        print(pre_file)
        raise TypeError('SSU: not recognize pre_file type.')

    # go to the SSH path
    wd = os.getcwd()
    os.chdir(pathSSH)

    pool_inputs = []
    js = [i for i in range(ipvoc+1)]
    for j in js:
        for i in range(n):
            pool_inputs.append((results[i], ilocations[ilocs[i]],
                             sim_times[i], #ipath,iloc,itimes
                             pathInit, rPaths[ilocs[i]], # pathInit, rPath
                             ssh_namelist, ssh_inds, # sshf, ssh_inds
                             iref_files[i], ipre_files[i], # iref_file, ipre_file
                             ioutType, isavorgs[ilocs[i]],j))
    with Pool(ncpu) as pool:
        errs = pool.starmap(run_one_sim, pool_inputs)

    # return to the pervious path
    os.chdir(wd)
    
    return results, errs
