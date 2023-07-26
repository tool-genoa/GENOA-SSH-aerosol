# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  AutoTrainingParallel.py runs the parallel reduction. 
#
# ================================================================================

import os
import time
import json
import shutil
import numpy as np

from copy import deepcopy
from datetime import datetime
from itertools import combinations
from multiprocessing import Pool

from Parameters import RunSets, Roptions, prefix, \
                       pathSSH_rdc, locs, Ttot, DeltaT, \
                       Tnow, pathInitFiles, pathNewChem, \
                       initfile, pathNewRes, pathSSH_sav, \
                       NokeepSp, primaryVOCs, tail, \
                       ncpu, ipvocs, nviz
from Functions import  isfloat, isint, get_info, \
                       create_folder, get_locations, \
                       count_cut, combinations_index, \
                       combinations_number, round_up, \
                       read_parameters_from_csv
from KineticMCMtoPython import update_kinetic_rate
from DataStream import read_chem_sets,to_SSH_sets
from ReductionStrategy import reduce_reaction_byLump, \
                              reduce_reaction_bySpecies, \
                              trim_scheme, reaction_seperate, \
                              reaction_merge
from SSHSetUp import run_SSH, get_namelist, set_SSH, ssh_items
from AutoTesting import auto_testing

# primary voc types
#ipvocs = 3 # 1 for only apinene, 2 for only bpinene, 3 for only limonene

def auto_training_para(Setups = RunSets[3], locs = locs, pathInitFiles = pathInitFiles, ATSetups = RunSets[4]):
    """
        Objectives: run GENOA reduction

        Inputs:
            
        Outputs:
    """
    global refconc_paths, concs_ave, RefConcRead, nlimtar, nlimcmb

    ### initialization
    tbeg = time.perf_counter() # starting time

    # computation of the lumping ratios
    lptypes = ['plain'] # options: 'plain','tau','0','100','25','50','75']

    # reuse tag: active - 2: reuse ref_conc & train_order 
    tag_reuse = Setups['Reuse_chem']
    
    # read from Parameters.py
    IDchemPre = Setups['IDchemPre']
    prePath = os.path.abspath(Setups['prePath'])
    pre_chem_path = Setups['preChemPath']

    IDchemRef = Setups['IDchemRef']
    refPath = os.path.abspath(Setups['refPath'])

    IDchemFake = Setups['IDchemFake']
    fakePath = os.path.abspath(Setups['fakePath'])

    # freeze compounds
    frozen_species = Setups['frozenspecies']
    kept_species = Setups['keptspecies']
    if frozen_species != []:
        print('Read frozen species that do not participate reduction:', frozen_species)
    if kept_species != []:
        print('Read kept species that has to be kept in the reduced mechanism:', kept_species)

    # settings for pre-testing
    tag_pre_testing0 = 0 # tag at the last cycle - default 0
    err_pre_testing = float(Setups['err_pre_testing'])
    err_stage_change = float(Setups['err_stage_change'])
    
    try_ave_ref = float(Setups['try_ave_ref']) # allowed ave err
    try_max_ref = float(Setups['try_max_ref']) # allowed max err
    err_limit = float(Setups['err_limit']) # limitation of error tolerances
    
    # reduction strategies
    strategy_types = Setups['strategy_types']
    # available strategies
    strategy_types_all = ['rm','rm1','rmp','jp','lp','rp','rs','rsn','da']
    for s in strategy_types:
        if s not in strategy_types_all:
            raise ValueError('Not found the required reduction strategy. ',s, strategy_types_all)
    
    # assignment from table
    input_file = Setups['training_parameter_table']
    # name in the input table: name of variable in the program
    err_ref = err_limit, 
    input_vars = { # error tolerances
                   'err_ref': 'err_ref', 'err_rav': 'err_rav',
                   'err_pre': 'err_pre', 'err_pav': 'err_pav',
                   'delta_err': 'delta_err',
                   # tags
                   'tag_pre_testing': 'tag_pre_testing',
                   'tag_stage': 'tag_stage',
                   'tag_rm1': 'tag_rm1', 
                   'tag_efficient': 'tag_efficient'  
                 }
    # read values
    if isinstance(input_file,str):
        print('Reading reduction parameters from table: ', input_file)
        input_table = read_parameters_from_csv(input_file)
        npara = len(input_table[list(input_table.keys())[0]]) # len of inputs
    else: # init
        input_table = {}
        npara = len(Setups['err_ref']) # init lists with the length of err_ref input
    
    input_params = {} # init parameters & options
    for key in input_vars.keys():
        # local name
        val = input_vars[key]
        # assign values
        # read from table
        if key in input_table.keys():
            input_params[val] = input_table[key]
            print(f'Read {key} for variable {val} from table: ', input_params[val])
            
        # read from ini file + Parameter.py
        elif val in Setups.keys():
            if isinstance(Setups[val],list):
                # check lenth
                if len(Setups[val]) != npara:
                    raise ValueError (f'length of list {val} is not consistent. Check config file. ',n,len(Setups[key]))
            # assign a list
            input_params[val] = [Setups[val]] * npara
            print(f'Build variable {val}: ', input_params[val])
            
        # No default
        else:
            raise NameError('key not found in table, config file, not Parameters.py', key, val)

    # if record reduction tree
    ntree_ctl = Setups['ntree']
    if ntree_ctl <= 1: tag_tree = False
    # save reduction tree nodes: apt_cas
    # save chems and records in reduction tree folder: ntree/iauto/ntry
    # tag_all = False # record all reduction trees

    # extract pre-testing locations from testing file
    nPreTest = Setups['nPreTest']
    
    tag_reach = Setups['tag_reach'] # 2 search, 1 active, 0 inactive, -1 record (2->-1)
    err_rav_reach = Setups['err_rav_reach']

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
            raise OSError('Not provied Test_file.')
    else:
        raise ValueError('nPreTest can not be read.', nPreTest)
        
    # time limitation: the next reduction cycle won't start if time exceeds the limitation
    tlim = Setups['time_limit'] # unit in second
    if tlim: print('Set time limitation: ',tlim)
    else: print('No time limitiation for training.')
    tag_tlim = False
        
    ### check/assgin number of processors used per strategy
    if ncpu < 16:
        print('Warning! Reduction can be very slow with this number of processors: ',ncpu)
        raise ValueError('ncpu should >= 16 for parallel reduction.')
    
    # read limitation for number of candidate reductions via removing reactions
    nlimtar = Setups['nrcn_limit']
    
    # assgin processros for lp and rm: rm - one processor * 6 reactions
    ncpus = {}
    s,n = 0,0 # num for other strategies
    for l in strategy_types:
        if l in ['rs','da','rsn']: 
            s += 1 # number of processors used
            ncpus[l] = 1
        elif l in ['jp','rp']: 
            s += 2
            ncpus[l] = 2
        elif l in ['rm','rm1','rmp']: 
            n += 1
        elif l == 'lp': 
            n += 3 # 3 times more processors

    # check numbers in case 
    # removing    
    tmp = list(set(['rm','rm1','rmp']) & set(strategy_types))
    if tmp != []: 
        ncpurm = int((ncpu-s)/n) # max number of processors for removing
        if nlimtar <= 0:
            # max try of combination times
            for i in range(10): # try maximum 10 times
                if combinations_number(i+1) > ncpurm * 6: # 6 rdc per processor
                    # max num of tar
                    nlimtar = i
                    # mas number of reduction attempts for removing at one time
                    nlimcmb = combinations_number(i)
                    break
        else: nlimcmb = combinations_number(nlimtar)
                
        if ncpurm >= 1:
            print('For removing reactions. Max number of processors used: ',ncpurm,' Max number of n',nlimtar, nlimcmb,'\n')
        else: 
            raise ValueError('Max number of processors for removing reactions is less than 1.',ncpurm,nlimtar)
        for l in tmp: ncpus[l] = ncpurm
    else: # no removing
        ncpurm = 0
        nlimtar = 0
    
    # lumping - use all rest cpu
    if 'lp' in strategy_types: ncpus['lp'] = max(0,ncpu - sum(ncpus.values()))
    print('Processors usage distribution: ',ncpus)

    # settings for parallel reductions
    train_key = Setups['train_key']
    if train_key not in ['reaction','species']:
        raise ValueError('train_key should be reaction or species. ', Setups['train_key'])

    # reduction order file
    if isinstance(Setups['train_order'],str):
        iforder = Setups['train_order']
        print('Read train_order list',iforder)
    else:
        iforder = None # init generate it later
        
        # serach files if reuse
        if tag_reuse == 2:
            # get chem number
            i = int(IDchemPre.split('-')[-1])
            # get chems forlder name
            s = '/'+IDchemPre.split('-')[0].replace(prefix,'')
            if s+'_chems' in pre_chem_path:
                l = '{:s}/reduction.use.{:d}'.format(s.replace(s+'_chems',s+'_recs'),i+1)
                if os.path.exists(l):
                    iforder = l
                    print('Reuse train_order list: ', iforder)

    ### prepare reduction - init parameters that can be updated by restart file
    # other reduction options - currently not read from Parameters.py
    tag_redo = 0 # tag_redo == 1: redo the current reduction cycle
    # output message
    tag_record = 0 # tag_record == 1: output detailed reduction infos

    iauto = 0 # initial reduced times - used to build mechanism name
    ntree = 0 # index of the reduction tree - not used if tag_tree == 0

    # errors
    nerr = 0 # index of error tolerances
    ifit_pre, ierr_pre = 0.0, 0.0 # ave errs of the previous case

    # load/update parameters with refile - for strat from a break-point
    refile = Setups['restart']
    if refile:
        if os.path.exists(refile):
            print('Reading restart file: ',refile)
            with open(refile, 'r') as f: re_data = json.loads(f.read())
            for i in ['iauto','ntree','nerr','ifit_pre','ierr_pre']:
                if i in re_data: # find changes
                    print('- Update ',i,': ',re_data[i])
                    exec('{0:s}=re_data["{0:s}"]'.format(i))
        else:
            raise FileNotFoundError('ATP: restart file is not found.',refile)

    # if ref case is not the pre case. run one time the pre case to ensure the accuracy of pre results concs
    if IDchemPre == IDchemRef: 
        tag_tchemP = False
    else:
        print('Input IDchemPre! = IDchemRef. Active tag_tchemP.', IDchemPre, IDchemRef) 
        tag_tchemP = True # run pre mechanism to check errs between pre and ref mechanisms

    ### Directories
    # absolute path
    pathSSH_rdc_abs = os.path.abspath(pathSSH_rdc)
    # sav chemistry/record files in pathNewChem
    path_sav_rec = '{0:s}/{1:s}_recs'.format(pathNewChem,tail)
    path_sav_chem = '{0:s}/{1:s}_chems'.format(pathNewChem,tail)
    path_sav_res = '{0:s}/{1:s}'.format(pathNewRes,tail)

    # sav soa ref/pre files in pathSSH_rdc_abs
    pathRef = '{:s}/Ref_{:s}'.format(pathSSH_rdc_abs,IDchemRef)
    pathPre = '{:s}/Pre_now'.format(pathSSH_rdc_abs)
    pathPre_sav = '{:s}/Pres'.format(pathSSH_rdc_abs)
    pathTree = '{:s}/tree'.format(pathSSH_rdc_abs)

    # save all changes
    path_ssh_chem = '{:s}/chems'.format(pathSSH_rdc_abs)

    # sav namelist
    pathNml = '{:s}/nmls'.format(pathSSH_rdc_abs)

    # create folder if need
    for i in path_sav_rec,path_sav_chem, \
             pathRef,path_sav_res,pathTree, \
             pathPre,pathPre_sav,pathNml, \
             path_ssh_chem:
        create_folder(i)

    # cp pre chem
    if os.path.exists('{:s}/{:s}'.format(pre_chem_path,IDchemPre)):
        os.system('cp -rf {:s}/{:s} {:s}/'.format(pre_chem_path,IDchemPre,path_sav_chem))
    else:
        raise FileNotFoundError('Can not find IDchemPre: {:s}/{:s}'.format(pre_chem_path,IDchemPre))

    # add toto folder in ssh/
    if os.path.exists('{:s}/ssh'.format(pathSSH_rdc_abs)):
        # save toto files
        os.system('mkdir -p {:s}/ssh/toto'.format(pathSSH_rdc_abs))
    else:
        raise FileNotFoundError('Can not find ssh folder:', '{:s}/ssh'.format(pathSSH_rdc_abs))

    # generate record files
    if iauto > 0: frec = open ('{:s}/CaseRdc_{:s}'.format(pathNewChem,tail),'a+') # add
    else: frec = open ('{:s}/CaseRdc_{:s}'.format(pathNewChem,tail),'w+') # renew

    # get coordinate of conditions
    nlocs = len(locs) # number of training dataset
    locations = get_locations(locs + prelocs)
    
    # get initconcs
    # regenerate ref concs
    if tag_reuse == 2: RefConcReadin = None
    else: RefConcReadin = path_sav_chem+'/s_'+tail
        
    concs_ave, concs_min, refconc_paths, RefConcRead = update_kinetic_rate(locs,
                                                          IDchem_fake = IDchemFake,
                                                          path_fake = fakePath,
                                                          pathInitFiles = pathInitFiles,
                                                          RefConcRead = RefConcReadin)
    # read from files
    if tag_reuse == 2: RefConcRead = path_sav_chem+'/s_'+tail + str(len(locs))

    ## build general namelist & ref/pre files
    ssh_namelist, ssh_inds = get_namelist()
    to_pre_paths = [] # copy results
    
    # get info for namelist
    ntnml,npnml, n = 0, 0, 0 # number of namelist
    tnmls,trslts,tcases = [],[],[] # for parallel of training dataset
    pnmls,prslts,pcases = [],[],[] # for parallel of pre-testing dataset

    for il,iloc in enumerate(locs + prelocs):
        # location index in the map
        y,x,m = iloc
        ilc = locations[il]
        if initfile == 'storage':
            lc = 'm{:d}y{:d}x{:d}'.format(m,y,x)
            rPath = 'm{:d}/y{:d}/x{:d}'.format(m,y,x)
        else:
            if initfile == 'month': lc = 'm{:d}y{:d}x{:d}'.format(m,y,x)
            elif initfile == 'Test': lc = '{:s}_{:d}'.format(inittail,il)
            else: lc = 'y{:d}x{:d}'.format(y,x)
            rPath = lc

        # time of the month
        imonth = (datetime(2015,m+1,1,0) - datetime(2015,1,1,0)).total_seconds()
        for inow in Tnow: # for different time
            for idx in range(len(Ttot)):
                # time settings
                iT=Ttot[idx] #total time
                iD=DeltaT[idx] #time step

                # save path
                icase = '{:s}/c_{:d}_{:d}_{:d}h'.format(lc,iT,iD,inow)
                ipath = 'Results/{:s}'.format(icase) # relative path
                ifile = '{:s}/SOAs_{:d}h.txt'.format(rPath,inow)
                iref = '{:s}/{:s}'.format(pathRef,ifile)

                if il < nlocs: ipre = '{:s}/{:s}'.format(pathPre,ifile)
                else: ipre = '---'

                # write new changes of simulation
                ssh_vals=[' {:},\n'.format(imonth + inow * 3600.),
                          ' {:},\n'.format(imonth + inow * 3600.+iT),
                          ' {:d},\n'.format(iD),
                          '"{:s}",\n'.format('{:s}/{:s}/init_gas_{:d}h.dat'.format(pathInitFiles,rPath,inow)),
                          '"{:s}",\n'.format('{:s}/{:s}/init_aero_{:d}h.dat'.format(pathInitFiles,rPath,inow)),
                          '"{:s}",\n'.format(ipath),
                          '"{:s}",\n'.format('{:s}/{:s}/gas.cst'.format(pathInitFiles,rPath)),
                          '"{:s}",\n'.format(iref),
                          '"{:s}",\n'.format(ipre),
                          '{:6.2f},\n'.format(ilc[0]),
                          '{:6.2f},\n'.format(ilc[1]),
                          '"{:s}",\n'.format('{:s}/{:s}/meteo.dat'.format(pathInitFiles,rPath)),
                          '"{:s}",\n'.format('{:s}/{:s}/aero.cst'.format(pathInitFiles,rPath))]

                # write new namelist
                inamelist = '{:s}/namelist_{:d}.ssh'.format(pathNml,n)
                for i,j in enumerate(ssh_inds):
                    ssh_namelist[j] = ssh_items[i] + ssh_vals[i]
                with open (inamelist, 'w+') as f:
                    for i in ssh_namelist: f.write(i)

                # build folder
                if idx == 0:
                    for i in pathRef, pathPre:
                        os.makedirs('{0:s}/{1:s}'.format(i,rPath), exist_ok=True)
                # copy ref files
                iorg = '{:s}/{:s}_{:s}'. \
                        format(refPath,icase,IDchemRef.replace(prefix,''))
                os.system('cp {:s}/aero/Organics_1.txt {:s}'.format(iorg,iref))
                for i in range(ipvocs):
                    os.system('cp {0:s}.{2:d}/aero/Organics_1.txt {1:s}.{2:d}'.format(iorg,iref,i+1))
                if il < nlocs:
                    # pre files only for locs
                    iorg = '{:s}/{:s}_{:s}'. \
                            format(prePath,icase,IDchemPre.replace(prefix,''))
                    os.system('cp {:s}/aero/Organics_1.txt {:s}'.format(iorg,ipre))
                    for i in range(ipvocs):
                        os.system('cp {0:s}.{2:d}/aero/Organics_1.txt {1:s}.{2:d}'.format(iorg,ipre,i+1))
                    to_pre_paths.append(['{:s}'.format(ipath),'{:s}/{:s}'.format(pathPre,ifile)])
                    tnmls.append(inamelist)
                    trslts.append(ipath)
                    tcases.append(icase)
                    ntnml += 1
                else:
                    pnmls.append(inamelist)
                    prslts.append(ipath)
                    pcases.append(icase)
                    npnml += 1
                n += 1

    # init training conditions
    nnml,nmls,rslts,cases = ntnml,tnmls,trslts,tcases 

    # backup pathPre
    os.system('cp -rf {:s} {:s}/{:s}'.format(pathPre, pathPre_sav, IDchemPre))

    # init parameters for reduction, not change during restart
    IDchem = None # chem id init
    tag_iauto = 1 # if 0 stop reduction
    emax_ind, emax, eave, eave_max = [], 0.0, 0.0, 0.0 # init for testing
    
    ## reduction start -- reduction cycle
    while tag_iauto:

        t0 = time.perf_counter()

        if IDchem and nval != 0 and not tag_redo: # update pre
            IDchemPre = IDchem
            prePath = curPath

        # init pre reaction/species list
        speciesfile = '{:s}/{:s}/{:s}.mol'.format(path_sav_chem,IDchemPre,IDchemPre)
        reactionfile = '{:s}/{:s}/{:s}.reactions'.format(path_sav_chem,IDchemPre,IDchemPre)

        # name to cycle
        if tag_tree: IDchem = '{:s}-{:d}-{:d}'.format(prefix+tail,ntree,iauto) 
        else: IDchem = '{:s}-{:d}'.format(prefix+tail,iauto) 

        print(iauto, 'load IDchemPre & prePath: ',IDchemPre,prePath, '. current ID: ',IDchem)

        # assign parameters
        # get index
        if nerr < npara: i = nerr
        else: i = -1
        
        #for key,value in input_params.items():
        #    exec('{:s} = {:}'.format(key,value[i]))
        err_ref = input_params['err_ref'][i]
        err_rav = input_params['err_rav'][i]
        err_pre = input_params['err_pre'][i]
        err_pav = input_params['err_pav'][i]
        delta_err = input_params['delta_err'][i]
        tag_pre_testing = input_params['tag_pre_testing'][i]
        tag_stage= input_params['tag_stage'][i]
        istage = 0
        tag_rm1= input_params['tag_rm1'][i]        
        tag_efficient= input_params['tag_efficient'][i]
         
        # condition changes
        if tag_pre_testing == 2: # need to change to 0 or 1 at this step
            # compare eave to set conditions
            if eave == 0: tag_pre_testing = tag_pre_testing0
            else:
                if eave >= err_pre_testing: 
                    tag_pre_testing = 1
                    # set all tags after this step
                    for i in range(len(input_params['tag_pre_testing'][nerr:])):
                        if input_params['tag_pre_testing'][nerr+i] == 2:
                            input_params['tag_pre_testing'][nerr+i] = 1
                        elif input_params['tag_pre_testing'][nerr+i] == 0:
                            raise ValueError('tag_pre_testing should not change from 1 to 0', nerr, input_params['tag_pre_testing'][nerr:])
                else:
                    tag_pre_testing = 0
                    
        # switch conditions
        if tag_pre_testing != tag_pre_testing0:
            # redo tchem to get ifit_pre & ierr_pre
            if not tag_tchemP:
                tag_tchemP = True 
                ifit_pre, ierr_pre = 0.,0.
            if tag_pre_testing == 1:
                print('Active reduction under pre-testing conditions.',IDchem)
                # set pre-testing datasets
                nnml,nmls,rslts,cases = npnml,pnmls,prslts,pcases
            elif tag_pre_testing == 0:
                print('Active reduction under training conditions.',IDchem)
                # set training datasets
                nnml,nmls,rslts,cases = ntnml,tnmls,trslts,tcases
        # save tag
        tag_pre_testing0 = tag_pre_testing

        # check error tolerance        
        if tag_pre_testing:
            # need to have inputs for 
            if err_ref > err_limit or err_rav > err_limit:
                raise ValueError(f'err_ref {err_ref} or err_rav {err_rav} given is larger than err_limit {err_limit}.')
            # errors not used for pre-testing
            err_pre, err_pav = err_limit, err_limit
            # read delta_err as a ratio of err_rav
            if delta_err > err_limit: 
                delta_err = err_rav/delta_err

            # active tag_reach
            if tag_reach:
                if err_rav != err_rav_reach and round(err_rav,2) == round(err_rav_reach,2):
                    tag_reach = 1
                    print('Active tag reach! err_rav: {:f} -> tyr_ave_ref: {:f}\n'.format(err_rav, err_rav_reach))

        else: # training condition
            # need to have input err_pre and err_ref < limit
            if err_pre > err_limit or err_ref > err_limit:
                raise ValueError(f'err_pre {err_pre} or err_ref {err_ref} given is larger than err_limit {err_limit}.')
            # read delta_err as a ratio of err_ref
            if delta_err > err_limit: delta_err = err_ref/delta_err
            if err_pav > err_limit: err_pav = err_pre/err_pav
            
        # prepare the record files
        recordfile = '{:s}/Rec_{:s}'.format(path_sav_rec,IDchem)
        fall = open (recordfile+'_all','w+')
        fuse = open (recordfile+'_use','w+')

        # record output info
        tmp_out = 'preIDchem: {:s}\tprePath: {:s}\nrefIDchem: {:s}\trefPath: {:s}\n'.format(IDchemPre,prePath,IDchemRef,refPath)

        if tag_pre_testing: 
            tmp_out += 'Under pre-testing dataset!\n'
        else:
            tmp_out += 'Under training dataset!\n'

        # error tolerance
        tmp_out += 'Error Tolerances: err_ave_ref <= {:f}, err_max_ref <= {:f}, err_ave_pre <= {:f}, err_max_pre <= {:f}, delta_err <= {:f}\n'.format(err_rav,err_ref,err_pav,err_pre,delta_err)

        # reduction options
        tmp_out += 'Options:\t'
        
        if tag_stage: 
            tmp_out += 'aerosol-oriented treatment - {:d};\t'.format(tag_stage)
            info1 = '' # init

        if tag_redo: 
            tmp_out += 'tag_redo;\t'
            tag_redo = 0 # reset
        
        if tag_efficient:
            tmp_out += 'efficient treatment <= {:d};\t'.format(tag_efficient)

        if tag_reach == 1:
            tmp_out += 'with tag_reach to errs {:f}->{:f};\t'.format(err_rav,err_rav_reach)

        # reduction tree ntree, record the root case of this sub branch
        if tag_tree and ntree: 
            tmp_out += 'ntree: {:d} & Pre-case: {:};\t'.format(ntree, save_data)

        for f in [fall,fuse,frec]:
            f.write(tmp_out+'\n')
            f.flush()

        # prepare Path
        ResultFolder = 'Results_{:s}'.format(IDchem)
        curPath = '{:s}/{:s}'.format(os.path.abspath(path_sav_res),ResultFolder)

        # read intial reaction/species list
        rc,sp=read_chem_sets(reactionfile,speciesfile,speciesType='SSH',reactionType='SSH')
        sps = [i.name for i in sp]

        # check frozen and kept species
        if frozen_species != []:
            tmp = []
            for s in frozen_species:
                if s not in sps: tmp.append(s)
            if tmp != []:
                raise NameError('Frozen species not found in the scheme:', speciesfile, tmp,frozen_species)
                
        if kept_species != []:
            tmp = []
            for s in kept_species:
                if s not in sps: tmp.append(s)
            if tmp != []:
                raise NameError('Kept species not found in the scheme:', speciesfile, tmp,kept_species)

        # active elementory-like treatment if nrcn changes
        if tag_rm1 == 1:

            i = len(rc)
            rc = reaction_seperate(rc,sp)
            
            if len(rc) != i:
                tag_tchemP = True # test in the next round
                print('Active elementory-like treatment. No. of reaction before and after: ', i,len(rc))
                frec.write('=+=+ rm1. No. of reaction changes from {:d} to {:d}.\n'.format(i,len(rc)))
            else:
                frec.write('=+=+ try rm1. Number of reactions {:d} is not changed. \n'.format(i))

        # record
        nrea, ngas, naer = get_info(rc,sp,'')
        frec.write('Initial scheme No.reaction: {:d}\tNo.gas: {:d}\tNo.aerosol: {:d}\n'.format(nrea, ngas, naer))

        # check if active tag to record reduction trees
        if ntree_ctl and not tag_tree:
            # active pretesting if
            if tag_pre_testing and (nerr == npara - 1 or naer <= 50):
                tag_tree = ntree_ctl
                next_tree = [] # init # record tree process
                frec.write('Record reduction subtrees !!!\n')

        # init valid attempt times
        nval,ntry = 0,0 # number of valid reductions and of tried reduction attempts
        frosps, frorcn = [],[] # storage species already tested
        npres = {} # record the valid number in the previous round
        for i in strategy_types: npres[i] = 0
        
        # items for searching for reduction
        a_item = []
        if train_key == 'reaction':
            pass
        elif train_key == 'species':
            if iforder: # read
                print('Read the reduction order from file: ',iforder,iauto)
                with open (iforder,'r') as f:
                    # tmp file, could be partial file
                    if 'rdc_order.tmp' in iforder:
                        a_item = json.loads(f.read())
                        print('Read reduction list from json file. Find number of species',len(a_item))
                    else:
                        tmp = f.read().splitlines()[1:] # remove the first comment line
                        for i in tmp:
                            if i in sps: 
                                if i not in a_item: a_item.append(i)
                            else: raise NameError('Species in order list not found in the current species list: ',i)
                        print('Read reduction list from file (Not read 1st line!):  Find number of species ',len(a_item),len(tmp))

                if len(a_item) == 0: 
                    raise ValueError('ATP: read len of a_item == 0',iforder)
                else:
                    print('Read ',len(a_item),' species from ',iforder)
                    frec.write('rdc order is read: {:}\t{:d}\n'.format(iforder,len(a_item)))
                    # update for next
                    if len(a_item) != len(sps): # a_item is a partial file
                        iforder = None

            else: # generate order 
                print('Generate the reduction order.',iauto)
                # save chems
                for i,j in enumerate(sps):
                    if not sp[sps.index(j)].organic: continue
                    if j in frozen_species: continue
                    if j in primaryVOCs: continue
                    if j not in a_item:
                        a_item.append(j) # save name

                # get paths
                chems = []
                # prepare for further reduction
                for l in ['rs','da','rsn']: # not use of cut
                    if l in strategy_types:
                        i = 1
                        for j in range(i):
                            chems.append('{:s}-{:d}'.format(l,j))
                i = [l for l in ['jp','rp','rm','rm1','rmp','lp'] if l in strategy_types]
                n1 = len(a_item)
                n2 = int(n1/sum([ncpus[l] for l in i])) + 1
                for j in range(n2): # max number
                    for l in i: # strategy
                        for k in range(ncpus[l]): # processor
                            chems.append('{:s}-{:d}_{:d}'.format(l,k+1,j))
                            if len(chems) >= n1: break
                        if len(chems) >= n1: break
                    if len(chems) >= n1: break
                            
                print('Generate chems to build reduction order.',n1,len(chems))
                if len(chems) < n1: raise ValueError('len a_item should not be smaller than len chems.',n1,len(chems))
                # get chems
                pool_inputs = [(rc,sp,j,chems[i],path_ssh_chem) for i,j in enumerate(a_item)]
                with Pool(ncpu) as pool:
                    pool_outs = pool.starmap(chem_from_rs,pool_inputs)                    
                print('Find number of chems: ', len(pool_outs),pool_outs[0:3],'. Next build SSH-aerosols.',pathSSH_rdc_abs)

                # build ssh
                pool_inputs = [(i, pathSSH_rdc_abs) for i in pool_outs]
                with Pool(int(ncpu/2)) as pool:
                    pool_outs = pool.starmap(build_ssh_from_chem, pool_inputs)
                    
                print('Finish build SSH-aerosols. Next run simulations.',len(pool_outs))
                pool_inputs = []
                # run simulations
                nps = len(pool_outs)
                for i in range(nps):
                    for j in range(nnml):
                        # path_ssh, namelist, path_ref_res
                        pool_inputs.append((j,pool_outs[i],nmls[j],rslts[j]))
                with Pool(ncpu) as pool:
                    aerrs = pool.starmap(run_sim_with_nml, pool_inputs)
                    
                # get errs
                aerrs = np.array(aerrs).reshape(nps,nnml,2)
                tmp = [np.average(aerrs[i], axis = 0)[0] for i in range(nps)]

                # sort a_item
                a_item = [x for _,x in sorted(zip(tmp,a_item))]
                
                with open ('{:s}/reduction_order.org.{:d}'.format(path_sav_rec,iauto),'w+') as f:
                    f.write('# Reduction order: species from {:s}\n'.format(IDchemPre))
                    for i in a_item: f.write(i+'\n')
                with open ('{:s}/reduction_order.org.err.{:d}'.format(path_sav_rec,iauto),'w+') as f:
                    for i in range(n1):
                        f.write('{:s} {:}\n'.format(a_item[i],tmp[i]))
                print('Get the reduction order: {:d}. Saved in the file: {:s}/reduction_order.org.{:d}'.format(n1, path_sav_rec,iauto))

                # check errs
                if min(tmp) == -1: raise ValueError('Error in reduction order. -1 exists.', path_sav_rec,iauto)
                
                # move fast species to the front of the list
                print('Move species with kinetic >= 1 in front of the list.')
                tmp = [[],[]] # name, index in a_item
                for i in rc:
                    tag = 0
                    if 'KDEC' in i.rate.str: tag = 1
                    else:
                        i.rate.update_value()
                        if i.rate.pyformat and i.rate.pyformat > 1E0: tag = 1
                    if tag:
                        for j in i.reactants:
                            if j in a_item: 
                                if j not in tmp[0]:
                                    tmp[0].append(j)
                                    tmp[1].append(a_item.index(j))
                tmp = [x for _,x in sorted(zip(tmp[1],tmp[0]))]
                print('Find fast species: ',len(tmp),tmp)
                a_item = tmp + [x for x in a_item if x not in tmp]

                # save list
                iforder = '{:s}/reduction_order.use.{:d}'.format(path_sav_rec,iauto)
                with open (iforder,'w+') as f:
                    f.write('# Reduction ordered with fast species. {:s}\n'.format(IDchemPre))
                    for i in a_item: f.write(i+'\n')
                frec.write('rdc order is generated: {:}\t{:d}\n'.format(iforder,len(a_item)))

                print('Get the reduction order: {:d}. Saved in the file: {:s}/{:s}'.format(len(a_item),path_sav_rec,iforder))
        else:
            raise ValueError('ATP: train_key not found. ', train_key)

        print('prepare time: ',time.perf_counter()-t0, flush=True)
        
        # only init if ntree == 0, else read the saved value
        #if ntree == 0 and tag_iauto: idx = 0 # traverse

        while a_item != []:
            # search from the end of reaction list
            t1 = time.perf_counter()

            # save a_item in case stop
            with open ('{:s}/rdc_order.tmp'.format(path_sav_rec),'w+') as f:
                json.dump(a_item, f)

            # target name
            sname = a_item[0]

            # reduced by reaction/species
            if train_key == 'reaction':
                sys.exit('Not used check !!!')
                # get index in reaction list
                for i in list(reversed(range(len(rc)))):
                    if rc[i].status and rc[i].toSSH() == sname: 
                        itm = i
                        rcn = rc[itm] # reversed order

                if not rcn.status:# or itm in frorcn:
                    idx += 1
                    continue
                else: # get related reactions/species
                    irel_rcn = [itm]
                    sname = rcn.toSSH()
                    rel_sps = []
                    for i in rcn.reactants+rcn.products:
                        if i in NokeepSp: continue
                        elif i in primaryVOCs: continue
                        elif i in frozen_species:continue
                        elif i in frosps: continue
                        elif not sp[sps.index(i)].status: continue
                        elif not sp[sps.index(i)].organic: continue
                        else: rel_sps.append(i)

            elif train_key == 'species':
            
                if sname not in sps:
                    print('pop species: not in list.',sname)
                    a_item.pop(0) # not find species 
                    continue

                rel_sps = [sname]
                irel_rcn = []
                for i in range(len(rc)):
                    if not rc[i].status: continue
                    elif i in frorcn: continue
                    elif sname in rc[i].products:
                        irel_rcn.append(i)
                    elif sname in rc[i].reactants:
                        # if no products
                        if [j for j in rc[i].products if sp[sps.index(j)].organic] == []:
                            irel_rcn.append(i)
                print('\nsearch reduction related to: ',nval, rel_sps, len(irel_rcn),' rcns. Left: ', len(a_item)) #ddd
                
                # prepare cuts values for rm
                if ncpurm:
                    n = combinations_number(min(nlimtar,len(irel_rcn)))
                    for i in range(1,10):
                        if n/(i+1) < ncpurm: # find the interval
                            tmp = count_cut(int(round_up(n/(i+1)))) #cuts['rm']
                            #print('rm - update val,tot,cut: ',i,n,n/(i+1),int(round_up(n/(i+1))),len(cuts['rm']))
                            break
                else: tmp = [None]
                
                # generate cuts
                cuts = {}
                for i in strategy_types:
                    if i in ['rm','rm1','rmp']: cuts[i] = tmp # for rms
                    elif i != 'lp':
                        n = ncpus[i]
                        if i == 1: cuts[i]=[None] # default
                        else: cuts[i] = count_cut(n) #['rp','jp']: cuts[i]=[[1,2],[2,2]] #count_cut(2)

                # cuts from lp
                if 'lp' in strategy_types:
                    n = ncpu
                    for i in list(cuts.keys()): n -= len(cuts[i]) # get the number left for lp
                    if n <= 0: raise ValueError('no cpu left for lp.',n, ncpu, ncpus, cuts, sname)
                    else: 
                        n = min(n,ngas) # cpu number need to less than ngas
                        cuts['lp'] = count_cut(n)
                        #print('lp - update cuts: ',n,ngas)

                # check total numbers
                n = 0
                for i in list(cuts.keys()): n += len(cuts[i])
                #print('total n: ',n,ncpu)#,cuts)
                if n > ncpu:
                    for i in list(cuts.keys()): print('cuts ',i,len(cuts[i]),ncpus[i],cuts[i])
                    raise ValueError('number of processors exceeds limit. n > ncpu',n,ncpu)
                
            # parallel: search possible reduction and compile ssh with new chem
            pool_inputs = [] # rc, sp, itype, rel_sps, irel_rcn, frozen
            for itype in strategy_types:
                # test all pb at one step
                #if tag_all and itype in ['lp','rp','jp']: isps = [-1]
                #else: isps = rel_sps
                isps = rel_sps
                # for test
                #pool_inputs.append((rc,sp,itype,isps,irel_rcn,frozen_species,kept_species,path_ssh_chem,lptypes,None))
                #continue
                # use cut
                if cuts[itype] == [None]:
                    pool_inputs.append((rc,sp,itype,isps,irel_rcn,frozen_species,kept_species,path_ssh_chem,lptypes,None))
                else: # check type
                    if itype in ['rm','rm1','rmp','rp','jp','lp']:
                        for cut in cuts[itype]:
                            pool_inputs.append((rc,sp,itype,isps,irel_rcn,frozen_species,kept_species,path_ssh_chem,lptypes,cut))
                    else: #['da','rs']
                        pool_inputs.append((rc,sp,itype,isps,irel_rcn,frozen_species,kept_species,path_ssh_chem,lptypes,None))

            #ttest0 = time.perf_counter()
            #print('--- Time for initialization: ',sname, ttest0-t1,irel_rcn,rel_sps,ntry, flush=True)

            with Pool(ncpu) as pool:
                paths_all = pool.starmap(search_reduction, pool_inputs)

            # process paths_all
            paths,chems,changes,stgys,ncur = [],[],[],[],{} # reshape
            for i in strategy_types: ncur[i] = 0
            for item in paths_all:
                if item[1] == []: continue # no path, no reduction found
                n = len(item[1])
                #if item[0] in ['jp','rp']: print('Find reduction: ',item[0],n)
                for i in range(n):
                    stgys.append(item[0]) # strategys
                    chems.append(item[1][i]) # chem paths
                    changes.append(item[2][i]) # changes
                    ncur[item[0]] += 1 # count
                    #if stgys[-1] in ['jp','rp']: print('     ',i+1,stgys[-1],chems[-1],changes[-1]) #ddd

            # need initial check pre case
            if tag_tchemP:
                path_chem='rm-pre'
                to_SSH_sets(path_ssh_chem,path_chem,rc,sp,'10')
                stgys.append('rm')
                print('Test current IDchem ',path_ssh_chem,path_chem)
                chems.append(path_ssh_chem+'/'+path_chem)
                changes.append(['No change.'])

            # find at least one possible reduction
            if len(chems) >= 1 : ntry += 1
            else: # no reduction found, go to next
                frosps += rel_sps
                frorcn += irel_rcn
                print('pop species: no rdc attempt.',sname)
                a_item.pop(0)
                continue

            #ttesta = time.perf_counter()
            #print('--- Time for searching: ',ttesta-ttest0, ncur, flush=True)

            # parallel: compile all ssh
            pool_sshs = []
            nps = len(chems)
            for i in range(nps):
                pool_sshs.append((chems[i],pathSSH_rdc_abs))
            with Pool(int(ncpu/2)) as pool:
                paths = pool.starmap(build_ssh_from_chem, pool_sshs)

            ttestb = time.perf_counter()
            #print('--- Time for preparing ssh-aerosols: ',ttestb-ttesta, flush=True)

            # parallel: 1st time: run all simulation n_path * n_namelist
            print('start 1st simulation',nps) #ddd
            pool_nmls = []
            for i in range(nps):
                for j in range(nnml):
                    # path_ssh, namelist, path_ref_res
                    pool_nmls.append((j,paths[i],nmls[j],rslts[j]))
            with Pool(ncpu) as pool:
                aerrs0 = pool.starmap(run_sim_with_nml, pool_nmls)
            aerrs0 = np.array(aerrs0).reshape(nps,nnml,2) # shape n_path, n_namelist, 2
            
            ttestb1 = time.perf_counter()
            print('--- Time for running 1st simulations: ',ttestb1-ttestb,flush=True)

            # check first time the reduction. simulations are with all vocs
            # check err no -1
            if tag_pre_testing: # check only ref
                tmp = np.where(aerrs0[:,:,0] == -1) 
            else: # check both pre and ref
                tmp = np.where(aerrs0 == -1)
            if len(tmp[0]) != 0:
                raise ValueError('Find #-1 in 1st sml: ',len(tmp[0]),' in errs: ',paths[tmp[0][0]],', toto',tmp[1][0])
               
            fit_paths, fit_errs = [],[]
            # check rdc fit err criteria
            for c in range(nps):
                iref0, ipre0 = np.max(aerrs0[c], axis = 0)
                if iref0 <= err_ref and ipre0 <= err_pre: # max err
                    # ave err
                    iref, ipre = np.average(aerrs0[c], axis = 0)
                    if iref <= err_rav and ipre <= err_pav: # ave exceed
                        if tag_pre_testing: ierr = iref # ave error
                        else: ierr = iref0 # max error
                        if ifit_pre == 0. or (ierr-ifit_pre) <= delta_err:# delta ave err
                            # meet all criteria
                            #print('1st find fit: ',c,iref0,err_ref,ipre0,err_pre,iref,err_rav,ipre,err_pav, (ierr-ifit_pre),delta_err)
                            fit_paths.append(c)
                            fit_errs.append(ierr)

            # check pre case if need
            if tag_tchemP:
                if len(fit_paths) == 0 or nps-1 != fit_paths[-1]: # prechem does not fit error
                    raise ValueError('preChem case does not fits the error criteria. ',nps-1,iref0,err_ref,ipre0,err_pre, iref,err_rav, ipre,err_pav,ifit_pre,delta_err)
                else:
                    # save info
                    tmp = 'Test preChem: max ref {:8.5f}%, pre {:8.5f}%; ave ref {:8.5f}%, pre {:8.5f}%\n'.format(iref0*100., ipre0*100., iref*100., ipre*100.)
                    print(tmp) #ddd
                    frec.write(tmp)
                    frec.flush()
                    if ipvocs <= 0: # reset here
                        tag_tchemP = False
                        fit_paths.pop(-1) # remove from used path
                        ifit_pre = fit_errs.pop(-1) #update ifit_pre
                        ierr_pre = ifit_pre
                        print('New ifit_pre: ',ifit_pre)
                        # re-check rdc with new ifit_pre
                        if ifit_pre > 0.:
                            print('Reexamine rdc case and errs:')
                            i = 0
                            while i < len(fit_paths):
                                c = fit_paths[i]
                                iref, ipre = np.average(aerrs0[c], axis = 0)
                                if tag_pre_testing: 
                                    ierr, ipre = np.average(aerrs0[c], axis = 0) # ave error
                                else: 
                                    ierr, ipre0 = np.max(aerrs0[c], axis = 0)# max error
                                if (ierr-ifit_pre) > delta_err: # remove
                                    fit_paths.pop(i)
                                    fit_errs.pop(i)
                                    print('\tremove ',i ,c, ierr, ierr-ifit_pre,delta_err)
                                else: i += 1
                    
            # go to next
            if fit_paths == []:
                frosps += rel_sps
                frorcn += irel_rcn
                print('pop species: 1st sml not rdc fits criteria.',nps,iref0,ipre0)
                a_item.pop(0)
                continue

            nfs = len(fit_paths)
            print('Find 1st sml fit chems: ',nps,nfs) #ddd            
            # record all fit case and error
            apt_cas = {'case':[],'err':[],'size':[],'score':[]}
            icase = None
            
            if ipvocs > 0:
                # second time simulations with diff primary voc
                print('start 2nd simulation',nfs) #ddd
                pool_nmls = []
                for i in fit_paths:
                    for j in range(nnml):
                        for k in range(ipvocs): # options
                            # path_ssh, namelist, path_ref_res
                            pool_nmls.append((j,paths[i],nmls[j],rslts[j],k+1))
                with Pool(ncpu) as pool:
                    aerrs = pool.starmap(run_sim_with_nml, pool_nmls)
                #print('aerrs: ',aerrs) 
                # reshape err
                aerrs = np.array(aerrs).reshape(nfs,nnml*ipvocs,2)
                #print('af aerrs: ',aerrs)

                ttestc = time.perf_counter()
                print('--- Time for running 2nd simulations: ',ttestc-ttestb1, flush=True)
                print('--- Time for running all simulations: ',ttestc-ttestb, flush=True)

                # postprocess: analysis errors and set chem for the next run; record changes.
                # check second time the reduction.
                # check err no -1
                if tag_pre_testing: # check only ref
                    tmp = np.where(aerrs[:,:,0] == -1) 
                else: # check both pre and ref
                    tmp = np.where(aerrs == -1)
                if len(tmp[0]) != 0:
                    raise ValueError('Find #-1 in 2nd sml: ',len(tmp[0]),' in errs: ',paths[fit_paths[tmp[0][0]]],', toto',tmp[1][0])
            
                # check err and save qualified cases
                for c in range(nfs):
                    iref0, ipre0 = np.max(aerrs[c], axis = 0)
                    if iref0 <= err_ref and ipre0 <= err_pre: # max err
                        # ave err
                        iref, ipre = np.average(aerrs[c], axis = 0)
                        if iref <= err_rav and ipre <= err_pav: # ave exceed 
                            if tag_pre_testing: ierr = iref # ave error
                            else: ierr = iref0 # max error
                            if ierr_pre == 0. or (ierr-ierr_pre) <= delta_err:# delta ave err
                                # update accepted err
                                #print('2nd find fit: ',c,iref0,err_ref,ipre0,err_pre,iref, err_rav,ipre, err_pav, ierr, (ierr-ierr_pre),delta_err)
                                apt_cas['case'].append(fit_paths[c]) # case
                                apt_cas['err'].append(ierr) # error

                print('find 2nd fit chems: ',nfs,len(apt_cas['err'])) #ddd 
                
                if tag_tchemP: # check chemPre
                    if len(apt_cas['case']) == 0 or nps-1 != apt_cas['case'][-1]: # prechem does not fit error
                        raise ValueError('W/ ipvocs: preChem case does not fits the error criteria. ',nps-1,iref0,err_ref,ipre0,err_pre, iref,err_rav, ipre,err_pav,ierr_pre,ierr,ierr-ierr_pre,delta_err)
                    else:
                        # save info
                        tmp = 'Test preChem: w/ ipvocs: {:d} - max ref {:8.5f}% n pre {:8.5f}%; ave ref {:8.5f}% n pre {:8.5f}%\n'.format(ipvocs, iref0*100., ipre0*100., iref*100., ipre*100.)
                        print(tmp) #ddd
                        frec.write(tmp)
                        frec.flush()
                        
                        # resrt
                        tag_tchemP = False
                        apt_cas['case'].pop(-1) # remove from used path
                        fit_paths.pop(-1)
                        
                        # assgin ifit_pre & ierr_pre
                        ifit_pre = fit_errs.pop(-1)
                        ierr_pre = apt_cas['err'].pop(-1)
                        print('New ifit_pre & ierr_pre: ',ifit_pre,ierr_pre)
                        
                        # recheck rdc with new ifit_pre
                        if ierr_pre > 0. or ifit_pre > 0.:
                            print('Re-examine rdc case and errs w/ ipvocs:')

                            tmp = [] # save remove item
                            for i,j in enumerate(fit_paths):
                                tag = 0 # default not remove
                                if tag_pre_testing:
                                    ierr, ipre = np.average(aerrs0[j], axis = 0)
                                else:
                                    ierr, ipre = np.max(aerrs0[j], axis = 0)
                                if ifit_pre != 0. and (ierr-ifit_pre) > delta_err: # remove
                                    tag = 1
                                    print('Remove case cuz ifit_pre: ',j,ierr,ierr-ifit_pre,delta_err)
                                    
                                # get index in aerrs
                                if j in apt_cas['case']:
                                    k = apt_cas['case'].index(j)

                                    if tag_pre_testing:
                                        ierr, ipre = np.average(aerrs[i], axis = 0)
                                    else:
                                        ierr, ipre = np.max(aerrs[i], axis = 0)
                                    if ierr_pre != 0. and (ierr-ierr_pre) > delta_err: # remove
                                        tag += 2
                                        print('Remove case cuz ierr_pre: ',i,j,k,ierr,ierr-ierr_pre,delta_err)
                                # remove
                                if tag: tmp.append(j)

                            # remove if need
                            if tmp != []:
                                for i in tmp: # i = case id
                                    if i in fit_paths:
                                        j = fit_paths.index(i)
                                        fit_paths.pop(j)
                                        fit_errs.pop(j)
                                    else: raise ValueError('Cant find remove item in fit_paths.',tmp,i,fit_paths)
                                    if i in apt_cas['case']: # index in apt_cas
                                        j = apt_cas['case'].index(i)
                                        apt_cas['case'].pop(j)
                                        apt_cas['err'].pop(j)
                                
                    
            else: # not ipvocs
                apt_cas['case'] = fit_paths
                apt_cas['err'] = fit_errs

            # check if any case fits
            if apt_cas['case'] != []: # find accepted case
                for i,j in enumerate(apt_cas['case']):
                    # get new size
                    nrea1, ngas1, naer1 = get_info('{:s}/src/include/CHEMISTRY/'.format(paths[j]),'BCARY',Type='read')
                    apt_cas['size'].append([nrea1, ngas1, naer1])
                    # get size diff
                    #if nrea1 > nrea or ngas1 > ngas or naer1 > naer:
                    #    print('size after reduction increases! ',paths[j],gpt_cas['size'][-1],[nrea,ngas,naer])
                    # compute the score: default 0.1 * derr
                    tmp = max(0.1, max(0,naer-naer1)*10.+max(0,nrea-nrea1)+max(0,ngas-ngas1))
                    if tmp == 0.1: print('score size not change: ',i,j,apt_cas['size'][-1],tmp)
                    #print('compute score: ',j,naer1,nrea,ngas,tmp,apt_cas['err'][i],(err_ref - apt_cas['err'][i]),tmp*(err_ref - apt_cas['err'][i]))
                    # get delta err
                    if tag_pre_testing:
                        tmp *= (err_rav - apt_cas['err'][i])
                    else:
                        tmp *= (err_ref - apt_cas['err'][i])
                    apt_cas['score'].append(tmp)

                # if err increase and naer not change, refuse
                if tag_stage:
                     i = 0
                     info1 = ''
                     while i < len(apt_cas['case']): # check naero and ierr
                        if apt_cas['size'][i][2] < naer: # reduce # aero
                            i += 1 # go to next
                        else:
                            # check derrs
                            #print('a',i,len(apt_cas['err']),len(apt_cas['score']), len(apt_cas['case']),len(apt_cas['size']), len(fit_errs))
                            #print('b',apt_cas['case'][i],apt_cas['err'][i], apt_cas['score'][i],apt_cas['size'][i])
                            j = apt_cas['case'][i] # get case id
                            if apt_cas['err'][i] > ierr_pre or fit_errs[fit_paths.index(j)] > ifit_pre:
                                # remove cuz tag_stage
                                if ipvocs:
                                    tmp = '!Refuse w/ stage: {:d}, {:} >= {:d} or ierr {:8.5f}% > {:8.5f}% or ifit {:8.5f}% > {:8.5f}%\n'.format(j, apt_cas['size'][i], naer, apt_cas['err'][i]*100., ierr_pre*100., fit_errs[fit_paths.index(j)]*100., ifit_pre*100.)
                                else:
                                    tmp = '!Refuse w/ stage: {:d}, {:} >= {:d} or ierr {:8.5f}% > {:8.5f}%\n'.format(j, apt_cas['size'][i], naer, apt_cas['err'][i]*100., ifit_pre*100.)
                                # remove
                                apt_cas['score'].pop(i)
                                apt_cas['size'].pop(i)
                                apt_cas['err'].pop(i)
                                apt_cas['case'].pop(i)
                                istage += 1 # record
                                print(tmp,istage)
                                info1 += tmp
                            else: # derrs decrease 
                                i += 1

                if apt_cas['case'] != []: 
                    # sort by score
                    apt_cas['score'],apt_cas['size'], apt_cas['err'], apt_cas['case'] = zip(*sorted(zip(apt_cas['score'], apt_cas['size'], apt_cas['err'],apt_cas['case']),reverse=True)) # sort err, cas
                    #apt_cas['size'], apt_cas['err'], apt_cas['case'] = zip(*sorted(zip(apt_cas['size'], apt_cas['err'],apt_cas['case']))) # sort err, cas

                    # update ierr, icase
                    #ierr = apt_cas['err'][0]
                    icase = apt_cas['case'][0] #index in nps
                    #iaer = apt_cas['size'][0][2]

            if icase is not None: # accepted
                # update infos
                isgy = stgys[icase]
                npres[isgy] += 1
                
                # update errors
                ierr_pre = apt_cas['err'][0]
                ifit_pre = fit_errs[fit_paths.index(icase)]
                
                if tag_reach == 1 and err_rav != err_rav_reach:
                    if ierr_pre <= err_rav_reach and ifit_pre <= err_rav_reach:
                        print('err_rav reaches err_rav_reach: ',err_rav, ierr_pre, ifit_pre, err_rav_reach)
                        for i in range(len(input_params['err_rav'][nerr:])):
                            if input_params['err_rav'][nerr+1] == err_rav:
                                input_params['err_rav'][nerr+1] = err_rav_reach
                        err_rav = err_rav_reach
                        tag_reach = -1 # need to record
                
                nval += 1
                frosps, frorcn = [],[] #reset

                # reload rc,sp
                ifile = '{:s}/{:s}'.format(chems[icase],chems[icase].split('/')[-1])
                rc,sp = read_chem_sets(ifile+'.reactions', ifile+'.mol',
                                      speciesType='SSH',reactionType='SSH')
                sps = [i.name for i in sp]
                nrea, ngas, naer = get_info(rc,sp,'') # get scheme

                # trim
                #rc, sp, trim_sp = trim_scheme(rc, sp)
                
                # record chem files
                if nval == 1:
                    path_sav_chem_all = '{:s}/{:s}_sav'.format(path_sav_chem,IDchem)
                    create_folder(path_sav_chem_all)
                # back up reaction/species lists
                for i in 'reactions','mol':
                    os.system('cp {0:s}.{1:s} {2:s}/{3:d}.{1:s}'.format(ifile,i,path_sav_chem_all,nval))
                if not tag_pre_testing: # update pre files
                    for i in to_pre_paths:
                        os.system('cp {:s}/{:s}/aero/Organics_1.txt {:s}'.format(paths[icase],i[0],i[1]))
                        for j in range(ipvocs):
                            os.system('cp {0:s}/{1:s}.{3:d}/aero/Organics_1.txt {2:s}.{3:d}'.format(paths[icase],i[0],i[1],j+1))
                # renew idx, itm if need
                if train_key == 'reaction':
                    pass
                elif train_key == 'species':
                    # update a_item: change the rank of merged species
                    if isgy in ['jp','rp','lp']:
                        if isgy in ['jp','rp']: # may in format A->B1,B2,B3
                            tmp = changes[icase][1].replace('->',',').split(',')
                        elif isgy == 'lp':
                            tmp = changes[icase][0] # lumped species
                        # get index for all changed species
                        isnm = 0
                        for i in tmp:
                            if i in a_item and i not in sps: # get index
                                j = a_item.index(i)
                                isnm = max(isnm, j) # get max index
                                a_item[j] = None # remove species
                        for i in tmp:
                            if i in a_item and i in sps: # get index
                                j = a_item.index(i)
                                if j < isnm: # update rank if current is smaller
                                    a_item[j] = None
                                    a_item.insert(isnm,i)
                        a_item = [i for i in a_item if i is not None]

                print('Find: ',isgy, nval, 'Next: ',sname,' to ',a_item[0], ifit_pre, ierr_pre,' left n:',len(a_item))
                #if trim_sp != []: print('Trim species: trim_sp')  
                #isgy in ['rm','rm1','da','rs']

            else: # no reduction is found
                # update frosps,frorcn cuz they hve been tested
                frosps += rel_sps
                frorcn += irel_rcn
                # go next
                print('pop species: 2nd sml no fits.',sname,len(chems))
                a_item.pop(0)
                #if tag_all: a_item = [] # try once

            #ttest2 = time.perf_counter()
            #print('finish check results: ',ttest2-ttestc, flush=True)

            # prepare output
            out = []
            for c, chn in enumerate(changes):
                if not tag_record and c != icase: continue
                # only output accepted changes if not tag_record
                auto_type = stgys[c]
                if auto_type == 'lp':
                    #lump,lumppd,trim_sp = chn
                    lump,lumppd = chn
                    tmp = [lump[0]]
                    for i in lumppd: tmp.append(i)
                    #if trim_sp != []:
                    #    tmp.append('trim_scheme: {:}'.format(trim_sp))
                    out.append(tmp) # remove \t

                elif auto_type in ['jp','rp']:
                    out.append([chn[1]])
                    #if chn[2] != []: 
                    #    out[-1].append('trim_scheme: {:}'.format(trim_sp))

                elif auto_type =='da':
                    if len(chn) == 1:
                        out.append(chn)
                    else: # in one line
                        tmp = ''
                        for i in chn: tmp += i+'\t'
                        out.append([tmp])
                elif auto_type in ['rs','rsn']:
                    tmp = [chn[0]]
                    if chn[1][0] != [[]]:
                        tmp.append('As reactants:')
                        tmp += chn[1][0]
                    if chn[1][1] != [[]]:
                        tmp.append('As products:')
                        tmp += chn[1][1]
                    if chn[1][2] != [[]]:
                        tmp.append('New reactions:')
                        tmp += chn[1][2]
                    out.append(tmp)
                elif auto_type in ['rm','rmp']: out.append(chn)
                elif auto_type =='rm1': out.append(chn[1])

            # record output
            t2 = time.perf_counter()
            info = '\n\n======No.try {:d}\tNo.val {:d}\tChem size: [{:d},{:d},{:d}]\ttime: {:.1f}s======\n'.format(ntry,nval,nrea,ngas,naer,t2-t1)

            if tag_reach == -1: # need to record
                info += 'err_rav reaches try_ave_ref: {:6.3f}%, ipre {:6.3f}% & ierr {:6.3f}%\n'.format(try_ave_ref*100.,ifit_pre*100.,ierr_pre*100.)
                tag_reach = 0

            info += 'Related to {:s}: {:s}\n'.format(train_key, sname)
            info += 'Find #{:d} chems in 1st sml\tcount: {:}\nFind #{:d} qualified for 2nd sml.\n'.format(nps,ncur,nfs)
            # find case that need standard
            if len(apt_cas['case']) > 0:
                # find accpted case
                if icase is not None:
                    if ipvocs:
                        tmp = '------accept: {:s}\tref_err: {:6.3f}% / {:6.3f}%\tnpres: {:}------\n'.format(isgy,ierr_pre*100.,ifit_pre*100.,npres)
                    else:
                        tmp = '------accept: {:s}\tref_err: {:6.3f}%\tnpres: {:}------\n'.format(isgy,ierr_pre*100.,npres)
                        
                    if len(a_item) > 0: 
                        tmp += '  Next {:s} to {:s}, #{:d} left.\n\n'.format(sname, a_item[0],len(a_item))
                    tmp += '{:s}\n'.format(chems[icase])
                    
                    # record 1st round errors
                    if ipvocs: 
                        tmp += '1st round errs, icase {:d} w/ ave err {:}:\n'.format(icase, np.average(aerrs0[icase], axis = 0)*100.)
                    iref0, ipre0 = np.max(aerrs0[icase], axis = 0)
                    j = int(np.where(aerrs0[icase]==iref0)[0][0])
                    tmp += 'max_err_ref: {:6.3f}%\tloc: {:s}\n'.format(iref0*100.,cases[j])
                    if not tag_pre_testing:
                        j = int(np.where(aerrs0[icase]==ipre0)[0][0])
                        tmp += 'max_err_pre: {:6.3f}%\tloc: {:s}\n'.format(ipre0*100.,cases[j])
                    # record 2nd round errors
                    if ipvocs:
                        ifit = fit_paths.index(icase)
                        tmp += '2nd round ave errs, ifit {:d} w/ ave err {:}\n'.format(ifit,np.average(aerrs[ifit], axis = 0))
                        for k in range(ipvocs):
                            iref = [aerrs[ifit][k+j*ipvocs][0]*100. for j in range(nnml)]
                            #print('iref: ',len(iref),iref)
                            if not tag_pre_testing:
                                ipre = [aerrs[ifit][k+j*ipvocs][1]*100. for j in range(nnml)]
                                tmp += '\t#{:d}\tref ave {:6.3f}%\tmax {:6.3f}%\tpre ave {:6.3f}%\tmax {:6.3f}%\n'.format(k+1,sum(iref)/len(iref),max(iref),sum(ipre)/len(ipre),max(ipre))
                            else: tmp += '\t#{:d}\tref ave {:6.3f}%\tmax {:6.3f}%\n'.format(k+1,sum(iref)/len(iref),max(iref))
                    # accepted change
                    if tag_record: j = icase
                    else: j = 0
                    for i in out[j]: tmp += '{:}\n'.format(i)
                    #if trim_sp != []: tmp += '! trim_scheme: find species {:}\n'.format(trim_sp)
                else:
                    tmp = '------no accepted case from valid reductions.-----\n'

                # rec in fall & fuse
                for f in [frec,fall,fuse]:
                    f.write(info) # general info
                    f.write(tmp)
                    if tag_stage and info1 != '': 
                        f.write(info1)
                        info1 = '' # reset
                    # record all
                    f.write('\nSelected from {:d} cases.\n All: {:}, {:}. icase: {:}\n'.format(len(apt_cas['case']),apt_cas,[stgys[i] for i in apt_cas['case']], icase))
                    f.flush()
            else:
                fall.write(info) # general info
                fall.write('------no reduction fits error thresholds.------\n')
                fall.flush()

            # plain record
            if tag_record: #ddd
                fall.write('\n------For all changes in {:d} nmls------\n'.format(nnml))
                for i in range(nps):
                    # changes
                    fall.write('---{:d}-{:d}\tStrategy: {:s}---{:s}\n'.format(ntry,i,stgys[i],chems[i]))
                    for j in out[i]:
                        fall.write('{:}\n'.format(j))
                    fall.write('---errors---\n')
                    # errors
                    iref0, ipre0 = np.max(aerrs0[i], axis = 0)
                    if tag_pre_testing:
                        j = int(np.where(aerrs0[i]==iref0)[0][0])
                        tmp = 'max: {:6.3f}%\tloc: {:s}\n'.format(iref0*100.,cases[j])
                        iref0, ipre0 = np.average(aerrs0[i], axis = 0)
                        tmp = 'ave: {:6.3f}%\t'.format(iref0*100.) + tmp
                        iref0, ipre0 = np.median(aerrs0[i], axis = 0)
                        tmp = 'median: {:6.3f}%\t'.format(iref0*100.) + tmp
                        fall.write(tmp)
                    else:
                        j = int(np.where(aerrs0[i]==iref0)[0][0])
                        tmp = 'ref: {:6.3f}%\tloc: {:s}\n'.format(iref0*100.,cases[j])
                        j = int(np.where(aerrs0[i]==ipre0)[0][0])
                        tmp += 'pre: {:6.3f}%\tloc: {:s}\n'.format(ipre0*100.,cases[j])
                        fall.write(tmp)
                        fall.write('---namelists---\n')
                        for j in range(nnml):
                            fall.write('loc: {:s}\tref: {:6.3f}\tpre: {:6.3f}\n'.format(cases[j], \
                                        aerrs0[i][j][0]*100.,aerrs0[i][j][1]*100.))
                        fall.write('\n')
                        if ipvocs:
                            for i in range(nfs):
                                for j in range(nnml):
                                    for k in range(ipvocs): # options
                                        s = k+j*ipvocs
                                        print('2nd err: ',i,j,k,s,aerrs[i][s][0],aerrs[i][s][1], paths[fit_paths[i]],nmls[j])
                fall.flush()

            # check time if need out
            if tlim and t2 - tbeg > tlim: # in second
                tmp = 'Stop cuz time run out. Total time: {:.1f}s, time limitation: {:.1f}s\n'.format(t2-tbeg,tlim)
                frec.write(tmp)
                print(tmp)
                tag_tlim = True
                break
                
            # save tree or record for break point
            if (tag_tree and len(apt_cas['case']) >= 1):
                # path: ntree-iauto-idx-ntry
                tree_path = '{:s}/{:d}/{:d}/{:d}'.format(pathTree,ntree,iauto,ntry)
                tmp = [] # save info
                for k,i in enumerate(apt_cas['case']):
                    if i == icase: continue # only record the one not used
                    # save chems
                    ifile = chems[i].split('/')[-1]
                    ipath = '{0:s}/{1:d}-{2:s}'.format(tree_path,k,ifile)
                    if tmp == []: # init save path
                        # generate folder
                        os.makedirs(tree_path, exist_ok=True)
                        # record the last branch: the reudction tree to try next
                        next_tree.append(ipath)
                        print('Add in next tree:',ipath,len(next_tree)) # test
                    # copy chem files
                    for j in ['reactions','mol']:
                        os.system('cp {0:s}/{1:s}.{2:s} {3:s}.{2:s}'.format(chems[i],ifile,j,ipath))
                    # prepare some output for info
                    auto_type = stgys[i]
                    chn = changes[i]
                    if auto_type == 'lp':
                        lump,lumppd = chn
                        out = [lump[0]]
                        for j in lumppd: out.append(j)
                    elif auto_type in ['jp','rp']:
                        out = [chn[1]]
                    elif auto_type == 'da':
                        if len(chn) == 1: out = chn
                        else: # in one line
                            out = ['']
                            for j in chn: out[0] += j+'\t'
                    elif auto_type in ['rs','rsn']:
                        out = [chn[0]]
                        if chn[1][0] != [[]]:
                            out.append('As reactants:')
                            out += chn[1][0]
                        if chn[1][1] != [[]]:
                            out.append('As products:')
                            out += chn[1][1]
                        if chn[1][2] != [[]]:
                            out.append('New reactions:')
                            out += chn[1][2]
                    elif auto_type in ['rm','rmp']: out = chn
                    elif auto_type == 'rm1': out = chn[1]

                    tmp.append([ipath, i, apt_cas['err'][k], apt_cas['size'][k], auto_type, out])

                if icase != None:
                    # add the case with no change
                    ifile = '99-{:s}-0'.format(stgys[icase])
                    to_SSH_sets(tree_path,ifile,rc,sp,'00')
                    ipath = '{:s}/{:s}'.format(tree_path,ifile)
                    tmp.append([ipath,99,stgys[icase]])
                
                # save info related
                if tmp != []:
                    save_data = {'idx':idx,'nval':nval,
                                 'nerr': nerr, 
                                 'npres': '{:}'.format(npres), # save as string
                                 'cases': tmp, # save cases
                                 'tag_stage':tag_stage,
                                 'ifit_pre':ifit_pre,
                                 'ierr_pre':ierr_pre}
                    # record icase info
                    if icase == None: save_data['icase'] = info1
                    else: save_data['icase'] = [chems[icase],stgys[icase],ierr,iaer]

                    with open('{:s}/log'.format(tree_path), 'w') as f:
                        json.dump(save_data, f, indent=4)
            
            #ttest3 = time.perf_counter()
            #print('output and finish this step: ',ttest3-ttest2, flush=True)

        # save chem
        nrea, ngas, naer = get_info(rc,sp,'')
        
        # 1/0 with/without viz file
        if naer <= nviz: 
            to_SSH_sets(path_sav_chem,IDchem,rc,sp,'20')
        else: 
            to_SSH_sets(path_sav_chem,IDchem,rc,sp,'10')

        # record
        t3 = time.perf_counter()
        tmp = '\nEND. Total run: {:d}\tValid rdc: {:d}\tTrain time: {:.1f}s\n'.format(ntry,nval,t3-t0)
        tmp += 'npres: {:}'.format(npres)
        tmp += '\nFinal scheme No.reaction: {:d}\tNo.gas: {:d}\tNo.aerosol: {:d}\n'.format(nrea, ngas, naer)
        tmp += '------------------------------------------\n'
        for f in [fall,fuse,frec]:
            f.write(tmp)
            f.flush()
            
        # pre_testing and record if need
        if nval or emax_ind == []:
            emax_ind, emax, eave, eave_max = auto_testing(Setups = ATSetups, 
                                                          IDchem = IDchem,
                                                          chempath = path_sav_chem,
                                                          ind = prelocs,
                                                          trd_max = 1,
                                                          out_file = False,
                                                          sav_soapath = False)
            # change format from [[y,x,m]] to str
            if emax_ind != []: emax_ind = 'm{:d}y{:d}x{:d}'.format(emax_ind[0][2], emax_ind[0][0], emax_ind[0][1])
        t4 = time.perf_counter()
        for f in [fall,fuse,frec]:
            f.write('-------------------------\n')
            f.write('Testing AT{:s}: err loc: {:s}\terr max: {:6.4f}\t err ave: {:6.4f}\t err ave max: {:6.4f}\ttest time: {:.1f}s\n'.format(IDchem, emax_ind, emax, eave, eave_max,t4-t3))
            f.flush()

        # save results
        if nval or tag_tlim: # obtain new results and put into path_sav_res
            # obtain new results and put into path_sav_res
            run_SSH(IDchem,path_sav_chem,
                    pathSSH_rdc,path_sav_res,
                    IDchem.replace(prefix,''),
                    ResultFolder,Ttot,DeltaT,
                    [locations[0:nlocs],locs],
                    Tnow,'fast_compile',initfile,
                    pathInitFiles,ipvoc=ipvocs)
            # update iforder
            iforder = None

        # reset tag_rm1
        if nval and tag_rm1 == 1:
            tag_rm1 = 0
            for i in range(len(input_params['tag_rm1'][nerr:])):
                if input_params['tag_rm1'][nerr+i] == 2:
                    input_params['tag_rm1'][nerr+i] = 0
            print('tag_rm1 reset to 0. nval: ', nval, nerr)
            print('Current tag_rm1: ', input_params['tag_rm1'][nerr:])

        # check if stop cuz error too large
        if emax > try_max_ref or eave > try_ave_ref:
            tag_iauto = 0
            frec.write('stop cuz errs too large: max/ave {:6.4f}/{:6.4f} & limitation: {:6.4f}/{:6.4f}\n'.format(emax,eave,try_max_ref,try_ave_ref))

        # check if redo cuz stage changes
        if tag_stage and eave >= err_stage_change:
            frec.write('redo cuz errs >= err_stage_change: max/ave {:6.4f}/{:6.4f} & limitation on average {:6.4f}\n'.format(emax,eave,err_stage_change))
            tag_redo = 1
            
        # check if stop cuz time limitation
        if tag_tlim:
            tmp = '!!!Time is up. Tlim: {:.1f}\n'.format(tlim)
            re_data = {}
            # save restart file
            # a_item, IDchemPre, prePath, iauto
            for i in ['ntree','nerr','ierr_pre','ifit_pre']:
                exec('re_data["{0:s}"]={0:s}'.format(i))
            re_data['iauto'] = iauto+1
            re_data['train_order'] = '{:s}/reduction_order.for.{:d}'.format(path_sav_rec,iauto+1)
            # save list
            with open (re_data['train_order'],'w+') as f:
                f.write('# Reduction ordered with fast species. {:s}\n'.format(IDchem))
                for i in a_item: f.write(i+'\n')
            # update chem
            re_data['IDchemPre'] = IDchem
            re_data['prePath'] = curPath
            re_data['preChemPath'] = path_sav_chem
            re_data['refile'] = refile

            # incase changed
            re_data['refPath'] = refPath
            re_data['IDchemFake'] = IDchemFake
            re_data['fakePath'] = fakePath

            # save refile
            n = [int(i.split('.')[1]) for i in os.listdir(path_sav_rec) if 'restart.' in i]
            if n == []: n = 0
            else: n = max(n)+1
            refile = '{:s}/restart.{:d}'.format(path_sav_rec,n)
            with open(refile, 'w') as f:
                json.dump(re_data, f, indent=4)
            print('For restart: save current rdc in ',refile)

            IDchem = None # not run testing
            tag_iauto = 0 # exit cuz time out

        elif tag_redo and nerr < npara -1 : # redo this reduction cyle with nerr+=1
            # get new nerr with the last value
            nerr = npara - 1
            tmp = ' Redo this reduction cycle with the last nerr {:d}.\n'.format(nerr)
            
        # go to next reduction cycle
        elif nval == 0 or (tag_efficient and nval <= tag_efficient): 

            # add elementory-like treatment
            if len(rc) <= tag_rm1:
                tmp = 'Go to next step. add rm1 cuz No. of reactions <= tag_rm1. {:d} <= {:d}.\n'.format(len(rc),tag_rm1)
                tag_rm1 = 1 # for next round

            # update w/o aerosol-oriented treatment if tag_stage = 2
            elif tag_stage == 2 and istage:
                tmp = 'Stay this step w/o tag_stage. istage: {:d}\n'.format(istage)
                tag_stage, istage = 0, 0
                
            # update error tolerance index
            elif nerr < npara-1: 
                tmp = 'Go to next step. nerr += 1, err index {:d} < {:d}.\n'.format(nerr,npara-1)
                nerr += 1

            # update reduction tree if activated !!!
            elif tag_tree and ntree < ntree_ctl and next_tree != []: # no record reduction tree
                tmp = 'No tree is recorded. ntree_ctl:{:d},tag_tree:{:d},next_tree:{:}.\n'.format(ntree_ctl,tag_tree,next_tree)
                # at least find one branch
                inext_tree = next_tree[-1]
                icase = inext_tree.split('/')[-1] # [n]-stgy-[c]
                ipath = inext_tree[:-len(icase)]  # ntree/iauto/ntry
                ind = int(icase.split('-')[0])

                # get all files in the same path
                tree_files = list(sorted(os.listdir(ipath)))

                # load new subtree root case if both log and reaction files are exist
                if icase+'.reactions' in tree_files and 'log' in tree_files:
                    ntree += 1 # update
                    IDchem = '{:s}-{:d}-init'.format(prefix+tail,ntree) # new IDchem name

                    # record
                    tmp += 'Load subtree: {:s} (No.{:d}) as {:s}.\n'.format(inext_tree,len(next_tree),IDchem)
                    # load/update save_data
                    with open('{:s}/log'.format(inext_tree), 'r') as f:
                        save_data = json.loads(f.read())
                    if ind == 99: save_data['cases'] = save_data['cases'][-1]
                    else: save_data['cases'] = save_data['cases'][ind]
                    save_data['IDchem'] = IDchem

                    # move chems from pathTree to path_sav_chem_all IDchem
                    create_folder('{:s}/{:s}'.format(path_sav_chem,IDchem))
                    for i in ['reactions','mol']:
                        shutil.move('{:s}.{:s}'.format(icase,i), '{0:s}/{1:s}/{1:s}.{2:s}'.format(path_sav_chem,IDchem,i))
                    # save log
                    with open('{:s}/init.log'.format(path_sav_chem_all), 'w') as f:
                        json.dump(save_data, f, indent=4)
                else:
                    raise FileNotFoundError('not find the icase file or log file in next_tree',icase,ipath,inext_tree,tree_files)

                # update next_tree
                tree_files = list(sorted([i for i in tree_files if '.reactions' in i and icase not in i]))
                print('current next_tree: ',next_tree[-1])
                if tree_files != []: next_tree[-1] = '{:s}/{:s}'.format(ipath,tree_files[0][:len('.reactions')]) # remove '.reactions'
                else: next_tree = next_tree[:-1] # remove the tree point
                print('next_tree is updated: ',next_tree[-1],ntree)

                # update parameters
                istage = 0, {'err':0,'aer':0} # reset stage
                iauto = -1
                idx = save_data['idx']
                # not repeat 
                if ind == 99: idx += 1
                #nval,npres = save_data['nval'],save_data['npres']
                nerr = save_data['nerr']
                tag_stage = save_data['tag_stage']
                ierr_pre = save_data['ierr_pre']
                ifit_pre = save_data['ifit_pre']
                nval = -1 # for update IDchem

            #elif not tag_all:
            #    tag_all = True
            #    tmp = 'tag_all = 1.\n'

            
            else:  # exit
                tmp = 'No next reduction cycle.\n'
                tag_iauto = 0

            if tag_reach > 0 :
                tmp += 'tag_reach still actived: {:d}. ave err: {:f} -> {:f}\n'.format(tag_reach, err_rav, try_ave_ref)

            # print & record
            print(tmp)
            frec.write('\n'+tmp)
            
        # record
        t5 = time.perf_counter() - t0
        
        if tag_tree:
            tmp = 'No.tree {:d}\tNo.cyl {:d}\tNo.try {:d}\tNo.val {:d}'.format(ntree,iauto,ntry,nval)
        else:
            tmp = 'No.cyl {:d}\tNo.try {:d}\tNo.val {:d}'.format(iauto,ntry,nval)
            
        for f in [fall,fuse,frec]:
            f.write('{:s}\ttime {:.1f}\n'.format(tmp,t5))
            if tag_iauto == 0: f.write('!!!Stop training.\n')
            f.write('===============================================\n')
            f.flush()

        for f in [fall,fuse]: f.close()
        
        # number of cycle + 1
        iauto += 1

    frec.close()

    return IDchem, path_sav_chem

def search_reduction(reactions0, species0, strategy, rel_sps, irel_rcn, frozen, kept, path_ssh, lumptypes, cut):

    # check no searching
    if irel_rcn == []: # no reaction
        if strategy in ['rm','rm1','rmp']:
            return [strategy,[],[]]
    elif rel_sps == []: # no species
        if strategy in ['da','rs','rsn','lp','rp','jp']:
            return [strategy,[],[]]

    # searching
    if strategy in ['da','rs','rsn','rm','rm1','rmp']:
        return search_simple_reduction(reactions0, species0,\
                                       strategy, rel_sps, irel_rcn, \
                                       path_ssh, cut, kept)
    elif strategy in ['lp','rp','jp']:
        if rel_sps == [-1]: # search from all species # caution when use this option !!
            return search_complex_reduction(reactions0, species0, strategy, [], \
                                            frozen, kept, path_ssh, cut, lumptypes = lumptypes)
        else:
            return search_complex_reduction(reactions0, species0, strategy, rel_sps, \
                                            frozen, kept, path_ssh, cut, lumptypes = lumptypes)
    
def search_simple_reduction(reactions0, species0, strategy, rel_sps, irel_rcn, path_sav, cut=None, kept = []):
    """1) search possible related reduction based on current strategy
       2) output in the ssh-aerosol
       3) compile ssh-aerosol with the new """

    #ttest0 = time.perf_counter()

    # check kept species
    if kept != []: rel_sps = [i for i in rel_sps if i not in kept]
    
    species = deepcopy(species0)
    reactions = deepcopy(reactions0)

    # get the entire species list
    sps = [i.name for i in species]

    # find all possibility
    if strategy == 'da':  # remove aerosol of provided species
        rm_sps = [i for i in rel_sps if species[sps.index(i)].condensed]
        tar = rm_sps
        
    elif strategy in ['rs','rsn']: # remove species
        rm_sps = [i for i in rel_sps if species[sps.index(i)].organic]
        tar = rm_sps
        
    elif strategy == 'rm': # remove reactions
        tar = irel_rcn
        
        # check primaryVOC, if has, only remove this product
        irel_rcn_pvoc, rel_rcn_pvoc = [],{}
        for i in irel_rcn:
            if set(primaryVOCs) & set(reactions[i].reactants):
                irel_rcn_pvoc.append(i)
        
    elif strategy == 'rmp': # remove reaction products
        tar = [] # rcn index, pd index
        for i in irel_rcn:
            for j in range(len(reactions[i].products)):
                if species[sps.index(reactions[i].products[j])].organic:
                    tar.append([i,j])
                    
    elif strategy == 'rm1': # remove partial reactions
        # change reactions to elementary reactions
        rm_rcn_ind, rm_rcn_spr = [], []
        # search
        for i in irel_rcn:
            n = len([j for j in reactions[i].products if species[sps.index(j)].organic])
            if n > 1:
                pds = []
                for j in reactions[i].reactants + reactions[i].products:
                    pds.append(species[sps.index(j)]) # get species list
                rm_rcn_ind.append(i) # reactions that need to be replaced
                rm_rcn_spr.append(reaction_seperate([reactions[i]],pds,Type = 'single'))
        tar = []
        if rm_rcn_ind != []:
            for i in range(len(rm_rcn_ind)): 
                for j in range(len(rm_rcn_spr[i])):
                    tar.append([rm_rcn_ind[i],j])
    else:
        print(strategy)
        raise NameError('ATP: search_simple_reduction: strategy not recognized.')

    # find all combination
    ntar = len(tar)
    if tar == []: # no element
        rm_items = []
    elif ntar == 1: # one possibility
        if isinstance(cut,list) and len(cut)==2 and cut[0] != 1:
            rm_items = []
        else:
            rm_items = [tar]
    else:
        if ntar <= nlimtar: #combinations_number(ntar) <= nlimcmb:
            # need to trim tar
            # get all possibility
            npb = combinations_index(tar)
        else: # too many possibility (ntar = 5,6: 80, 192)
            # try 1 and all
            if cut is None or cut[0] == 1: print('ntar>nlimtar',ntar,strategy) #ddd
            npb = [tar]+[[i] for i in tar]

            # add other possibility
            i = 1
            while ntar-i > 1:
                n = list(combinations(tar,ntar-i))
                if len(n) < nlimcmb - len(npb):
                    npb += n
                else: break #npb += n[0:nlimcmb-len(npb)]
                i += 1
                
            # add rm with the same reactants
            if strategy == 'rm':
                tmp = [[],[]] # reactants, inds
                for i in tar:
                    j = list(sorted(reactions[i].reactants))
                    if j in tmp:
                        tmp[1][tmp[0].index(j)].append(i)
                    else: # add new
                        tmp[0].append(j) #rcn
                        tmp[1].append([i]) #ind
                for i in range(len(tmp[0])):
                    if len(tmp[1][i]) != 1 and len(tmp[1][i]) != ntar:
                        if tmp[1][i] not in npb:
                            npb.append(tmp[1][i])
                            print('Add rm test w same reactants ',tmp[0][i],tmp[1][i])
                        
        nc = len(npb) # total

        # cut
        if isinstance(cut,list) and len(cut)==2:
            i,j = cut
            if nc <= j: # one for one
                if i <= nc: rm_items = [npb[i-1]]
                else: rm_items = []
            else: # need distribute jobs
                n = int(nc/j) # interval
                l = nc%j # rest
                if i <= l: rm_items = npb[(n+1)*(i-1):(n+1)*i]
                else: rm_items = npb[l+n*(i-1):l+n*i]
        else: rm_items = npb
        #print('rm + cut: ',cut, isinstance(cut,list), len(cut), rm_items) #ddd

    paths, changes = [], [] # init
    for c in range(len(rm_items)):

        # generate ssh if not contain
        #path_now = '{:s}/{:d}'.format(path_ssh,c)
        #set_SSH(None,None, path_now, pathSSH_sav=path_sav_ssh, mode='fast')

        # update chem files
        tar = rm_items[c] # targets
        if strategy == 'da':
            for i in tar:
                species[sps.index(i)].condensed = False
        elif strategy in ['rs','rsn']:
            out_info = [[],[],[]] # as reactants, as products, new as products
            for i in tar:
                species[sps.index(i)].status = False
                # update reaction list
                for j in range(len(reactions)):
                    if i in reactions[j].reactants: 
                        out_info[0].append([j,reactions[j].toSSH()])
                        reactions[j].status = False 
                    if i in reactions[j].products:
                        # save original reactions
                        out_info[1].append([j,reactions[j].toSSH()])
                        n = 0
                        rt = 0. # ratios
                        while n < len(reactions[j].products):
                            js = reactions[j].products[n]
                            if js == i:
                                reactions[j].products.pop(n)
                                rt += reactions[j].ratiosPD.pop(n)
                            else: n += 1
                        if strategy == 'rsn': # normolize the products
                            if rt > 0.: # need to change rate
                                if rt >= 1.: # negative rate
                                    reactions[j].status = False
                                else: 
                                    for k in range(len(reactions[j].ratiosPD)):
                                        reactions[j].ratiosPD[k] /= (1.-rt)
                                    # update rt
                                    reactions[j].rate.str += ('*{:6.4E}'.format(1.-rt))
                                    reactions[j].rate.update()
                        out_info[2].append([j,reactions[j].toSSH()])
                            
        elif strategy == 'rm':
            for i in tar:
                if i in irel_rcn_pvoc:
                    rel_rcn_pvoc[i] = [reactions[i].products,reactions[i].ratiosPD]
                    # find indexs
                    for j in rel_sps:
                        if j in reactions[i].products:
                            # remove
                            k = reactions[i].products.index(j)
                            reactions[i].products.pop(k)
                            reactions[i].ratiosPD.pop(k)
                else:
                    reactions[i].status = False
                
        elif  strategy == 'rmp':
            out_info = []
            for i in tar:
                j,k = i
                pd = reactions[j].products.pop(k)
                rt = reactions[j].ratiosPD.pop(k)
                out_info.append([pd,rt]) # product and ratio
        
        elif strategy == 'rm1':
            # get index for reactions that need to be changed
            ind0, ind1, ind_new = [], [], [] # removed, added, real_added
            for j in tar:
                if j[0] not in ind0: # record rcn need to be disabled
                    ind0.append(j[0])
                    ind1.append([])
                k = ind0.index(j[0])
                if j[1] not in ind1[k]:
                    ind1[k].append(j[1]) # need add
            ind0, ind1 = zip(*sorted(zip(ind0, ind1))) # sort
            for i,j in enumerate(ind0):
                reactions[j].status = False # disable old rcn
                for s,k in enumerate(ind1[i]):
                    reactions.insert(j+s+1,rm_rcn_spr[rm_rcn_ind.index(j)][k]) # add new rcn after disabled rcn
                    ind_new.append(j+s+1) # record read added index

        # save chem
        if cut: path_chem = '{:s}-{:d}_{:d}'.format(strategy,cut[0],c)
        else: path_chem = '{:s}-{:d}'.format(strategy,c)
        
        #reactions, species, trim_sp = trim_scheme(reactions, species)
        imark = to_SSH_sets(path_sav,path_chem,reactions,species,'10',kept_species = kept)

        # restore species/reaction lists
        chn = []
        
        if strategy == 'da':  # remove aerosol of provided species
            for i in tar: # species
                species[sps.index(i)].condensed = True
                chn.append(i)
                
        elif strategy in ['rs','rsn']: # remove species
            for i in tar:
                species[sps.index(i)].status = True
                chn.append(i)
            # resume reaction list
            for i in out_info[0] + out_info[1]:
                reactions[i[0]].loadSSH(i[1])
            # limit size of products
            if len(out_info[1]) >= 5:
                out_info[1] = ['with {:d} reactions'.format(len(out_info[1]))]
            chn = [chn, out_info]
            
        elif strategy == 'rm': # remove reactions
            for i in tar:
                if i in irel_rcn_pvoc:
                    reactions[i].products = rel_rcn_pvoc[i][0]
                    reactions[i].ratiosPD = rel_rcn_pvoc[i][1]
                else:
                    reactions[i].status = True
                chn.append([reactions[i].toSSH(Type='all'),i]) # str, id
                
        elif strategy == 'rmp': # remove reactions
            for i in range(len(tar)):
                # locate changes
                j,k = tar[i]
                pd,rt = out_info[i]
                # for output
                chn.append([reactions[j].toSSH(Type='all'),j,pd,rt]) # str, id
                
                # resume
                reactions[j].products.insert(k,pd)
                reactions[j].ratiosPD.insert(k,rt)

        elif strategy == 'rm1': # remove partial reactions
            chn0, chn1 = [], []
            # restore old reaction list
            for i in list(sorted(ind_new, reverse=True)): # remove new reactions
                chn1.append(reactions[i].toSSH(Type='all'))
                reactions.pop(i)
            for i in ind0: # add back old reactions
                reactions[i].status = True
                chn0.append(reactions[i].toSSH(Type='all'))
            chn = [chn0, chn1]

        # check reaction list
        if imark == 0:
            paths.append('{:s}/{:s}'.format(path_sav,path_chem))
            changes.append(chn)

    return [strategy,paths,changes]

def search_complex_reduction(reactions, species, strategy, rel_sps, frozen, kept, path_sav, cut=None, Nmax=999, lumptypes=['plain']):

    #ttest0 = time.perf_counter()

    fro = [] # frozen species & relationship
    # update frozen species if needs
    if strategy == 'lp': 
        # kept species not undergo lumping
        if frozen+kept != []: fro = [i for i in frozen+kept]
    elif frozen != []: # for other strategy
        fro = frozen

    Ntry = 0 # total try numbers
    paths,changes,lchn = [],[],[]

    while (Ntry < Nmax):
        # copy
        rc = deepcopy(reactions)
        sp = deepcopy(species)


        if strategy == 'jp':# jumping
            lump, lumppd = reduce_reaction_bySpecies(rc,sp,'jump', 1,
                                                       frozen = fro,
                                                       target = rel_sps,
                                                       cut_order = cut[0])
        elif strategy == 'rp': # replacing inside reaction
            lump, lumppd = reduce_reaction_bySpecies(rc,sp,'replace', 1,
                                                       frozen = fro,
                                                       target = rel_sps,
                                                       cut_order = cut[0])
        elif strategy == 'lp': # lumping
            lump, lumppd, lchn = reduce_reaction_byLump(rc,sp, 'DU_1', refconc_paths, 
                                                        concs_ave, 
                                                        frozen = fro, 
                                                        RefConcRead = RefConcRead, 
                                                        target = rel_sps,
                                                        lumptype = lumptypes,
                                                        rcSpsCut = cut)

        # get add times
        if lump == []: break # stop searching
        #elif len(lump) > 1: raise ValueError('ATP: len lump > 1!',lump,len(lump),lumppdmstrategy)
        elif lchn == []: ntimes = 1
        else: ntimes = len(lchn)

        # generate ssh-aerosol and output infos
        for n in range(ntimes):

            if strategy == 'lp':
                if lump[0] not in fro: fro.append(lump[0])
                chn = [lump, lumppd[n]]
            else:
                if strategy == 'jp': # jump
                    itmp,jtmp = '',''
                    for s in lump: itmp+=(s+',') # jumped species
                    for s in lumppd[0][0]: jtmp += (s+',') # new species
                    ifro = '{:s}->{:s}'.format(itmp[:-1], jtmp[:-1])
                    #print('check ifro jp.', ifro)
                elif strategy == 'rp': # replace
                    ifro = lump[0]+'->'+lumppd[0][0][0]
                if ifro not in fro: fro.append(ifro)
                chn = [lump, ifro]
                
            # change rc,sp
            if n != 0:
                #rc, sp = deepcopy(reactions), deepcopy(species)
                lp_sps, lp_rcn = lchn[n]
                for s in list(lp_sps.keys()): sp[s] = lp_sps[s]
                for s in list(lp_rcn.keys()): rc[s] = lp_rcn[s]
            # sav chem
            if cut: path_chem = '{:s}-{:d}_{:d}'.format(strategy,cut[0],Ntry)
            else: path_chem = '{:s}-{:d}'.format(strategy,Ntry)
            Ntry += 1
            
            # trim scheme
            #rc, sp, trim_sp = trim_scheme(rc, sp)
            #chn.append(trim_sp)
            
            imark = to_SSH_sets(path_sav,path_chem, rc, sp, '10', kept_species = kept)
            
            # record
            if imark == 0:
                paths.append('{:s}/{:s}'.format(path_sav,path_chem))
                changes.append(chn)

    return  [strategy,paths,changes]

def run_sim_with_nml(ind, path_ssh, namelist, path_ref_res, ipvoc=0, file_savorg=False):

    """
        ipvoc: type for initial voc. if 0 then no checking on initial concentration. Diff types are embedded in ssh-aerosol.
    """
    wd = os.getcwd()
    os.chdir(path_ssh)

    # run simulation
    if ipvoc: 
        toto = 'toto/toto{:d}-{:d}'.format(ind,ipvoc)
        os.system('./ssh-aerosol {:s} {:d} 1>{:s}  2>&1'.format(namelist,ipvoc,toto))
    else: 
        toto = 'toto/toto{:d}'.format(ind)
        os.system('./ssh-aerosol {:s} 1>{:s}  2>&1'.format(namelist,toto))

    # copy organic concs. if need
    if file_savorg:
        os.system('cp {:s}/aero/Organics_1.txt {:s}'.format(path_ref_res, file_savorg))

    # analysis errors
    errs = [-1.,-1.]
    with open (toto) as f:
        tmp = f.read()
        # check error
        if 'sshError' in tmp or 'sumXk' in tmp:
            print('Error is found in simulation. message: ',tmp, path_ssh,ind,namelist)
            errs = [1.,1.] # not use this rdc
        else:
            for i in tmp.splitlines():
                if 'err_ref' in i and ':' in i: 
                    errs[0] = float(i.split(':')[1])
                    #max([float(j) for j in i.split(':')[1].split(' ') if isfloat(j)])
                elif 'err_pre' in i and ':' in i: 
                    errs[1] = float(i.split(':')[1])
                #elif 'sumXk' in i:
                #    raise ValueError('sumXk is zero is found.',path_ssh,ind,namelist,i)
    #if ipvoc: print('errs: ',errs, path_ssh, toto)
    os.chdir(wd)

    return errs

def build_ssh_from_chem(path_chem, path_ssh):
    # get id for chem
    ichem = path_chem.split('/')[-1]
    strategy, ind = ichem.split('-')
    # get path for ssh
    path_now = '{:s}/train/{:s}/{:s}'.format(path_ssh, strategy, ind)
    # check if exists
    if not os.path.exists(path_now):
        # copy from ssh
        create_folder('{:s}/train/{:s}'.format(path_ssh, strategy))
        os.system('cp -rf {:s}/ssh {:s}'.format(path_ssh, path_now))
    # copy chem
    cmd = ''
    for i in ['aer.vec','reactions','species','RO2']:
        cmd += 'cp {0:s}/{1:s}.{2:s} {3:s}/src/include/CHEMISTRY/BCARY/BCARY.{2:s};'.format(path_chem,ichem,i,path_now)
    os.system(cmd)
    # compile and check result
    wd = os.getcwd()
    os.chdir(path_now)
    for n in range(6):
        if n == 2: os.system('./clean >tmptoto 2>&1; ./compile >tmptoto 2>&1')
        elif n > 4: raise ValueError('Can not compile. Have tried 5 times.',path_now,path_chem) 
        else: os.system('./compile >tmptoto 2>&1')
        
        # check tmptoto
        with open ('tmptoto','r') as f:
            if 'scons: done building targets.' in f.read(): 
                tag = 0
                break
    os.chdir(wd)
    return path_now

def chem_from_rs(rc,sp,isp,path_chem,path_ssh_chem):
    # copy and build new chem
    rc0,sp0 = deepcopy(rc), deepcopy(sp)
    for l,k in enumerate(sp0):
        if k.name == isp: sp0[l].status = False
    for l,k in enumerate(rc0):
        if isp in k.reactants: rc0[l].status=False
        if isp in k.products:
            n = k.products.index(isp)
            k.products.pop(n)
            k.ratiosPD.pop(n)
    to_SSH_sets(path_ssh_chem+'/order',path_chem,rc0,sp0,'10')
    #pool_sshs.append(('{:s}/{:s}/{:s}'.format(path_ssh_chem,'order',path_chem),pathSSH_rdc_abs))
    return '{:s}/{:s}/{:s}'.format(path_ssh_chem,'order',path_chem)

if __name__ == '__main__':
    auto_training_para()
