# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  AutoTrainingSeries.py runs the series reduction.
#
#  To be checked woth GENOA v2.0. 
#
# ================================================================================

import sys
import os
import json
import shutil
import time
import numpy as np

from string import ascii_letters, digits

from Parameters import RunSets, Roptions, prefix, \
                       primaryVOCs, pathSSH_rdc, \
                       locs, Ttot, DeltaT, Tnow, \
                       initfile, pathInitFiles, \
                       pathNewRes, pathNewChem, tail
from Functions import isfloat, isint, get_info, \
                      create_folder, get_locations
from DataStream import read_chem_sets,to_SSH_sets
from ReductionStrategy import reduce_reaction_byRemove, \
                              reduce_reaction_byLump, \
                              reduce_reaction_bySpecies, \
                              trim_scheme, reaction_seperate, \
                              reaction_merge
from SSHSetUp import run_SSH
from KineticMCMtoPython import update_kinetic_rate
from AutoTesting import auto_testing
from ChemRelation import tree, get_species_from_reactions

# default
auto_type_name = {'rm': 'Removing reaction',
                  'rm1':'Removing elementary-like reaction',
                  'lp': 'Lumping',
                  'da': 'Removing gas-particle partitioning',
                  'rs': 'Removing species',
                  'jp': 'Jumping',
                  'rp': 'Replacing' }

auto_names = list(auto_type_name.keys())
letters = ascii_letters + digits

def auto_training_srs(Setups = RunSets[3], locs = locs, pathInitFiles = pathInitFiles, err_ref_pre = 0.50, ATSetups = RunSets[4]):
    """
        Objectives: run GENOA reduction

        Inputs:
            
        Outputs:
    """

    ## initialization
    # read this part of info from Parameters.py

    IDchemPre = Setups['IDchemPre']
    IDchemVad = IDchemPre
    prePath = os.path.abspath(Setups['prePath'])
    pre_chem_path = Setups['preChemPath']

    IDchemRef = Setups['IDchemRef']
    refPath = os.path.abspath(Setups['refPath'])

    IDchemFake = Setups['IDchemFake']
    fakePath = os.path.abspath(Setups['fakePath'])

    err_refs = Setups['err_ref']
    err_pres = Setups['err_pre']

    auto_types = Setups['strategy_types']
    BranchRatios = Setups['BranchRatio']

    # frozen compounds
    NoLumps = Setups['frozenspecies']

    # tags
    tag_quick_mode = Setups['tag_quick_mode']
    tag_pre_rdc = Setups['tag_pre_reduction']
    tag_check_median = True
    tag_check_average = True
    tag_testing = Setups['tag_testing']

    tag_check_redo = Setups['tag_redo']
    tag_pipe = Setups['tag_hybrid']

    # number of aerosol that indicates late stage of reduction
    naer_lim = Setups['naer_late-stage']

    tag_check_testing = False
    # assign ref criteria
    try_ave_ref = Setups['try_ave_ref']
    try_max_ref = Setups['try_max_ref']
    try_at_err = Setups['try_at_err']

    # assign try_at_err if not assigned
    if try_at_err > 1. : try_at_err = try_ave_ref + 0.01
    try_ave_max_ref = try_ave_ref + 0.005 # ref err for try_ave_ref
    try_ave_all = try_ave_ref
    try_ave_max_all = try_ave_max_ref + 0.01
    try_max_all = try_max_ref + 0.10
    try_max_ind_use = [] # init

    if try_at_err <= max(err_refs): # check error settings
        print('Read errors for pre-testing.',try_ave_ref,try_max_ref,try_at_err)

    # pre locs
    nPreTest = Setups['nPreTest']
    if ATSetups['Test_file'] is not None:
        if not os.path.exists(ATSetups['Test_file']):
            print('Not find Test_file.', ATSetups['Test_file'])
            raise FileNotFoundError(ATSetups['Test_file'])
        print('read nPreTest: ', nPreTest, 'from file ', ATSetups['Test_file'])
        with open (ATSetups['Test_file'],'r') as f:
            prelocs = json.loads(f.read())[0:nPreTest]
    else:
        raise OSError('Not provied Test_file.')

    iauto0 = 0 # initial

    # add pre-reduction
    if tag_pre_rdc:
        for i,j in enumerate(err_refs):
            if j > tag_pre_rdc: 
                tag_pre_rdc = i # set to compare to nerr in order to stop pre-reduction
                print('pre-reduction when err_ref <=', err_refs[tag_pre_rdc], tag_pre_rdc)
                break

    # generate record files
    if iauto0 > 0: frec = open ('{:s}/CaseRdc_{:s}'.format(pathNewChem,tail),'a+') # add
    else: frec = open ('{:s}/CaseRdc_{:s}'.format(pathNewChem,tail),'w+') # renew

    # directory
    # absolute path
    pathSSH_rdc_abs = os.path.abspath(pathSSH_rdc)
    # sav chemistry/record files in pathNewChem
    path_sav_rec = '{0:s}/{1:s}_recs'.format(pathNewChem,tail)
    path_sav_chem = '{0:s}/{1:s}_chems'.format(pathNewChem,tail)
    path_sav_res = '{0:s}/{1:s}'.format(pathNewRes,tail)

    for i in path_sav_rec, path_sav_chem, path_sav_res: create_folder(i)

    # cp pre chem
    if os.path.exists('{:s}/{:s}'.format(pre_chem_path,IDchemPre)):
        if not os.path.exists('{:s}/{:s}'.format(path_sav_chem,IDchemPre)):
            shutil.copytree(pre_chem_path+'/'+IDchemPre,path_sav_chem+'/'+IDchemPre)
    else:
        raise FileNotFoundError('Can not find IDchemPre: {:s}/{:s}'.format(pre_chem_path,IDchemPre))

    # init parameters
    tag_iauto = 1 # if 0 stop
    iauto = iauto0 # reduced times

    ntotal = 0 # init totoal valid run in one turn
    nerr = 0 # index of error
    npres = {} # record the valid number in the previous round
    rdcs = [999,[],[]] # record reductions with ave_err for late stage reduction # naer, [name], [errs]
    tag_rm_out = 0
    tag_rm1 = False # default no rm1
    tag_redo = 0
    tag_stage = 0
    nauto = len(auto_types)
    emax_ind, emax, eave, eave_max = [], 0.0, 0.0, 0.0 # init for testing

    if 'rm' in auto_types: nBRT = 0 # index of branching ratio

    # check input length in case of quick mode
    if tag_quick_mode:
        i, j, k, l = len(err_refs), len(err_pres), len(auto_types), len(BranchRatios)
        if max(i, j, k, l) != min(i, j, k, l):
            print('ATS: in quick mode, check the length of err_refs, err_pres, auto_types, and BranchRatios.',i,j,k,l)
            raise ValueError

    if tag_pipe: # init for pipe reduction
        rcd_pipe = {}
        npipe = 0
        pipe_list = auto_types

    # get coordinate
    locations = get_locations(locs)
    # get concs
    concs_ave, concs_min, refconc_paths, RefConcRead = update_kinetic_rate(locs,
                                                          IDchem_fake = IDchemFake,
                                                          path_fake = fakePath,
                                                          pathInitFiles = pathInitFiles,
                                                          RefConcRead = path_sav_chem+'/s_'+tail)
    # prepare ref cases at different location and different time
    paths_pre, paths_ref, cases  = [], [], []
    n = 0
    for il in locs:
        if initfile in ['month','storage']:
            lc = 'm{:d}y{:d}x{:d}'.format(il[2],il[0],il[1])
        else:
            lc = 'y{:d}x{:d}'.format(il[0],il[1])
        for inow in Tnow:
            for idx in range(len(Ttot)):
                iT=Ttot[idx] # total time
                iD=DeltaT[idx] # time step
                icase = '{:s}/c_{:d}_{:d}_{:d}h'.format(lc,iT,iD,inow)
                cases.append(icase) # for output
                iref = '{:s}/{:s}_{:s}'.format(refPath,icase,IDchemRef.replace(prefix,''))
                if os.path.exists(iref): paths_ref.append(iref)
                else:
                    raise FileNotFoundError('ATS: ref file not found. '+iref)
                ipre = '{:s}/{:s}_{:s}'.format(prePath,icase,IDchemPre.replace(prefix,''))
                if os.path.exists(ipre): 
                    # copy files
                    shutil.copy('{:s}/aero/Organics_1.txt'.format(ipre),'{:s}/{:d}.txt'.format(pathSSH_rdc,n))
                    paths_pre.append('{:d}.txt'.format(n))
                    n += 1
                else:
                    raise FileNotFoundError('ATS: pre file not found. '+ipre)
    ncases = len(cases)

    ## reduction start
    while tag_iauto:
        run_time = time.perf_counter()
        if not tag_redo: # do not need to redo

            # set pre case
            if iauto != iauto0:
                if tag_pipe and npipe: # if pre-recorded
                    IDchemPre = IDchemPre_pipe
                    prePath = prePath_pipe
                    # reset
                    if npipe == -1: npipe = 0
                elif nval != 0: # update pre
                    IDchemPre = IDchem
                    prePath = curPath

            # prepare for the new run: IDchem, IDchemPre, prePath
            if iauto < len(letters): midname = ''
            else: 
                midname = digits[int((iauto - len(letters))/len(letters))]
                #print('Use midname: ', midname, iauto)
            IDchem = prefix + tail + midname + letters[iauto % len(letters)]

            auto_type = auto_types[iauto % nauto] # strategy

            # in quick mode one err/ratio for one type
            if tag_quick_mode: nerr, nBRT = iauto, iauto
            # assign errors
            if nerr < len(err_pres): err_pre = err_pres[nerr]
            else: err_pre = err_pres[-1]

            if nerr < len(err_refs): err_ref = err_refs[nerr]
            else: err_ref = err_refs[-1]

        # check if add pre testing
        if not tag_check_testing:
            if err_ref >= try_at_err: tag_check_testing = True

        # prepare the record files
        recordfile = '{:s}/Record_{:s}_{:s}'.format(path_sav_rec,IDchem,auto_type)
        # output info
        tmp = 'Training. Reduction strategy: {:s}\nIDchemPre: {:s}\tprePath: {:s}\nIDchemRef: {:s}\trefPath: {:s}\nError Tolerance: err_ref <= {:f}, err_pre <= {:f}.\nCurrent reduction step: {:d} with tail: {:s}. Reduction stage: {:d}\n'.format(auto_type_name[auto_type],IDchemPre,prePath,IDchemRef,refPath,err_ref,err_pre,iauto+1,tail,tag_stage)
        # special treatments
        if auto_type == 'rm': 
            BranchRatio = BranchRatios[nBRT]
            tmp += 'Branching ratio for removing reactions: {:f}\tnBRT: {:d}\n'.format(BranchRatio,nBRT)
        if tag_quick_mode: 
            tmp += 'In Quick Mode!\n'
        if tag_pipe and auto_type in pipe_list:
            tmp +=  'In Pipe Mode!\n'

        # reopen for the output
        fall = open (recordfile+'_all','w+')
        fuse = open (recordfile+'_use','w+')
        for f in [fall,fuse,frec]:
            f.write(tmp)
            f.flush()

        # prepare folder and tmp folder
        ResultFolder = 'Results_{:s}'.format(IDchem)
        curPath = '{:s}/{:s}'.format(os.path.abspath(path_sav_res),ResultFolder)
        speciesfile = '{:s}/{:s}/{:s}.mol'.format(path_sav_chem,IDchemPre,IDchemPre)
        reactionfile = '{:s}/{:s}/{:s}.reactions'.format(path_sav_chem,IDchemPre,IDchemPre)
        speciesfile_tmp = '{:s}/{:s}_tmp/{:s}_tmp.mol'.format(path_sav_chem,IDchem,IDchem)
        reactionfile_tmp = '{:s}/{:s}_tmp/{:s}_tmp.reactions'.format(path_sav_chem,IDchem,IDchem)

        # get the targets for certain auto_type
        # read input MCM files and compute the properties of species
        rc,sp=read_chem_sets(reactionfile,speciesfile,speciesType='SSH',reactionType='SSH')

        # organise reaction list for rm1
        if auto_type == 'rm1':
            tag_rm1 = True # active if no rdc with rm
            rc = reaction_seperate(rc,sp)
            BranchRatio = 1.0 # try all reactions

        if iauto == iauto0 or (tag_rm1 and auto_type == 'rm1'):
            nrea, ngas, naer = get_info(rc,sp,'') # record scheme
            frec.write('Initial scheme No.reaction: {:d}\tNo.gas: {:d}\tNo.aerosol: {:d}\n'.format(nrea, ngas, naer))
            if naer <= naer_lim and not tag_stage: tag_stage = 1
            if tag_testing and IDchemPre != IDchemRef: # record the pre-testing result of the input pre case
                ATSetups['err_out'] = 1.
                emax_ind, emax, eave, eave_max = auto_testing(Setups = ATSetups, 
                                                              IDchem = IDchemPre,
                                                              ind = prelocs,
                                                              chempath = path_sav_chem,
                                                              trd_max = 1,
                                                              out_file = False,
                                                              sav_soapath = False)
                if emax_ind != []: emax_ind = 'm{:d}y{:d}x{:d}'.format(emax_ind[0][2], emax_ind[0][0], emax_ind[0][1])
                frec.write('Pre-Testing on IDchemPre {:s}: err loc: {:}\terr max: {:6.4f}\t err ave: {:6.4f}\t err ave max: {:6.4f}\n'.format(IDchemPre, emax_ind, emax, eave, eave_max))
            frec.flush()

        if tag_testing and iauto == iauto0: # check pre case
            if eave > try_ave_all or emax > try_max_all or eave_max > try_ave_max_all: # error exceed limit
                print(iauto,emax_ind,eave,try_ave_all,emax,try_max_all,eave_max,try_ave_max_all)
                raise ValueError('ATS: error in the pre case exceed the limit.')
            else:
                try_ave_now = max(eave, try_ave_ref)
                try_ave_max_now = max(eave_max, try_ave_max_ref)
                try_max_now = max(emax, try_max_ref)

        # copy in tmp folder
        to_SSH_sets(path_sav_chem,IDchem+'_tmp',rc,sp,'10')

        # for those strategies that do not use the functions in ReductionStrategy
        rkaer = []

        if auto_type == 'da':  # remove aerosol

            # obtain aerosol list
            with open ('{:s}/{:s}/{:s}.aer.vec'.format(path_sav_chem,IDchemPre,IDchemPre),'r') as f:
                info_aer = f.read().splitlines()
                tmp = info_aer[10:-1]
                ntmp = len(tmp) #reverse order
                for i in range(ntmp): rkaer.append(tmp[ntmp-i-1].split('\t')[0])

        elif auto_type == 'rs': # remove species

            nsp = len(sp)
            for s in range(nsp):
                i = sp[nsp-s-1]
                if i.name in primaryVOCs or not i.organic: continue
                else: rkaer.append(i.name)

        elif auto_type == 'rm1' and tag_rm1 == False:
            rkaer = []

        # add new treatment for rm, set removing list
        elif tag_stage == 2 and auto_type in ['rm','rm1'] and BranchRatio >= 1.0:
            t1 = time.perf_counter()
            # get species list
            tmp_sp = [i.name for i in sp]
            # get tree list
            tmp_p, tmp_c = tree(rc)
            # get primary list
            tmp0 = {'OH':0,'NO3':0,'O3':0}
            for i,j in enumerate(rc):
                if set(primaryVOCs) & set(j.reactants):
                    for k in ['OH','NO3','O3']:
                        if k in j.reactants: tmp0[k] += 1
            # get the removed species
            for i in range(len(rc)):
                ind = len(rc) - i - 1
                j = rc[ind] # reversed order
                # check primaryVOCs
                if set(primaryVOCs) & set(j.reactants):
                    for k in ['OH','NO3','O3']:
                        if k in j.reactants:
                            tmp = tmp0[k] # get oxidants
                    if tmp <= 1: continue # at least keep one initial reaction
                # get organic products
                tmp = [k for k in j.products if sp[tmp_sp.index(k)].organic]
                # check if it is the only way to get the removed species
                for k in rc:
                    if k == j: continue
                    for s in k.products: # other generates the products
                        if s in tmp: # scheme contains other way to get s
                            for l in k.reactants: # not a loop
                                if l in tmp: 
                                    break
                            tmp.remove(s)
                if tmp != []: # found the removed species
                    for k in tmp:
                        if sp[tmp_sp.index(k)].condensed: # find condensed species
                            rkaer.append(ind)
                            break
                        else: # radical
                            if k not in tmp_p: continue
                            for l in tmp_c[tmp_p.index(k)]: # check if further removed species can condense
                                if sp[tmp_sp.index(l)].condensed:
                                    tag = 0
                                    for s in tmp_c: # check if any other path
                                        if l in tmp_c: tag += 1
                                    if tag == 1: 
                                        rkaer.append(ind) # only one production
                                        break
                                if ind in rkaer: break
                        if ind in rkaer: break
            # record other items
            for i in range(len(rc)):
                j = len(rc) - i - 1
                if j not in rkaer: 
                    # check primary VOC
                    if set(primaryVOCs) & set(rc[j].reactants):
                        for k in ['OH','NO3','O3']:
                            if k in rc[j].reactants: # get oxidants
                                tmp = tmp0[k]
                        if tmp <= 1: continue
                    rkaer.append(j)
        else: rkaer = [-1]


        # init valid attempt times
        nval = 0
        # valid removed reactions
        vals = []
        # reactions that should not be removed
        frozen = []
        
        # add NoLumps in rm
        if auto_type in ['rm','rm1'] and NoLumps != []:
            for i in rc:
                for j in NoLumps:
                    if j in i.reactants or j in i.products: frozen.append(i)
            print('ATS: frozen', NoLumps)

        for i in NoLumps:
            if i in rkaer: 
                print('ATS: frozen species from rkaer',i)
                rkaer.remove(i)
            elif 'P'+i in rkaer:
                print('ATS: frozen species from rkaer','P'+i)
                rkaer.remove('P'+i)
            else: print('ATS: not found species in rkaer',i)

        nrun = 0 # init run times
        while rkaer != []:
            nrun += 1 # run times

            if auto_type == 'lp':
                # reduction by lumping species
                lump,lump_prop,lchn = reduce_reaction_byLump(rc,sp,'DU_1', refconc_paths, concs_ave, 
                                                    frozen = frozen, RefConcRead = RefConcRead)
                # check lumping
                if lump != []:
                    if len(lump)==1: tar = lump[0] # for output
                    else:
                        raise ValueError('lp: lumped multiple times: '+str(len(lump)))
                else: break # stop lumping

            elif auto_type in ['rm','rm1']:
                if rkaer == [-1]:
                    # reduction by removing reactions
                    lump,rcn = reduce_reaction_byRemove(rc,sp,'all_1',BranchRatio,concs_min,frozen)
                else:
                    tmp1 = get_species_from_reactions(rc)
                    tmp_sp = [i.name for i in sp]
                    rcn = rkaer[0]
                    rc[rcn].status = False
                    tmp2 = get_species_from_reactions(rc)
                    for i in tmp1:
                        if i not in tmp2 and sp[tmp_sp.index(i)].organic : 
                            sp[tmp_sp.index(i)].status=False
                    lump = [rc[rcn].toSSH(Type='all')]

                # check removing
                if lump != []:
                    if len(lump)!=1:
                        print('ATS: removing multiple times not only once.', lump, auto_type)
                        raise ValueError
                else: break # stop removing
                # output
                tar = rcn

            elif auto_type == 'jp':
                # reduction by removing reactions
                lump, lumppd = reduce_reaction_bySpecies(rc,sp,'jump',1, frozen)
                if len(lump) == 0: break # no rp found, stop
                else: # for print and frozen
                    tar = ''
                    if len(lump) == 1 : # 'A1 A2 A3 -> B' or 'A1 -> B' 
                        for s in lumppd[0][0]: tar+=(s+',')
                        tar = '{:s}->{:s}'.format(lump[0], tar[:-1])
                    else: # 'A -> B1 B2 B3'
                        for s in lump: tar+=(s+',')
                        tar = '{:s}->{:s}'.format(tar[:-1],lumppd[0][0][0])
                
            elif auto_type == 'rp':
                # reduction by removing reactions
                lump, lumppd = reduce_reaction_bySpecies(rc,sp,'replace',1, frozen)
                if len(lump) == 0: break
                else:
                    tar = lump[0]+'->'+lumppd[0][0][0]

            else: # for those not use functions from ReductionStrategy
                tar = rkaer[0]

                if auto_type == 'rr':
                    for i in sp:
                        if i.name == rkaer[0]: 
                            i.status = False
                    for i in range(len(rc)):
                        # remove products
                        n = 0
                        while n < len(rc[i].products):
                            js = rc[i].products[n]
                            if js == rkaer[0]:
                                rc[i].products.pop(n)
                                rc[i].ratiosPD.pop(n)
                            else: n += 1
                        # remove reaction
                        for j in rc[i].reactants: 
                            if j == rkaer[0]: 
                                rc[i].status = False

                elif auto_type == 'da':
                    tag_da = 0
                    for i in sp:
                        if 'P'+i.name == rkaer[0]: 
                            i.condensed = False
                            tag_da = 1
                            break
                    if not tag_da:
                        print('ATS: not found aerosol that need to be removed: ',rkaer[0],auto_type)
                        raise NameError

                elif auto_type == 'rs':
                    out_info = [[],[]] # as reactants, as products
                    for i in range(len(rc)):
                        if rkaer[0] in rc[i].reactants: 
                            out_info[0].append(rc[i].toSSH())
                            rc[i].status = False 
                        if rkaer[0] in rc[i].products:
                            out_info[1].append(rc[i].toSSH())
                            n = 0
                            while n < len(rc[i].products):
                                js = rc[i].products[n]
                                if js == rkaer[0]:
                                    rc[i].products.pop(n)
                                    rc[i].ratiosPD.pop(n)
                                else: n += 1 

            # clean up
            rc,sp,trim_sp = trim_scheme(rc,sp)
            # output into SSH-aerosol
            to_SSH_sets(path_sav_chem,IDchem,rc,sp,'10')

            # init for simulations
            tag_sim = 1 # 0/1/2 for no train/ train/ train with wider errs
            trd_max = 1 # for pre-testing/testing
            tag_err = 0 # for evaluation # default 0 accepted
            fout = [] # record results

            nrea0, ngas0, naer0 = get_info(rc,sp,'')
            # add special case that not check training dataset
            if tag_stage and tag_check_testing:
                if auto_type in ['lp','rp','jp'] or tag_stage == 2: # late stage
                    if naer0 < naer: tag_sim = 0

            if tag_sim:
                # simulation
                paths_now, err = run_SSH(IDchem,path_sav_chem,pathSSH_rdc,False,
                                        IDchem.replace(prefix,''),ResultFolder,
                                        Ttot,DeltaT,[locations,locs],Tnow,
                                        'fast_compile',initfile,pathInitFiles,
                                        ref_file = paths_ref, pre_file = paths_pre)
                # record
                err0, err1 = [i[0] for i in err], [i[1] for i in err]
                for i in range(ncases):
                    fout.append('Case: {:s}\tSOA_ref: {:6.4f}\tSOA_pre: {:6.4f}\n'.format(cases[i],err0[i],err1[i]))
                # get ave error
                err_ave = np.average(err0)
                # check max error
                if max(err0) > err_ref: tag_err = 1 # refused
                elif max(err1) > err_pre: tag_err = 1 # refused
                # check mediam error
                if not tag_err and tag_check_median and (np.median(err0) > err_ref * 0.5 or np.median(err1) > err_pre * 0.5):
                    tag_err = 1
                # check avergae error
                if not tag_err and tag_check_average and (np.average(err0) > err_ref * 0.5 or np.average(err1) > err_pre * 0.5):
                    tag_err = 1
                    
            # add pre-testing to test
            if tag_check_testing or (tag_check_redo and tag_redo):
                try_testing = 0
                if not tag_err: 
                    ATSetups['err_out'] = try_max_all
                    try_max_ind, try_max, try_ave, try_ave_max = auto_testing(Setups = ATSetups, 
                                                                              IDchem = IDchem,
                                                                              chempath = path_sav_chem,
                                                                              ind = prelocs,
                                                                              trd_max = trd_max,
                                                                              out_file = False,
                                                                              sav_soapath = False)
                    # change format from [[y,x,m]] to str
                    try_max_ind = 'm{:d}y{:d}x{:d}'.format(try_max_ind[0][2], try_max_ind[0][0], try_max_ind[0][1])
                    # init
                    try_testing = 1
                    # evaluate the pre_testing results
                    if try_ave > try_ave_now or try_max > try_max_ref: tag_err = 1
                    elif tag_stage == 2 and try_ave_max > try_ave_max_ref: tag_err = 1
                    if not tag_err: # sav results
                        try_max_ind_use, try_max_use, try_ave_use, try_ave_max_use = try_max_ind, try_max, try_ave, try_ave_max
                        # update try_ave_now
                        if try_ave > try_ave_now:
                            try_ave_now = try_ave
                        elif try_ave_now > try_ave_ref:
                            try_ave_now = max(try_ave_ref, try_ave)
                        # update try_max_now
                        if try_max > try_max_now:
                            try_max_now = try_max
                        elif try_max_now > try_max_ref:
                            try_max_now = max(try_max_ref, try_max)
                        # update try_ave_max_now
                        if try_ave_max > try_ave_max_now:
                            try_ave_max_now = try_ave_max
                        elif try_ave_max_now > try_ave_max_ref:
                            try_ave_max_now = max(try_ave_max_ref,try_ave_max)

                # prepare pre-testing info for recording
                if try_testing: 
                    if tag_check_testing: 
                        if tag_sim == 0:  # no simulate on training conditions
                            try_text = 'NoSimCuzNaer\t'
                        elif tag_sim == 1:
                            try_text = 'Sim\t'
                        else:
                            raise TypeError('No recognize tag_sim: ',tag_sim)
                        # size
                        try_text += '[{:d}, {:d}, {:d}]\t'.format(nrea0, ngas0, naer0)
                    elif tag_redo: try_text = 'Redo\t'
                    # record errors
                    try_text += 'err loc: {:s}\terr max: {:6.4f} <= {:6.4f}\terr ave: {:6.4f} ({:6.4f}) <= {:6.4f}\trefuse: {:d}\n'.format(try_max_ind, try_max, try_max_now, try_ave_max, try_ave, try_ave_all,tag_err)

            if tag_err: # refused 
                if auto_type  in ['rm','rm1'] and rkaer == [-1]:
                    for i in rcn: frozen.append(rc[i])
                else: frozen.append(tar)
                # reload chem
                rc,sp = read_chem_sets(reactionfile_tmp,speciesfile_tmp,speciesType='SSH',reactionType='SSH')
                # for record
                files = [fall]

            else: # accepted
                nval += 1
                # save changes
                if auto_type in ['rm','rm1']:
                    for i in lump: vals.append(i)
                else: vals.append(tar)

                # update paths_pre
                if not tag_sim: # get paths_now if not simulated before
                    paths_now, err = run_SSH(IDchem,path_sav_chem,pathSSH_rdc,False,
                                                IDchem.replace(prefix,''),ResultFolder,
                                                Ttot,DeltaT,[locations,locs],Tnow,
                                                'fast_compile',initfile,pathInitFiles)
                for i,j in enumerate(paths_now):
                    shutil.copy('{:s}/{:s}/aero/Organics_1.txt'.format(pathSSH_rdc,j),
                                '{:s}/{:d}.txt'.format(pathSSH_rdc,i))
                # update chem in tmp folder
                to_SSH_sets(path_sav_chem,IDchem+'_tmp',rc,sp,'10')
                # update size
                #nrea, ngas, naer = nrea0, ngas0, naer0
                # for record
                files = [fall,fuse]

            # record
            if (tag_check_testing or (tag_check_redo and tag_redo)) and try_testing:
                frec.write(try_text)
                frec.flush()
            for f in files:
                tmp = '\n-----------No.{:d}---------No.valid {:d}-------'.format(nrun,nval)
                if tag_sim: tmp += 'Err_ave {:6.4f}'.format(err_ave)
                f.write(tmp+'\n')
                # species
                if auto_type == 'lp':
                    for i in tar: f.write(i+'\t')
                    for i in lump_prop[0]: f.write('\n'+i) # use one lump ratio
                    f.write('\n')

                elif auto_type  in ['rm','rm1']:
                    for i in lump: f.write(i+'\n')
                else:
                    f.write(tar+'\n')

                if auto_type == 'rs':
                    if out_info[0] != []:
                        f.write('As reactants\n')
                        for i in out_info[0]:f.write(i+'\n')
                    if out_info[1] != []:
                        f.write('As products\n')
                        for i in out_info[1]:f.write(i+'\n')

                # record trimmed species if any
                if trim_sp != []:
                    f.write('Trim species:')
                    for i in trim_sp: f.write(' '+i)
                    f.write('\n')

                # results of testing
                if (tag_check_testing or (tag_check_redo and tag_redo)) and try_testing: f.write(try_text)

                # results
                for i in fout: f.write(i)
                # separator
                f.write('------------------------------------------\n')
                f.flush()

            if rkaer != [-1]: 
                if tar in rkaer:
                    rkaer.remove(tar)
                else:
                    print(tar, rkaer)
                    raise NameError('ATS: item not in the remove list')
            
            # special treat for rm: only remove once
            if auto_type in ['rm','rm1']: 
                if tag_stage and nval >= 1: 
                    rkaer = []
                    for f in [fall,fuse,frec]: f.write('Only reduced once.\n')

        # update scheme
        nrea, ngas, naer = get_info(rc,sp,'')

        if auto_type == 'rm1':
            frec.write('--- rm1 scheme: [{:d}, {:d}, {:d}]\n'.format(nrea, ngas, naer))
            rc = reaction_merge(rc,sp) # remerge
            nrea, ngas, naer = get_info(rc,sp,'') # get size

        if naer < 20: to_SSH_sets(path_sav_chem,IDchem,rc,sp,'20') #output viz file
        else: to_SSH_sets(path_sav_chem,IDchem,rc,sp,'10')

        # final print out
        tmp = 'END, total run: {:d} times, valid run: {:d}\n'.format(nrun,nval)
        tmp += 'Current scheme No.reaction: {:d}\tNo.gas: {:d}\tNo.aerosol: {:d}\n'.format(nrea, ngas, naer)
        tmp += '------------------------------------------\n'
        if nval: tmp += 'Valid:\n'
        else: tmp += 'No vaild reduction (nval = 0).\n'

        for f in [fall,fuse,frec]:
            f.write(tmp)

            if auto_type == 'lp':
                for i in vals:
                    for j in i: f.write('\t'+j)
                    f.write('\n') #'\t|\t')
            else:
                for i in vals: f.write(i+'\n')
            f.flush()

        # pre_testing and record if need
        if tag_testing or tag_pipe:
            if nval or (tag_check_redo and tag_redo):
                # check if pre exceed
                #if eave_max > err_ref and not tag_redo

                if (tag_check_testing or tag_redo) and try_max_ind_use != []:
                    emax_ind, emax, eave, eave_max = try_max_ind_use, try_max_use, try_ave_use, try_ave_max_use
                else:
                    ATSetups['err_out'] = 1.0 # no out in pre-testing
                    emax_ind, emax, eave, eave_max = auto_testing(Setups = ATSetups, 
                                                                  IDchem = IDchem,
                                                                  chempath = path_sav_chem,
                                                                  ind = prelocs,
                                                                  trd_max = 1,
                                                                  out_file = False,
                                                                  sav_soapath = False)
                    # change format from [[y,x,m]] to str
                    if emax_ind != []: emax_ind = 'm{:d}y{:d}x{:d}'.format(emax_ind[0][2], emax_ind[0][0], emax_ind[0][1])

                # check the results
                if eave_max > err_ref:  #eave
                    if tag_check_redo and eave_max > try_ave_all:
                        if tag_redo == 1:
                            tag_redo += 1
                            raise RuntimeError('ATS: need to redo twice.') # should not happed
                        else: 
                            tag_redo = 1 # if dont want to redo, comment this line
                    else: # not check redo by pre-testing results
                        tag_redo = 0
                else: tag_redo = 0 # comment if not want to active tag_redo

            # update stage
            if naer <= naer_lim: 
                if not tag_stage: tag_stage = 1
                if not tag_pipe:
                    if rdcs[0] > naer: # reset
                        rdcs = [naer,[],[]]
                    elif rdcs[0] < naer:
                        raise ValueError('naer > rdcs[0], check values.', rdcs, naer)
                    # update 
                    rdcs[1].append([IDchem,curPath]) # name and path
                    rdcs[2].append(eave) # error

            for f in [fall,fuse,frec]:
                f.write('Pre-Testing IDchem {:s}: err loc: {:s}\terr max: {:6.4f}\t err ave: {:6.4f}\t err ave max: {:6.4f}\n'.format(IDchem, emax_ind, emax, eave, eave_max))
                f.flush()

        # save results
        tag_rm_out = 0

        # record pipe info
        if tag_pipe and auto_type in pipe_list:
            if rcd_pipe == {}: # save pre case when first 
                IDchemPre_pipe = IDchemPre
                prePath_pipe = prePath

            rcd_pipe[auto_type] = [nval, [nrea, ngas, naer], [emax, eave], [IDchem,curPath]] # NO.rdc, scheme, error
            if nval: npipe += 1

        if not tag_redo and nval:
            # obtain new results and put into path_sav_res
            paths_now, err = run_SSH(IDchem,path_sav_chem,pathSSH_rdc,path_sav_res,
                                        IDchem.replace(prefix,''),ResultFolder,
                                        Ttot,DeltaT,[locations,locs],Tnow,'fast_compile',
                                        initfile,pathInitFiles)
            # update valid IDchem
            IDchemVad = IDchem
            for i,j in enumerate(paths_now):
                shutil.copy('{:s}/{:s}/aero/Organics_1.txt'.format(path_sav_res,j),
                            '{:s}/{:d}.txt'.format(pathSSH_rdc,i))
        elif nval == 0:
            # special treatment for 'rm'
            if auto_type == 'rm' and not tag_rm_out and len(npres) == nauto and sum(npres.values()) == npres['rm']:
                if 'rm1' not in auto_types: tag_rm_out = 1
            if auto_type == 'rm1' and not tag_rm_out and len(npres) == nauto and sum(npres.values()) == npres['rm1']:
                tag_rm_out = 1

        # remove previous tmp files
        shutil.rmtree('{:s}/{:s}'.format(path_sav_chem,IDchem+'_tmp'),ignore_errors=True)

        if not tag_redo:

            # record previous round info
            npres[auto_type] = nval

            # check tag_iauto
            if iauto % nauto == 0: ntotal = nval # start new round
            else: ntotal += nval # compute total reduced times

            # for pre-reduction
            if tag_pre_rdc and (iauto + 1) % nauto == 0: ntotal = 0

            # if pipe finish one round
            if tag_pipe and sorted(pipe_list) == sorted(list(rcd_pipe.keys())):
                # assign IDchem_pipe
                # find more than one reduction
                ipipe = None # init
                if npipe > 1:
                    tmp_tag, k = 1, 2
                    while tmp_tag and k > -1:
                        tmp_val, tmp_id, tmp_err = [], [], []
                        for i,j in enumerate(pipe_list):
                            # check
                            if rcd_pipe[j][1][k] in tmp_val:
                                f = tmp_val.index(rcd_pipe[j][1][k])
                                tmp_id[f].append(j)
                                tmp_err[f].append(rcd_pipe[j][2][1])
                            else:
                                tmp_val.append(rcd_pipe[j][1][k])
                                tmp_id.append([j])
                                tmp_err.append([rcd_pipe[j][2][1]])
                        # check values
                        if len(tmp_val) == 1: # all the same
                            k -= 1 # try next
                        else:
                            i = tmp_val.index(min(tmp_val))
                            if len(tmp_id[i]) > 1: # check by errors
                                j = tmp_err[i].index(min(tmp_err[i]))
                            else:
                                j = 0
                            ipipe = tmp_id[i][j]
                            tmp_tag = 0
                elif npipe == 1: # only one reduction
                    for i,j in enumerate(pipe_list):
                        if rcd_pipe[j][0]:
                            ipipe = j
                            break

                if ipipe is not None:
                    IDchemPre_pipe = rcd_pipe[ipipe][3][0]
                    prePath_pipe = rcd_pipe[ipipe][3][1]

                    # record
                    tmp = '!!! Need pipe selection: {:d}\n'.format(npipe)
                    tmp += json.dumps(rcd_pipe)+'\n'
                    tmp += 'Found! type: {:s}, pipeID: {:s}, pipePath: {:s}\n'.format(ipipe,IDchemPre_pipe,prePath_pipe)
                    frec.write(tmp)
                    # reset
                    npipe = -1
                    # update scheme
                    nrea, ngas, naer = rcd_pipe[ipipe][1]
                    frec.write('Pipe scheme {:}\n'.format(rcd_pipe[ipipe][1]))

                else: frec.write('!!! No pipe selection: \n')

                rcd_pipe = {} # reset

            # check if stop
            if tag_quick_mode:

                if iauto == nauto - 1: tag_iauto = 0 # only one round in quick mode

            elif ((iauto + 1) % nauto == 0 and ntotal == 0) or tag_rm_out: # reach the end of 1 round
                # 1st update branching ratio
                if 'rm' in auto_types and nBRT < (len(BranchRatios) - 1) :
                    nBRT += 1
                    npres = {} # reset after 1 round
                # 2nd check err index
                elif nerr < (len(err_refs) - 1):
                    if naer <= 1: # termination criteria defined by user # change later
                        tag_iauto = 0
                    else:
                        nerr += 1 # previous run no valid
                        npres = {} # reset after 1 round
                        if 'rm' in auto_types: 
                            nBRT = 0
                            if tag_stage: # set for late stage
                                if BranchRatios != [1.]: BranchRatios = [1.]
                            # set only for BCARY reduction
                            elif len(BranchRatios) >= 3 and BranchRatios != [1E-1, 5E-1, 1.] and err_refs[nerr] == 0.03 : 
                                BranchRatios = [1E-1, 5E-1, 1.]

                # meet the condition to stop
                else:
                    tag_iauto = 0 # stop reduction if in 1 round no reduced found

                # update auto_types
                if tag_iauto == 0 and 'rm1' not in auto_types:
                    tag_iauto = 1
                    tag_rm1 = 1

            elif (iauto + 1) >= len(letters) * (len(digits) + 1): # reach the maximum times (682)
                tag_iauto = 0

        # record
        run_time = time.perf_counter() - run_time

        for f in [fall,fuse,frec]:
            #f.write('No.reduced\t{:d}\tNo.round\t{:d}\tNo.ValidTillNow - {:d}\t{:d}\tnerr_next\t{:d}\ttime: {:.1f}\n'.format(iauto+1, int(iauto/nauto)+1, (iauto+1)%nauto, ntotal, nerr, run_time))
            f.write('No.reduced\t{:d}\tNo.round\t{:d}\tUsed time: {:.1f}\n'.format(iauto+1, int(iauto/nauto)+1, run_time))
            if tag_pre_rdc: f.write('!!!Pre-reduction process.\n')
            if tag_rm_out: f.write('!!!tag_rm_out active. Jump to the next set of error tolerance.\n')
            if tag_redo: 
                if tag_check_redo: f.write('!!!tag_redo active!!! cuz ave_err too large.\n')
                else: f.write('!!!tag_redo active!!! try not set limitation on lump Psat.\n')
            if tag_iauto == 0: f.write('!!!Stop running. The final reduced mechanism is {:s}.\n'.format(IDchemVad))
            f.write('===============================================\n')

        for f in [fall,fuse]: f.close()

        # set only for pre-reduction without lp
        if tag_pre_rdc and 'lp' not in auto_types and nerr == tag_pre_rdc: # pre-reduction finished
            auto_types = ['rm','jp','lp','rp','rs','da']
            nauto = len(auto_types)
            # change iauto from
            iauto = (int(iauto/nauto)+1) * nauto
            print('Finish pre-reduction. current iauto: ',iauto)
            tag_pre_rdc = False
        else:
            if not tag_rm_out and not tag_redo: iauto += 1 # record reduced time
            elif tag_iauto != 0: # backup
                if tag_redo: isp = 'redobk'
                elif tag_rm_out: isp = 'rmbk'
                shutil.move('{:s}/{:s}'.format(path_sav_chem,IDchem),
                            '{:s}/{:s}_{:s}'.format(path_sav_chem,IDchem,isp))
                for f in ['all','use']: 
                    shutil.move(recordfile+'_'+f,
                                '{:s}_{:s}'.format(recordfile+'_'+f,isp))

        # update reduction strategy
        if tag_rm1 and 'rm1' not in auto_types and 'rm' in auto_types:
            auto_types.insert(auto_types.index('rm')+1,'rm1')
            #auto_types = ['rm','rm1','jp','lp','rp','rs','da']
            nauto = len(auto_types) # reser nauto
            BranchRatios = [1.]
            nBRT = 0
            if iauto%nauto: # reset iauto
                iauto += nauto - iauto%nauto
            #frec.write('+=+= change auto_types: add rm1\n')
            frec.write('+=+= change auto_types: add rm1. Now: {:}\n'.format(auto_types))
            # update stage
            tag_stage = 2
            # set pre case from rdcs cases
            if len(rdcs[1]) > 1:
                err = min(rdcs[2])
                i = rdcs[2].index(err)
                IDchem,curPath = rdcs[1][i][0],rdcs[1][i][1]
                nval = 1 # for reset IDchem
                frec.write('Select pre cases {:s} with err {:6.4f} from the list: {:}\n'.format(IDchem,err,rdcs))
    frec.close()

    return IDchemVad, path_sav_chem # final IDchem


if __name__ == '__main__':

    auto_training_srs()
