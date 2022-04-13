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

from Parameters import SSHSpeciesInit,NokeepSp,cst,AeroDict,prefix
from Functions import isfloat,compare,isint
from ChemRelation import get_species_from_reactions, generation
from KineticMCMtoSSH import KineticSSHtoStr, mcm_photolysis
import Module as md

# if add fake species for radicals
tag_First = False # first time from MCM to SSH, no input SOAPorg file
tag_letter = False # output viz with letter if possible
tag_seperate = False # seperate reactions

def read_chem_sets(reactionfile,speciesfile,speciesType='MCM',reactionType='MCM'):
    """read MCM reactions (merged if kinetic is the same) and return reactions (Reaction Objects) and species (string) lists
    available species input file Type: MCM, SSH
    available reaction input file Type: MCM
"""
    """
        Objectives:

        Inputs:
            
        Outputs:
    """
    # get reaction list

    # record species appear in the reaction list
    sps = []
    # store info as md.Reaction() instances 
    reactions=[]

    if reactionType == 'MCM': 
        # process MCM reaction in html format
        with open (reactionfile) as f: fin=f.read().splitlines()

        for info in fin:

            if '\xe2\x86\x92' not in info: continue #'\xe2\x86\x92'

            # "Goto MCM\tBCARY + NO3\t\xe2\x86\x92\tNBCO2\t1.90D-11"
            # a_info looks like ['Goto MCM', 'BCARY + NO3', '\xe2\x86\x92', 'NBCO2', '1.90D-11']
            a_info=[i.replace(' ','') for i in info.split('\t')]

            # get species index in a_info
            if len(a_info)==5: #'Goto MCM' exists, index is from 1 
                ind=1
            else:
                ind=0

            # reactants
            reactants=[i.replace(' ','') for i in a_info[ind].split('+')]
            # MMM: some photolysis reactions: change sequence if the first reactant is a NokeepSp species
            # may can be removed
            if reactants[0] in NokeepSp:
                reactants=reactants[1:]+[reactants[0]]
                #print reactants

            # products
            products=[i.replace(' ','') for i in a_info[ind+2].split('+')]

            # raw kinetic rate
            raw_rate=a_info[-1]
            # change D notion to E notion
            for i in '-','+':
                if 'D'+i in raw_rate: raw_rate=raw_rate.replace('D'+i,'E'+i)

            # MMM: check if there is a branch ratio in raw rate
            # if not single product in the reaction
            if tag_seperate: nb = 1
            else:
                a_rt = raw_rate.split('*0.') # get rate and ratio
                nb=len(a_rt)

            # check the part before *0.
            if nb==1: # no ratio, only one product
                ratio=1
                rate=raw_rate
            elif nb>2:
                print(a_rt)
                raise ValueError('DS: pb1 multiple *0. in rate: '+raw_rate)
            else:
                # check the part after *0.
                # '1D-16*0.35'
                if isfloat('0.'+a_rt[1]):
                    ratio='0.'+a_rt[1]
                # '1D-16*0.35*KDEC'
                elif isfloat('0.'+a_rt[1].split('*')[0]):
                    ratio='0.'+a_rt[1].split('*')[0]
                # other cases
                else:
                    print(a_rt)
                    raise ValueError('DS: pb2 *0. not for ratio '+raw_rate)
                # get rate without ratio
                rate=raw_rate.replace('*'+ratio,'')
                ratio=float(ratio)

            ## record chemical reaction
            if not tag_seperate and ind == 0 and len(reactions[-1]) > 1 and reactions[-1].rate.str==rate and compare(reactions[-1].reactants,reactants):
                rcn=reactions[-1] # merge if rate, reactants are same with the previous one 
            else:
                # build new reaction
                reactions.append(md.Reaction())
                rcn=reactions[-1]
                # add reactants and rate
                for rc in reactants: rcn.record_reactant(rc,1)
                rcn.rate.update(rate)

            # add new products and ratio
            for pd in products: rcn.record_product(pd,ratio)

            # add new species in sps
            for i in reactants,products:
                for j in i:
                    if j not in sps: sps.append(j)

        # in case: NBCO2 + HO2 -> NBCOOH KRO2HO2*0.975
        for i in reactions:
            if 'HO2' in i.reactants and len(i.products)== 1 and i.ratiosPD[0] != 1.0:
                i.rate.update(i.rate.str+'*'+str(i.ratiosPD[0]))
                i.ratiosPD[0] = 1.0

    elif reactionType == 'SSH':
        with open (reactionfile) as f: fcn = f.read().splitlines()[2:-1] # remove tabulation lines

        rcn = None # init
        # read reactions
        for i in range(len(fcn)):
            line = fcn[i]

            # remove comment line
            if '===' in line: continue

            # reactants and products
            elif '->' in line:
                # check
                if line[0] == ' ':
                    print(line)
                    raise ValueError('DS: has -> and empty in the beginning for the line.')

                # add a new reaction: '->' is the symbol to add new reaction
                reactions.append(md.Reaction())
                rcn=reactions[-1]

                # get all items
                items = [i for i in line.split(' ') if i != '']

                for j in items:
                    if j == '+' or j == '': continue
                    elif '->' == j :
                        # n index of products
                        n = items.index(j) 
                        break
                    # record reactants
                    else:
                        rcn.reactants.append(j)
                        rcn.ratiosRC.append(1.0)

                # add products
                # branch ratio
                tmp = 1.0
                for j in items[n+1:]: # check: n might be referenced before assignment
                    if j == '+' or j == '//': continue
                    elif isfloat(j):
                        tmp = float(j)
                    else:
                        rcn.products.append(j)
                        rcn.ratiosPD.append(tmp)
                        tmp = 1.0
            # products
            elif line[0] == ' ':
                # rest of products
                items = [i for i in line.split(' ') if i != '']
                tmp = 1.0
                for j in items:
                    if j == '+' or j == '//': continue
                    elif isfloat(j): tmp = float(j)
                    else:
                        # record products
                        rcn.products.append(j)
                        rcn.ratiosPD.append(tmp)
                        tmp = 1.0

            # kinetic rate
            elif line[0] == '%': #!!! may need other symbol
                # kinetic rate in MCM format
                rcn.rate.str = line[1:]
                if 'J' in line: rcn.rate.Photolysis = True
            # kinetic rate
            elif 'KINETIC' in line:
                # kinetic rate in SSH format
                rcn.rate.SSH = line
                if 'PHOTOLYSIS' in line:  rcn.rate.Photolysis = True

                if rcn.rate.str == None: # generate str base on SSH
                    rcn.rate.str = KineticSSHtoStr(rcn.rate.SSH)

        # check
        for rcn in reactions:
            if not rcn.status: continue
            # check previous reaction before creating a new reactions
            for j in rcn.reactants,rcn.ratiosRC: #rcn.products,rcn.reactants,rcn.ratiosPD,rcn.ratiosRC:
                if j == []:
                    raise ValueError('DS: empty in transformation')
            if len(rcn.products) != len(rcn.ratiosPD) or len(rcn.reactants) != len(rcn.ratiosRC):
                print(len(rcn.products),len(rcn.ratiosPD),len(rcn.reactants),len(rcn.ratiosRC))
                raise AttributeError('DS: length not fit')
            if rcn.rate.str is None: 
                print(rcn.reactants, rcn.products, rcn.rate.SSH)
                raise ValueError('DS: rate.str is empty')
            # record species in the sps list
            for i in rcn.products,rcn.reactants:
                for j in i:
                    if j not in sps: sps.append(j)

    elif reactionType == 'KPP':
        # process MCM reaction/species in KPP format PRAM
        with open (reactionfile) as f: fin=f.read().splitlines()

        rates = [] # restore rates
        for info in fin:

            if '{' in info: # find reactions
                # {74.}     C10H15O2O2 + HO2 = C10H16O4iso1 :             KRO2HO2 ;
                # reactants
                reactants, products = [], []
                for i in info.split('=')[0].split('}')[1].split('+'): reactants.append(i.replace(' ','').replace('\t',''))

                # products
                for i in info.split('=')[1].split(':')[0].split('+'): products.append(i.replace(' ','').replace('\t',''))

                # raw kinetic rate
                raw_rate=info.split(':')[1].split(';')[0].replace(' ','')

                # check kinetic rate
                if raw_rate == '0.0': continue # not record cuz rate == 0.0

                # change D notion to E notion
                if 'D' in raw_rate and 'KDEC' not in raw_rate: raw_rate=raw_rate.replace('D','E')

                # sepearte numbers and other info
                rate = ''
                ratio = 1.
                for i in raw_rate.split('*'):
                    if '0.' in i and isfloat(i): ratio *= float(i) # ratio
                    else: rate += i +'*'
                rate = rate[:-1] # remove the last '*'

                # check existed reaction if merge reactions
                ind = -1
                if not tag_seperate: 
                    for i in range(len(reactions)):
                        if compare(reactions[i].reactants,reactants) and rates[i][0]==rate:
                            ind = i # merge if rate, reactants are same with the previous one 
                            break

                # build new reaction
                if ind == -1: 
                    reactions.append(md.Reaction())
                    rates.append([rate,[ratio]])
                    rcn=reactions[-1]

                    # add reactants and rate
                    for rc in reactants: rcn.record_reactant(rc,1)

                    rcn.rate.update(raw_rate) # update the raw rate

                    # add new products and ratio
                    for pd in products: rcn.record_product(pd,1)

                else:
                    # merge into existed reactions
                    rcn=reactions[ind]
                    old_ratio = sum(rates[i][1])

                    rates[i][1].append(ratio)
                    new_ratio = sum(rates[i][1])

                    if round(new_ratio,5) == 1.0: new_rate = rates[i][0]
                    else: new_rate = '{:s}*{:6.4f}'.format(rates[i][0], new_ratio)
                    # update rate
                    rcn.rate.update(new_rate)

                    # change old branching ratio
                    for pd in range(len(rcn.products)): rcn.ratiosPD[pd] *= (old_ratio / new_ratio) 
                    # add new branching ratio
                    for pd in products:
                        if pd in rcn.products: rcn.ratiosPD[rcn.products.index(pd)] += ratio/new_ratio
                        else: rcn.record_product(pd, ratio/new_ratio)

                # add new species in sps
                for i in reactants,products:
                    for j in i:
                        if j not in sps: sps.append(j)

    else:
        raise TypeError('Unknown type in DS:'+reactionType)

    # species list
    if speciesType == 'MCM': 
        soapfile = './files/SOAP_'+speciesfile.split('/')[-1]
        if not os.path.exists(soapfile):
            print('DS: not find soapfile: ', soapfile)
            soapfile = './files/SOAP_'+prefix
        species = sps_to_species(sps,speciesfile,'MCM',soapfile)
    elif speciesType == 'SSH': species = read_species(sps,speciesfile)
    elif speciesType == 'KPP':
        # for current scheme
        species0 = sps_to_species(sps,speciesfile,'KPP')
        sps0 = []
        for i in species0: sps0.append(i.name)
        sps1 = list(set(sps)-set(sps0))
        species1 = sps_to_species(sps1,'../../PRAM/mcm_monoterpene_mass.txt','MCM')

        for i in species1: sps1.append(i.name)
        species = species0 + species1

    else:
        raise TypeError('DS: speciesType unknown '+speciesType)

    return reactions,species

def sps_to_species(sps,speciesfile,Type = 'MCM', soapfile = './files/SOAP_'+prefix):
    """
        Objectives:

        Inputs:
            
        Outputs:
    """

    species=[] # elements are md.Species() instances

    if Type == 'MCM':
        # init using data from Constant.py
        for i in list(SSHSpeciesInit.keys()):
            species.append(md.Species(i))
            species[-1].mass=SSHSpeciesInit[i]

        # process species properties
        spLists=[[],[],[],[]]
        with open (speciesfile) as f: fsp=f.read().splitlines()
        for i in fsp:
            if 'InChI' not in i: continue

            # get a list of [species_name, SMILES, InChI, MWs]
            sls=i.split('\t')
            slist=[]
            for j in sls:
                if j.strip()!='':slist.append(j.strip())
            if len(slist) !=4:
                print(soapfile)
                raise OSError('DS: not 4 split by \\t: '+i)
            for j in range(4): spLists[j].append(slist[j])

        # record in species
        for i in sps:
            if i in list(SSHSpeciesInit.keys()): continue
            # record/process info if name is a match
            if i in spLists[0]:
                species.append(md.Species(i))
                j = species[-1]
                k = spLists[0].index(i)
                j.SMILES=spLists[1][k]
                j.InChI=spLists[2][k]
                j.mass=float(spLists[3][k])
                j.organic=True
                j.source.append('MCM')
                j.update_advanced(AeroDict)

            else: # record not organic species
                species.append(md.Species(i))

            # read gamma from SSH output file
            #read_properties_SSH(species,'./files/gammainfSSH')
        # read SOAP structure
        if 'CRI' in speciesfile or 'cri' in speciesfile : update_AIOMFACformat(species,'./files/SOAP_CRI')
        elif not tag_First: update_AIOMFACformat(species, soapfile)

    elif Type == 'KPP':
        with open (speciesfile) as f: fsp=f.read().splitlines()[1:]
        n0, n1 = 0, 0
        for i in fsp:
            # get a list of [species_name, No. x 4, MWs, functional grounps x 29, log10(p0 (Pa)) at 298 K SIMPOL,log10(H (M/atm)) SIMPOL + AIOMFAC]
            sls=i.split(',')
            sname = sls[0].replace(' ','') # species name
            if sname in sps:
                if sname in list(SSHSpeciesInit.keys()): continue
                n1 += 1
                # record/process info if name is a match
                species.append(md.Species(sname))
                j = species[-1]
                j.formula='C{:d}H{:d}N{:d}O{:d}'.format(int(sls[1]),int(sls[2]),int(sls[3]),int(sls[4]))
                j.mass=float(sls[5]) # MWs
                j.organic=True

                # functional groups
                for n,s in enumerate(['C','H','N','O']):
                    j.functionalGroups[s] = int(sls[n+1])

                # ratios
                tmp = j.functionalGroups
                j.ratios = {'OM/OC':round(j.mass/(tmp['C']*12.),3),
                    'H/C':round(tmp['H']*1./tmp['C'],3),
                    'O/C':round(tmp['O']*1./tmp['C'],3),
                    'N/C':round(tmp['N']*1./tmp['C'],3)}
                j.DU = int(0.5*(2+2*tmp['C']-tmp['H']+tmp['N']))

                # init j.RO2,j.condensed,j.Psat_atm,j.Psat_torr
                if sname[-1] == 'O': j.Radical=True
                elif sname[-2:] == 'O2':
                    j.RO2=True
                    j.Radical=True
                else: j.condensed = True

                j.psat_atm = 10 ** float(sls[-2]) / 101325.0 # pa to atm
                j.psat_torr = j.psat_atm * 760.
                #j.dHvap_KJ = 0.0
                j.gamma= None
                j.henry= 10 ** float(sls[-1])
                j.source.append('PRAM')
            else: 
                n0 +=1
                #print('Not used PRAM species: ', sname)

        print('read species, in total ',n0+n1, 'use: ',n1, 'not use: ',n0)
        # read structure
        update_AIOMFACformat(species, filename = './files/AIMOFAC_PRAM_SOAP.csv',Type = 'Table')

    else:
        raise TypeError('DS: sps_to_species. Not recognized Type')

    return species

def to_latex_sets(path,chem,reactions,species,soapfile = None):
    """
       return a reaction list and a species list in latex format
       Do not output inorganic products/
       inputs:
        chem: chem name
        path: save path
    """    
    # get path
    path = '{:s}/{:s}/'.format(path,chem)
    os.makedirs(path, exist_ok=True)

    # get species list
    sps = []
    for i in species: sps.append(i.name)

    if soapfile is not None:
        print('read new properites from soapfile: ',soapfile)
        with open (soapfile,'r') as f: info = f.read().splitlines()
        for i in info:
            if 'Ginf' in i: 
                tmp = i.split(' ')# add gamma/henry
                print(tmp)
                if tmp[1] in sps:
                    s = species[sps.index(tmp[1])]
                    print('Update: ',s.name, s.henry, s.gamma,float(tmp[-1]), float(tmp[-2]) )
                    s.henry, s.gamma = float(tmp[-1]), float(tmp[-2])
                else: print('to_latex_sets: ', tmp[1], 'not found in species list')
    # seperator
    spr = '\\middlehline\n'

    with open (path+chem+'.latex','w') as f:
        # reactions
        titles = 'No. & Reactions & Kinetic Rate constant$^{a}$\\\\ \n'
        f.write(titles)
        n = 0
        for rcn in reactions:
            if not rcn.status: continue
            n += 1
            text = '{:d}&'.format(n)
            # reactants
            for i,j in enumerate(rcn.reactants):
                if species[sps.index(j)].lump != []: tsp = 'm'+j
                else: tsp = j
                if rcn.ratiosRC[i] == 1.0: text += tsp
                else:  text += '{:5.3E} {:s}'.format(rcn.ratiosRC[i],tsp)
                text += ' + '
            text = text[0:-3] + ' -> '# add ->

            # products
            if rcn.products != []:
                for i,j in enumerate(rcn.products):
                    # check organics
                    if species[sps.index(j)].organic:
                        if species[sps.index(j)].lump != []: tsp = 'm'+j
                        else: tsp = j
                        if rcn.ratiosPD[i] == 1.0: text += tsp
                        else:  text += '{:5.3f} {:s}'.format(rcn.ratiosPD[i],tsp)
                        text += ' + '
                if text[-2] == '+': text = text[0:-3] + ' & '
                else: text += ' & '

            # rate
            rts = rcn.rate.SSH.split(' ')

            if 'TB' in rts: text += '[{:s}] * '.format(rts[rts.index('TB')+1]) # add TB species
            if 'ARR1' in rts: text += rts[rts.index('ARR1')+1]
            elif 'ARR2' in rts: 
                # get the negative value
                if rts[-1][0] == '-': val = rts[-1][1:]
                else: val = '-' + rts[-1]
                text += '{:s} * EXP({:s}/TEMP)'.format(rts[rts.index('ARR2')+1],val.replace('.00',''))
            elif 'ARR3' in rts:
                # get the negative value
                if rts[-1][0] == '-': val = rts[-1][1:]
                else: val = '-' + rts[-1]
                text += '{:s}*TEMP**{:s}*EXP({:s}/TEMP)'.format(rts[rts.index('ARR3')+1],rts[rts.index('ARR3')+2],val.replace('.00',''))
            elif 'MCM3' in rts: text += 'J({:s}, {:s}, {:s})'.format(rts[-3],rts[-2],rts[-1])#photolysis
            else: text += '!!!'

            f.write(text + '\\\\' + spr)
        f.write('\n\n\n')
        # species
        #titles = 'Surrogate & Type & Molecular formula & MW$^a$ & P$_{sat}^b$ & $\Delta$H$_{vap}^c$ & K$_p^d$ & H$^e$ & $\gamma^f$\\\\ \n'
        for s in species:
            if s.organic:
                text = ''
                # name
                if s.lump == []: name = s.name
                else: name = 'm' + s.name
                # formula
                formula = ''
                for i in ['C','H','N','O']: 
                    if i in s.functionalGroups:  
                        if round(s.functionalGroups[i],2) == int(s.functionalGroups[i]):
                            formula+='{:s}$_{{{:.0f}}}$'.format(i, s.functionalGroups[i])
                        else: formula+='{:s}$_{{{:.2f}}}$'.format(i, s.functionalGroups[i])
                if s.condensed:
                    if s.psat_atm <= 1E-13: stype = 'LVOC'
                    else: stype = 'SVOC'
                    text += '{:s}&{:s}&{:s}&{:6.1f}&{:5.2E}&{:5.2E}&{:5.2E}&{:5.2E}'.format(name,stype,
                             formula,s.mass,s.psat_atm,s.dHvap_KJ,s.gamma,s.henry)
                else:
                    if s.Radical: stype = 'Radical'
                    else: stype = 'VOC'
                    text += '{:s}&{:s}&{:s}&{:6.1f}&&&&'.format(name,stype,formula,s.mass) 
                f.write(text + '\\\\' + spr)

def to_SSH_sets(path,chem,reactions,species,tag_complete = 1,tag_fake = False):
    """return files that need in SSH"""
    """
        Objectives:

        Inputs:
            
        Outputs:
    """
    ### files need in SSH
    # filename.reactions rout
    # filename.species sps

    # photolysis-filename pys
    # filename.RO2 RO2
    # filename.cst

    # species-list-gas-filename: from species.spack.dat <- generate by SSH
    # species-list-aerosol-filename

    # namelist.ssh
    # input data
    
    # if tag_complete: only output essential files for SSH-aerosol
    # not output files: Parameter.sav, Reduction.sav, ReductionStrategy.sav, BCARYm3.soap, BCARYm3.gen, BCARYm3.mol, BCARYm3.aer

    # if add fake radicals in reactions
    if tag_fake: # symbol = 'FA'
        sps=[]
        for i in species: sps.append(i.name)
        for n in range(len(reactions)):
            i=0
            while i < len(reactions[n].products):
                if reactions[n].products[i] in sps:
                    if species[sps.index(reactions[n].products[i])].Radical:
                        # add fake species
                        reactions[n].products.append('FA'+reactions[n].products[i])
                        reactions[n].ratiosPD.append(reactions[n].ratiosPD[i])
                i+=1

    # default files for init aerosol list and species list in SSH
    species_list_aer_init='./files/species-list-aer-genoa.dat'

    # output
    path = '{:s}/{:s}/'.format(path,chem)
    os.makedirs(path, exist_ok=True)

    if tag_complete == 2:
        os.system('cp Parameters.py '+path+'/Parameter.sav')
        os.system('cp ReductionStrategy.py '+path+'/ReductionStrategy.sav')

    # reactions head lines
    rout=['SET UNIT GAS MOLCM3','SET TABULATION 11 DEGREES 0. 10. 20. 30. 40. 50. 60. 70. 78. 86. 90.']

    num=0 # reaction index

    for i in range(len(reactions)):

        if not reactions[i].status: continue
        rout.append('%========================='+str(num+1)+'============================')

        rout.append(reactions[i].toSSH())
        num+=1

    # write files
    # filename.reactions rout
    isMark = 0 # check if there is marked reactions (with !!!)
    with open (path+chem+'.reactions','w') as f:
        for i in rout: 
            f.write(i+'\n')
            if '!!!' in i: isMark+=1
        f.write('END\n')

    # filename.species species
    fsps=[]
    fsps.append('File for '+chem+' (name and molar mass)\n')
    fsps.append('# gaseous species # aqueous species\n')
    fsps.append('{:d} 0\n'.format(len(species)))
    fsps.append('---Gas-phase----\n')

    # filename.RO2 RO2
    fRO2=[]
    fRO2.append('# RO2 species in '+chem+' mechanism\n')

    # filename.aer
    faer=[]
    faer_vec = []
    # add init in species-list-aer
    with open (species_list_aer_init,'r') as f: faer0 = f.readlines()
    faer.extend(faer0[:10])
    faer_vec.extend(faer0[:10])

    # record species
    # restore all the record species
    norg = 0 # record the number of orgs
    sname=[]
    for i in species:
        if i.status:
            fsps.append('{:s}    {:6.2f}\n'.format(i.name,i.mass))
            sname.append(i.name)
            if i.organic: norg += 1
            #if i.organic and not i.Radical: # for extract smiles structure in SOAP
            if tag_First and i.organic and not i.Radical:
                faer.append(i.to_species_list_aer(Type = 'smiles', firstRun = True))
            elif i.condensed:
                # with/without strucure
                if i.SMILES=='-':faer.append(i.to_species_list_aer(Type = 'smiles', henry = i.henry))
                else: faer.append(i.to_species_list_aer(Type = 'smiles'))
                if i.SOAPStructure==None:faer_vec.append(i.to_species_list_aer(Type = 'vectors', henry = i.henry))
                else: faer_vec.append(i.to_species_list_aer(Type = 'vectors'))
            if tag_fake and i.Radical: # fake i.RO2 
                sname.append('FA'+i.name)
                fsps.append('FA{:s}    {:6.2f}\n'.format(i.name,i.mass))
            if i.RO2 and i.SMILES != '-' : fRO2.append(i.name+'\n')
    
    # check sps if not with fake species
    if not tag_fake:
        sps = get_species_from_reactions(reactions)
        for i in sps:
            if i not in sname: 
                raise NameError('Datastream species {:s} is not recorded in species list.'.format(i))

    # add water
    faer.append(faer0[-1])
    faer_vec.append(faer0[-1])

    # change number of species if need
    fsps[2] = '{:d} 0\n'.format(len(sname))

    # filename.species species
    with open (path+chem+'.species','w') as f:
        for i in fsps: f.write(i)

    # filename.RO2 RO2
    with open (path+chem+'.RO2','w') as f:
        for i in fRO2: f.write(i)
    if tag_complete == 2: print('SSH output: number of RO2: ', len(fRO2)-1)

    # filename.aer
    if tag_complete or tag_First:
        with open (path+chem+'.aer','w') as f:
            for i in faer: f.write(i)

    with open (path+chem+'.aer.vec','w') as f:
        for i in faer_vec: f.write(i)

    if tag_complete:
        # filename.cst
        with open (path+chem+'.cst','w') as f:
            for i in cst: f.write(i+'\n')

        # get species list
        sps = []
        for i in species: sps.append(i.name)

        # MolProperty
        #out_species(species,path+chem+'.mol',Type='mol')

        # species-list-gas-filename: from species.spack.dat <- generate by SSH

        # output flow chart available generated with graphviz in viz-js.com
        with open (path+chem+'.viz','w') as f:
            # 3 for [products, branch ratio, reactants]
            # shape for products color for reactants # colorscheme = set19
            # colors see: http://www.graphviz.org/doc/info/colors.html
            rccols = {'O3': 'blue', 'OH': 'red', 'NO3': 'green',
                      'NO2': 'yellow3', 'NO': 'orange', 'RO2': 'gray65', 'HO2': 'cyan3',
                      'CO': 'black', 'SO2': 'black', 'H2O': 'orchid3',
                      'J': 'red4', 'CL': 'black'} # darkolivegreen4 # yellowgreen # deeppink2 # antiquewhite3
            if tag_letter: rccols_tmp = {}
            # init
            snote = []
            for i in range(len(sps)): snote.append([[],[],[]]) # record links with others

            for rc in reactions:
                for i in rc.reactants:
                    s = sps.index(i)
                    if species[s].organic:
                        # the other reactants
                        side = None
                        for j in rc.reactants:
                            if j in list(rccols.keys()): side = j
                        if side == None: # try other possibility
                            if 'RO2' in rc.rate.str: side = 'RO2'
                            elif 'H2O' in rc.rate.str: side = 'H2O'
                            elif 'J<' in rc.rate.str: side = 'J'

                        # products and brt
                        for j in range(len(rc.products)):
                            if rc.products[j] not in NokeepSp:
                                snote[s][0].append(rc.products[j])
                                snote[s][1].append(rc.ratiosPD[j])
                                snote[s][2].append(side)
            # output
            f.write('digraph G {\n')
            # write node
            tmp = [[],[],[],[],[]] # svoc, lvoc, elvoc, voc, radical
            rcshape=['box','box','diamond','ellipse','plain']
            msps = []

            for i in range(len(sps)):
                if species[i].status and species[i].organic:
                    if species[i].lump != []: msps.append(species[i].name)
                    if species[i].condensed:
                        if species[i].psat_atm <= 1E-13: tmp[2].append(sps[i]) #elvoc
                        elif species[i].psat_atm <= 1E-9: tmp[1].append(sps[i]) #lvoc
                        else: tmp[0].append(sps[i]) # svoc
                    elif species[i].Radical: tmp[4].append(sps[i]) # radical
                    else: tmp[3].append(sps[i]) # voc

            for i in range(len(tmp)):
                for j in range(len(tmp[i])):
                    if j == 0: f.write('\tnode [shape={:s}];'.format(rcshape[i]))
                    if tmp[i][j] in msps: f.write(' "m{:s}";'.format(tmp[i][j]))
                    else: f.write(' "{:s}";'.format(tmp[i][j]))
                f.write('\n')

            for i in range(len(sps)):
                if sps[i] in msps: sp = 'm'+sps[i]
                else: sp = sps[i]
                if snote[i][0] != []:
                    for j in range(len(snote[i][0])):

                        if snote[i][0][j] in msps: tmp = '\t"{:s}"->"m{:s}"[len=1.00'.format(sp,snote[i][0][j])
                        else: tmp = '\t"{:s}"->"{:s}"[len=1.00'.format(sp,snote[i][0][j])

                        # format
                        if snote[i][1][j] != 1.0: tmp +=', label = "{:.1f}"'.format(float(snote[i][1][j])*100.)
                        if snote[i][2][j] is not None: # products
                            psp = snote[i][2][j]
                            tmp +=',color = {0:s}, fontcolor = {0:s}'.format(rccols[psp])
                            # add letter
                            if tag_letter and snote[i][1][j] == 1.0 and psp not in list(rccols_tmp.keys()): 
                                tmp += ', label = "{:s}"'.format(psp)
                                rccols_tmp[psp] = True
                        f.write(tmp+']\n')
            f.write('}\n')
            if tag_letter: print('rccols_tmp',rccols_tmp)



    # save properties
    out_species(species,path+chem+'.mol',Type='mol')

    if tag_complete == 2: print('SSH output: directory in ',path)

    # check if there is unfinished reactions
    if isMark: print('\n!!!! in output reaction file: ',path+chem+'.reactions',', total ',isMark,' times!!! check!!!')

def read_properties_SSH(species,filename):
    """read property information from files"""
    """
        Objectives:

        Inputs:
            
        Outputs:
    """
    with open (filename) as f: info=f.read().splitlines()
    # species name list
    sps=[]
    for i in species: sps.append(i.name)

    num=0
    for i in info:
        tmp = i.split(' ')
        # found name
        if tmp[1] in sps:
            item=species[sps.index(tmp[1])]
            if item.gamma!=float(tmp[2]):
                item.gamma=float(tmp[2])
                item.henry=float(tmp[3])
                num+=1
    print('Read property information from SSH files. Total: ',len(info),' not same: ',num)

def read_species(sps,filename):
    """
        Objectives:

        Inputs:
            
        Outputs:
    """
    # load files
    with open (filename,'r') as f: info = f.read().splitlines()
    inds=[[],[]] # record species name and index for searching [[names],[ind1:ind2]]
    tmp = []
    for i in range(len(info)):
        # find a species
        if info[i] == '':
            # if not record before
            if tmp == []:
                # first index
                tmp.append(i)
                # get species name
                sname = info[i+1].split('\t')[1].replace('"','')
            elif len(tmp) == 1 :
                # read second index
                tmp.append(i)
                # output in inds
                inds[0].append(sname)
                inds[1].append(tmp)
                # reset
                tmp = []
    # check sname
    if len(set(inds[0])) != len(inds[0]):
        raise AttributeError('DS: species is repeated in file '+filename)

    # transfer sps to species
    species = []

    # init sps if it is empty
    if sps == []: 
        sps = inds[0]
        print('init sps from', filename,len(inds[0]))
        tag_sps = 1
    else: tag_sps = 0

    # read input namelist sps
    for i in sps:
        if i in inds[0]:
            ind = inds[1][inds[0].index(i)]
            species.append(md.Species())
            species[-1].fromTest(info[ind[0]:ind[1]])
            # check status
            if not species[-1].status: 
                print('DS: read species '+i+' is not activated from '+filename)
            # check if is SSHSpeciesInit
            if i in list(SSHSpeciesInit.keys()) and species[-1].mass != SSHSpeciesInit[i]:
                print("!!! DS: read SSHSpeciesInit does not have the same mass",i,species[-1].mass,SSHSpeciesInit[i])
                species[-1].mass = SSHSpeciesInit[i]

        # init using data from Constant.py
        elif i in list(SSHSpeciesInit.keys()):
            species.append(md.Species(i))
            species[-1].mass=SSHSpeciesInit[i]
            print('DS: adding sps existing SSHSpeciesInit', i)
        else:
            if i[0:2] != 'FA':
                raise NameError('DS: species '+i+' not found in the list '+filename)
    # add SSHSpeciesInit
    for i in set(SSHSpeciesInit.keys())-set(sps):
        species.append(md.Species(i))
        species[-1].mass=SSHSpeciesInit[i]

    if tag_sps: return species, sps
    else: return species

def out_species(species,filename,Type='mol'):
    # save in the json file
    with open (filename,'w') as f: f.close()
    if Type == 'mol':
        for i in species: i.toText(filename)

    elif Type == 'SOAP':
        with open ('./files/species_cxx_head','r') as f: infoHead = f.read()
        with open ('./files/species_cxx_end','r') as f: infoEnd = f.read()
        with open ('./files/species_single_cxx','r') as f: info = f.read().splitlines()
        with open (filename,'a+') as f:
            f.write(infoHead)
            for i in species:
                if i.status and i.SOAPStructure is not None: 
                    tmp = i.toSOAP(info)
                    if tmp: f.write(tmp)
                    else:
                        print(i.name, filename)
                        raise OSError('DS: not SOAP output.')
            f.write(infoEnd)

    else: 
        raise TypeError('DS: out_species unknown type: '+ Type)

def update_AIOMFACformat(species,filename, Type = 'SOAP'):

    print('read file in update_AIOMFACformat: ',filename)
    # get names
    sps=[]
    for i in species: sps.append(i.name)

    if Type == 'SOAP':
        with open (filename,'r') as f: info = f.read().splitlines()
        # trim the initialisation in info
        for i in range(len(info)):
            if '====2====' in info[i]:
                # read from the next line
                info = info[i:]
                print('file SOAPorg is trimmed from the beginning to line ',i+1)
                break

        newsp = [[],[],[]] # species name, index, value
        tag = 0
        # read SOAP structure
        for i in info:
            tmp = i.split(' ')
            if len(tmp) > 1:
                # record new species
                if tmp[1] == 'is' : # BCARY is constructed from smiles: CC1=CCCC(=C)C2CC(C)(C)C2CC1
                    if tmp[0] in sps:
                        tag = 1
                        newsp[0].append(tmp[0])
                        newsp[1].append([])#index
                        newsp[2].append([])#values
                    else: tag = 0
                # record index and value
                elif tag and isfloat(tmp[0]):
                    newsp[1][-1].append(int(tmp[0]))#index
                    newsp[2][-1].append(float(tmp[-3]))#value

        # record SOAP structure
        for i in range(len(newsp[0])):
            species[sps.index(newsp[0][i])].SOAPStructure=[newsp[1][i],newsp[2][i]]

    elif Type == 'Table':
        with open (filename,'r') as f: info = f.read().splitlines()[1:]
        title = info[0]

        for i in info:
            items = i.split(',')
            # check name
            if items[0] not in sps: continue
            # check length
            if len(items) != 62:
                raise ValueError('DS: Table length not right. should be 62. '+len(items)) # 56 + name + MWs + CHNO
            newsp = [[],[]]
            # get SOAPStructure
            for j in range(56):
                if items[j+6] not in ['0','']: # exist value
                    newsp[0].append(j)
                    newsp[1].append(float(items[j+6]))
            species[sps.index(items[0])].SOAPStructure = newsp
            species[sps.index(items[0])].dHvap_KJ = dHvap_simpol(newsp)
            print(items[0], newsp, species[sps.index(items[0])].dHvap_KJ)
    else: 
        raise TypeError('DS: update_AIOMFACformat not recognize input type '+Type)

def reaction_merge(reactions,species):
 
    reas = [[],[],[]] # store reactants, types, index
    # sort
    for i in range(len(reactions)):
        rcn = reactions[i]
        if not rcn.status: continue
        tmp = sorted(rcn.reactants)
        if tmp in reas[0]: # find the same reactants
            s = reas[0].index(tmp) # get index in reas
            rate = check_ssh_type(rcn)
            if rate in reas[1][s]:
                reas[2][s][reas[1][s].index(rate)].append(i)                
            else:
                reas[1][s].append(rate)
                reas[2][s].append([i])
        else: # not find, build new
            reas[0].append(tmp)
            reas[1].append([check_ssh_type(rcn)])
            reas[2].append([[i]])

    # merge
    new_rcns = []
    for s,i in enumerate(reas[2]):
        #print(s,reas[0][s],reas[1][s],reas[2][s])
        for p,k in enumerate(i):
            if len(k) == 1: # one reaction, no merge
                new_rcns.append(reactions[k[0]])
            elif len(k) >= 2: # need to merge
                tmp = [reactions[j] for j in k] # merged reactions
                new_rcn = merge_reactions(tmp, reas[1][s][p])
                #print(new_rcn)
                if new_rcn: new_rcns.append(new_rcn)
                else:
                    for j in k: new_rcns.append(reactions[j])
            else:
                raise ValueError('len(reas[2][s]) <= 1.',s,reas[0][s],reas[1][s],reas[2][s])      
    return  new_rcns
    
def reaction_seperate(reactions,species,Type = 'nokdec'):
    """rewrite reactions with seperqte products"""
    new_rcns = []
    sps = []
    # get species list
    for i in species:
        sps.append(i.name)

    # process reactions
    for i in reactions:
        if not i.status: continue
        sp = [] # store species that need to be seperate
        for j in i.products:
            if species[sps.index(j)].organic: #sp.append(j)
                if j not in NokeepSp: 
                    sp.append(j)
        # cases that not seperate
        if len(sp) <= 1 or 'J<' in i.rate.str: # not change reaction
            new_rcns.append(i)
        else:
            # check mass balance # get ratios
            rt = []
            for j in sp:
                rt.append(i.ratiosPD[i.products.index(j)])
            tmp = round(sum(rt),3) # check ratio for kinetic rate

            if tmp < 1.0: # add a destruction if sum(rt) < 1.0
                sp.append('')
                rt.append(1.0 - sum(rt))
                sumrt = 1.
            elif tmp > 1.0:
                sumrt = tmp #sum(rt)
                # need to alter branching ratio
                for j in range(len(rt)):
                    rt[j] /= sumrt
                #print(rt,sumrt,sum(rt))
            else: 
                sumrt = 1.0
                #print('1.0: ',rt,sum(rt))
            
            for j in range(len(rt)): rt[j] = round(rt[j],3)

            if Type == 'kdec':
                tag = 0
                for s in i.products:
                    if species[sps.index(s)].lump != [] : 
                        tag = 1
                        break
                if tag == 0: 
                    new_rcns.append(i)
                    continue

            # build new reactions
            for j,s in enumerate(sp):
                if rt[j] <= 0.0: continue # remove species with very small ratio
                # new reaction for each species
                new_rcns.append(md.Reaction())
                rcn = new_rcns[-1]
                # copy reactants with their ratios
                rcn.reactants = i.reactants
                rcn.ratiosRC = i.ratiosRC

                # trim rate
                rcn.rate.str = i.rate.str+'*{:.3E}'.format(rt[j])# new rate with ratio
                rcn.rate.update() # update rcn.species/rcn.type

                # add products
                if s != '':
                    rcn.products.append(s)
                    rcn.ratiosPD.append(sumrt)
                for k in i.products:
                    if k not in sp: # modify ratio for common species
                        rcn.products.append(k)
                        rcn.ratiosPD.append(i.ratiosPD[i.products.index(k)]/sumrt)
    return  new_rcns

def merge_reactions(rcns, Type):
    """merge reactions if they have the same reactants"""

    if len(rcns) < 2 or not Type: # check number of input reactions
        print('merge_reactions input number of reactions is less than 2! or type is not recognized.',len(rcns),Type)
        return None
    else: # check rate.SSH
        isps,irt = rcns[0].reactants, rcns[0].ratiosRC
        for i in rcns[1:]:
            if not (set(isps) == set(i.reactants) and set(irt) == set(i.ratiosRC)):
                print('reactants or ratiosRC not the same.',set(isps),set(i.reactants),set(irt),set(i.ratiosRC))
                return None
        if set(Type) & set(['MCM1','MCM2','SPEC']): return None
    # new reaction for each species
    new_rcn = md.Reaction()
    new_rcn.reactants = isps
    new_rcn.ratiosRC = irt
    # get lists of rates
    ilts = []
    for i in rcns: ilts.append([x for x in i.rate.SSH.split(' ') if x != ''])

    # get / check index
    inds = []
    for s in Type: 
        tmp = [i.index(s) for i in ilts]
        if len(set(tmp)) == 1: # same index for all reactions
            inds.append(tmp[0])
        else:
            print('index of rate type is not the same.',s,tmp,ilts)
            return None

    # new rate
    tps = ['ARR1','ARR2','ARR3','MCM3']
    lens = {'ARR1':1,'ARR2':2,'ARR3':3,'MCM3':3}
    ns = 0 # merge time
    for s in Type:
        if s in tps:
            if ns == 0: ns += 1
            else: raise ValueError('Find multiple ssh types.',ns,Type,tps)
            ind = inds[Type.index(s)] + 1 # get R1: the constant
            l = lens[s]
            for i in ilts: # check length
                if len(i[ind:]) != l:
                    raise ValueError('len of ssh type does not fit.',s,len(i[ind:]),ilts)
            rts = [float(i[ind]) for i in ilts]
            if s == 'ARR1':
                new_rcn.rate.str = '{:5.2E}'.format(sum(rts))
            elif s == 'ARR2':
                iexp = [x[ind+1] for x in ilts]
                if len(set(iexp)) != 1:
                    print('ARR2 exp value not the same',iexp,ilts)
                    return None
                else: # get negative value
                    new_rcn.rate.str = '{:5.2E}*EXP({:s}/TEMP)'.format(sum(rts),get_negative_str(iexp[0]))
            elif s == 'ARR3':
                itemp = [x[ind+1] for x in ilts]# T**a
                iexp = [x[ind+2] for x in ilts]# exp(b)    
                if len(set(itemp)) == 1 and len(set(iexp)) == 1:
                    new_rcn.rate.str = '{:5.2E}*TEMP**{:s}*EXP({:s}/TEMP)'.format(sum(rts),itemp[0],get_negative_str(iexp[0]))
                else:
                    print('ARR3 c2 and/or c3 value not the same',itemp,iexp,ilts)
                    return None
            elif s == 'MCM3':
                # check photolysis type
                num = [i.rate.str.split('J<')[1].split('>')[0] for i in rcns]
                if len(set(num)) == 1 and isint(num[0]):
                    num = int(num[0])
                    tmp = mcm_photolysis(num)
                    if len(tmp) == 3:
                        if tmp[1] != []: 
                            new_rcn.rate.str = 'J<{:d}>*{:5.2E}'.format(num,sum(rts))
                    else: return None
                return None

    if 'TB' in Type: # treat TB at the end
        ind = inds[Type.index('TB')] + 1
        tmp = [x[ind] for x in ilts]
        if len(set(tmp)) == 1: # one TB value
            new_rcn.rate.str += '*'+tmp[0] # add TB
        else:
            print('multiple TB.',tmp,ilts)
            return None
    # update
    new_rcn.rate.update()

    # get new branching ratios
    for j in range(len(rcns)):
        for i,s in enumerate(rcns[j].products):
            if s == '': continue # no product
            k = rcns[j].ratiosPD[i]*rts[j]/sum(rts) # ratio
            if s not in new_rcn.products:
                new_rcn.products.append(s)
                new_rcn.ratiosPD.append(k)
            else:
                ind = new_rcn.products.index(s)
                new_rcn.ratiosPD[ind] += k

    new_rcn.products = [i for _,i in sorted(zip(new_rcn.ratiosPD, new_rcn.products),reverse=True)]
    new_rcn.ratiosPD = list(sorted(new_rcn.ratiosPD,reverse=True))
    
    return new_rcn

def check_ssh_type(rcn):
    """read the input reaction and return the type in SSH"""
    if '!!!' in rcn.rate.SSH: 
        print('find !!! in rcn.rate.SSH.')
        return None
    else:
        itype = [i for i in rcn.rate.SSH.split(' ') if (not isfloat(i) and  i not in ['KINETIC','']) ]
        if itype == []: raise NameError('No keywords (e.g. SPEC, ARR1, ARR2, ARR3, MCM, TB) is found.',rcn.rate.SSH)
        return itype

def get_negative_str(val):
    if val[0] == '-': return val[1:]
    else: return '-' + val

def dHvap_simpol(SOAPStructure):
    simpol1 = {
                'bo':-9.0677E+02,
                'C':	-2.3229E+02,
                'C=C':	5.9336E+01,
                'OH':	-8.7901E+02,
                'aldehyde':	-5.2267E+02,
                'ketone':	1.9917E+01,
                'COOH':	-1.1963E+03,
                'nitrate':	-7.8253E+02,
                'peroxide':	4.4567E+02,
                'hydroperoxide':	-7.9762E+02,
                'aromatic ring':	-1.3635E+02,
                'ether':	-2.2814E+02,
                'phenol':	-4.2981E+02,
                'nitrophenol':	2.8685E+02
             }
    simpol0 = {
                'C':[0,1,2,3],
                'C=C':[16,17,18,19,20],
                'OH':[26],
                'aldehyde':	[31],
                'ketone':[29,30],
                'COOH':	[37],
                'nitrate':	[39,40,41],
                'peroxide':	[45,46,47,48,49,50,51,52,53],
                'hydroperoxide':[42,43,44],
                'aromatic ring':[21,22],
                'ether':[34,35,36],
                'phenol':[28],
                'nitrophenol':[38]
             }
    ksim = list(simpol0.keys())
    dH = 0.0
    for i,j in enumerate(SOAPStructure[0]):
        tag = 0
        for k in ksim:
            if j in simpol0[k]:
                dH +=  SOAPStructure[1][i] * simpol1[k]
                tag = 1
                break
        if tag == 0: print('!!!!not find strucutre: ', j)
    if dH != 0.0: return -1 * (dH + simpol1['bo']) * (2.303*8.314)/1000
    else: return 0.0

if __name__ == '__main__':
    rc,sp=read_chem_sets('../fromMCM/BCARY/reaction.txt','../fromMCM/BCARY/mcm_subset_mass.txt')
    to_SSH_sets('../toSSH/','BCARY',rc,sp)

