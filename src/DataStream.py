# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  DataStream.py provides functions to input and output of chemical mechanisms.
#
# ================================================================================

import os
from copy import deepcopy

import Module as md
from Parameters import SSHSpeciesInit,NokeepSp, \
                       AeroDict,prefix, \
                       species_list_aer_init, primaryVOCs
from Functions import isfloat,compare,isint, \
                      get_negative_str, \
                      convert_number_format_in_string, \
                      reformat_number_in_string
from ChemRelation import get_species_from_reactions, generation
from KineticMCMtoSSH import KineticSSHtoStr, mcm_photolysis

# if add fake species for radicals

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
    #sps = []
    # store info as md.Reaction() instances 
    reactions, rcn_ratios, rcn_pds = [], [], []

    if reactionType == 'MCM': 
        # process MCM reaction in html format
        with open (reactionfile) as f: fin=f.read().splitlines()

        # test spr
        spr = '→' # '\xe2\x86\x92' for python2

        for info in fin:

            if spr not in info: continue

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
            else:
                if nb > 2: print(a_rt,'DS: pb1 multiple *0. in rate: '+raw_rate)
                # check the part after *0.
                # '1D-16*0.35'
                if isfloat('0.'+a_rt[-1]):
                    ratio='0.'+a_rt[-1]
                # '1D-16*0.35*KDEC'
                elif isfloat('0.'+a_rt[-1].split('*')[0]):
                    ratio='0.'+a_rt[-1].split('*')[0]
                # other cases
                else:
                    print(a_rt)
                    raise ValueError('DS: pb2 *0. not for ratio '+raw_rate)
                # get rate without ratio
                if '*'+ratio in raw_rate and len(raw_rate.split('*'+ratio)) == 2:
                    rate=raw_rate.replace('*'+ratio,'')
                    ratio=float(ratio)
                else:
                    raise ValueError('multiple *ratio','*'+ratio,raw_rate)
                    
            ## record chemical reaction
            tag = -1
            if not tag_seperate and reactions != []:# and ratio != 1:
                for i in range(len(reactions)):
                    if reactions[i].rate.str==rate and compare(reactions[i].reactants,reactants):
                        if ratio >= 1:
                            print('sim rcn has ratio >= 1.',i,ratio,reactants,products,rate)
                        else:
                            if sum(rcn_ratios[i]) + ratio > 1:
                                print('new ratio exceed 1.',i,ratio,sum(rcn_ratios[i])+ratio,reactants,products,rate)
                                print(rcn_ratios[i],rcn_pds[i])
                            tag = i
            if tag != -1:
                rcn=reactions[tag] # merge if rate, reactants are same with the previous one
                rcn_ratios[tag].append(ratio)
                rcn_pds[tag].append(products)
            else:
                # build new reaction
                reactions.append(md.Reaction())
                rcn_ratios.append([ratio])
                rcn_pds.append([products])
                rcn=reactions[-1]
                # add reactants and rate
                for rc in reactants: rcn.record_reactant(rc,1)
                rcn.rate.update(rate)

            # add new products and ratio
            for pd in products: rcn.record_product(pd,ratio)

            # add new species in sps
            #for i in reactants,products:
            #    for j in i:
            #        if j not in sps: sps.append(j)

        # check ratio to make sure it is 1
        i = 0
        for k in range(len(rcn_ratios)):
            if sum(rcn_ratios[k]) == 1: i+=1
            else:
        # in case: NBCO2 + HO2 -> NBCOOH KRO2HO2*0.975
        #for i in range(len(reactions)):
                if round(sum(rcn_ratios[k]),4)==1: 
                    #print('trim ratio before',i,reactions[i].reactants,reactions[i].ratiosPD)
                    for j in range(len(reactions[i].ratiosPD)):
                        reactions[i].ratiosPD[j]=round(reactions[i].ratiosPD[j],3)
                    #print('trim ratio after',reactions[i].ratiosPD)
                    i+=1
                elif len(rcn_ratios[k]) == 1:
               #if 'HO2' in reactions[i].reactants and len(reactions[i].products)== 1 and reactions[i].ratiosPD[0] != 1.0:
                    reactions[i].rate.update(reactions[i].rate.str+'*'+str(reactions[i].ratiosPD[0]))
                    reactions[i].ratiosPD[0] = 1.0
                    i+=1
                else:
                    print('---rcn_ratio != 1', sum(rcn_ratios[k]),reactions[i].reactants,reactions[i].rate.str)
                    print(k,i,rcn_ratios[k],rcn_pds[k])
                    # need to add reactions
                    for j in range(len(rcn_ratios[k])):
                        if j == 0: tmp = deepcopy(reactions[i].rate.str)
                        else:
                            reactions.insert(i+j,md.Reaction())
                            reactions[i+j].reactants=deepcopy(reactions[i].reactants)
                            reactions[i+j].ratiosRC=deepcopy(reactions[i].ratiosRC)
                        #print(tmp+'*'+str(rcn_ratios[k][j]))
                        reactions[i+j].rate.update(tmp+'*'+str(rcn_ratios[k][j]))
                        reactions[i+j].products = rcn_pds[k][j]
                        reactions[i+j].ratiosPD = [1.]*len(rcn_pds[k][j])
                        print('Add new reaction',i,j,reactions[i+j].reactants,reactions[i+j].rate.str,reactions[i+j].products)
                    i+=len(rcn_ratios[k])

       # sps = get_species_from_reactions(reactions)
        if len(reactions) < 1:
            raise ValueError('saperator not found in reaction file. Please check.', spr,reactionfile)

        #sys.exit('here')
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
                if 'J<' in line: rcn.rate.Photolysis = True
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
            #for i in rcn.reactants,rcn.products:
            #    for j in i:
            #        if j not in sps: sps.append(j)

    elif reactionType == 'KPP':
        # process MCM reaction/species in KPP format PRAM
        with open (reactionfile) as f: fin=f.read().splitlines()

        rates = [] # restore rates
        for info in fin:
            if info[0:2] == '//': continue #find comments
            elif '{' in info: # find reactions
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
                #for i in reactants,products:
                #    for j in i:
                #        if j not in sps: sps.append(j)

    else:
        raise TypeError('Unknown type in DS:'+reactionType)

    # species list
    sps = get_species_from_reactions(reactions)
    if speciesType == 'MCM': 
        species = sps_to_species(sps,speciesfile,'MCM')
    elif speciesType == 'SSH': 
        species = read_species(sps,speciesfile)
    elif speciesType == 'KPP':
        # for current scheme
        species0 = sps_to_species(sps,speciesfile,'KPP')
        sps0 = [i.name for i in species0]

        sps1 = list(set(sps)-set(sps0))
        print('species need to be generated from MCM: ',sps1)
        if sps1 != []:
            #for i in ['CARENE']:
            #    if i in sps1: sps1.remove(i)
            #species1 = sps_to_species(sps1,'../../PRAM/mcm_monoterpene_mass.txt','MCM')
            species1 = read_species(sps1,'../toSSH/MPo/MPo.mol')
            sps1 = sps1 + [i.name for i in species1]
            species = species0 + species1
        else: species = species0
    else:
        raise TypeError('DS: speciesType unknown '+speciesType)

    return reactions,species

def sps_to_species(sps,speciesfile,Type = 'MCM'):
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
                print(speciesfile)
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
                # commant for generating PRAM species
                j.update(AeroDict)
                #j.update_advanced(AeroDict)

            else: # record not organic species
                species.append(md.Species(i))

            # read gamma from SSH output file
            #read_properties_SSH(species,'./files/gammainfSSH')

        # read SOAP structure
        inputfile = AeroDict['soapfile']
        #if 'CRI' in speciesfile or 'cri' in speciesfile : update_AIOMFACformat(species,'./files/SOAP_CRI')
        #elif not tag_First: update_AIOMFACformat(species, inputfile)
        if inputfile and os.path.exists(inputfile):
            # read SOAP output structure
            #inputfile = './files/SOAP_'+speciesfile.split('/')[-1]
            print('SOAP output file is used to build aerosol list: ',inputfile)
            update_AIOMFACformat(species, inputfile,'SOAP')

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
                    j.RO2=True # only for MCM species ??
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
    sps=[i.name for i in species]

    if soapfile is not None:
        print('Update new properites from soapfile: ',soapfile)
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
    spr = '\\hline\n'

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
            text = text[0:-3] + ' $\\rightarrow$ '# add ->

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
            textrt = ''

            if 'TB' in rts: textrt += '[{:s}] * '.format(rts[rts.index('TB')+1]) # add TB species
            if 'ARR1' in rts: textrt += rts[rts.index('ARR1')+1]
            elif 'ARR2' in rts: 
                # get the negative value
                if rts[-1][0] == '-': val = rts[-1][1:]
                else: val = '-' + rts[-1]
                textrt += '{:s} * EXP({:s}/T)'.format(rts[rts.index('ARR2')+1],val.replace('.00',''))
            elif 'ARR3' in rts:
                # get the negative value
                if rts[-1][0] == '-': val = rts[-1][1:]
                else: val = '-' + rts[-1]
                textrt += '{:s}*T**{:s}*EXP({:s}/T)'.format(rts[rts.index('ARR3')+1],rts[rts.index('ARR3')+2],val.replace('.00',''))
            elif 'MCM3' in rts: textrt += 'J({:s}, {:s}, {:s})'.format(rts[-3],rts[-2],rts[-1])#photolysis
            else: 
                # check specific types
                if rcn.rate.str in ['KFPAN', 'KBPAN']:
                    textrt += rcn.rate.str
                else:
                    textrt += ('!!!'+rcn.rate.str)

            # change format
            textrt = convert_number_format_in_string(textrt)
            f.write(text + textrt + '\\\\' + spr)
        f.write('\n\n\n')
        
        # species
        #titles = 'Surrogate & Type & Molecular formula & MW$^a$ & P$_{sat}^b$ & $\Delta$H$_{vap}^c$ & K$_p^d$ & H$^e$ & $\gamma^f$ & C* \\\\ \n'
        if soapfile is not None: f.write('Update from SOAP file: '+soapfile+'\n')
        
        # order species for output
        for i in species:
            if not i.condensed: i.psat_atm = 0.
        species1 = sorted(species, key=lambda x: (x.condensed,x.Radical, -x.psat_atm, x.mass, x.name), reverse=False)
            
        for s in species1:
            if s.organic:
                text = ''
                
                # get name
                if s.lump == []: name = s.name
                else: name = 'm' + s.name
                
                # get formula
                formula = ''
                for i in ['C','H','N','O']: 
                    if i in s.functionalGroups:
                        n = int(s.functionalGroups[i])
                        if n != 0:
                            formula += i # add atom
                            if abs(s.functionalGroups[i] - n) <= 0.01:
                                if n != 1: # add int
                                    formula+='$_{{{:d}}}$'.format(n)
                            else: # add float, 2 digit  # remove last zero
                                formula+='$_{{{:.2f}}}$'.format(s.functionalGroups[i])
                                

                # get type
                if s.condensed:
                    if s.psat_atm < 1E-13: stype = 'ELVOC'
                    elif s.psat_atm < 1E-9: stype = 'LVOC'
                    else: stype = 'SVOC'
                    # compute effec c
                    # the calculation was performed for a mean molar mass of 200 g/mol for ideal conditions
                    # C*=1.e6*Mow*gamma*Psat (en torr)/(760*8.202e-5*temperature)
                    # gamma=1 and MOW=200 g/mol and t = 298k
                    effc = 1E6 * 200 * s.gamma * s.psat_atm / (8.202e-5 * 298.)
                    
                    # aerosol properties
                    text0 = '{:s}&{:s}&{:s}&{:6.1f}&{:5.2E}&{:5.2E}&{:5.2E}&{:5.2E}&{:5.2E}'.format(name,stype,
                             formula,s.mass,s.psat_atm,s.dHvap_KJ,s.gamma,s.henry,effc)
                             
                else:
                    if s.Radical: stype = 'Radical'
                    else: stype = 'VOC'
                    
                    # gas-phase species properties
                    text0 = '{:s}&{:s}&{:s}&{:6.1f}&&&&'.format(name,stype,formula,s.mass) 
                    
                text += reformat_number_in_string(text0).replace('10^{0','10^{').replace('10^{-0','10^{-')
                f.write(text + '\\\\' + spr)

def to_SSH_sets(path,chem,reactions,species0,out_mode = '21', tag_conserved = True, tag_fake = False, kept_species = []):
    """return files that need in SSH"""
    """
        Objectives:

        Inputs:
            path: path to save the mechanism
            chem: saved name
            
            reactions: input reaction list 
            species0: input species list
            
            tag_Fake: active if generates fack mechansims from the input scheme
            
            tag_conserved: 1: mass conserve for mainly inorganic species (NokeepSps)
            
            kept_species: a list of species should be kept in the scheme
                  if not kept, print message & stop generation
                   
            out_mode (two digits): output mode for the mechanism
                first digit:
                0: only output .reactions & .mol files
                1: output all essential files for a SSH-aerosol simulation in a folder
                2: output all files related to the mechanism (viz file, latex file)
                second digit:
                
                second digit
                0: strict mode - stop if not kept species or !!! found in the list
                1: non strict mode
                
        Outputs:
            output files related to the input mechanism
    """

    # change species list and reaction list
    sps0 = [i.name for i in species0]
    sps = get_species_from_reactions(reactions)
    species, spRO2 = [],[]

    # check sps
    for i in sps:
        if i in sps0 and species0[sps0.index(i)].status:
            species.append(species0[sps0.index(i)])
        else:
            if i in sps0:  raise NameError('Datastream species {:s} is inactived in species list.'.format(i),path,chem)
            else: raise NameError('Datastream species {:s} is not recorded in species list.'.format(i),path,chem)

    # check kept species
    if kept_species != []:
        kept_sps = [] # check if kept species is reserved in the scheme
        for s in kept_species:
            if s not in sps: kept_sps.append(s)
        if kept_sps != []: 
            print('DS to_SSH_sets: !!! kept species not found: ',kept_sps)
            if out_mode[1] == '0': return len(kept_sps)

    # check rcn conservation: 
    # e.g., A + OH -> B + OH
    if tag_conserved:
        for i in reactions: i.check_conservation()
        
    # add back SSHSpeciesInit
    for i in list(SSHSpeciesInit.keys()):
        if i not in sps:
            species.append(md.Species(i))
            species[-1].mass=SSHSpeciesInit[i]
            sps.append(i)

    # if add fake radicals in reactions
    if tag_fake: # symbol = 'FA'
        for n in range(len(reactions)):
            i=0
            while i < len(reactions[n].products):
                if reactions[n].products[i] in sps:
                    if species[sps.index(reactions[n].products[i])].Radical:
                        # add fake species
                        reactions[n].products.append('FA'+reactions[n].products[i])
                        reactions[n].ratiosPD.append(reactions[n].ratiosPD[i])
                i+=1

    # output folder
    if out_mode[0] == '0': path = '{:s}/'.format(path)
    else: path = '{:s}/{:s}/'.format(path,chem)
    
    # creat folder if needs
    os.makedirs(path, exist_ok=True)

    ### save properties in .mol file
    out_species(species,path+chem+'.mol',Type='mol')
    
    ### prepare for reaction list: .reactions
    # reactions head lines
    rout=['SET UNIT GAS MOLCM3','SET TABULATION 11 DEGREES 0. 10. 20. 30. 40. 50. 60. 70. 78. 86. 90.']
    num=0 # reaction index

    for i in range(len(reactions)):

        if not reactions[i].status: continue
        rout.append('%========================='+str(num+1)+'============================')
        rout.append(reactions[i].toSSH())
        # check ro2 species
        if "TB RO2" in rout[-1]:
            for s in reactions[i].reactants:
                # find RO2 species with RO2 reactions
                if s in sps0 and species0[sps0.index(s)].organic and s not in spRO2:
                    spRO2.append(s)
        num+=1

    # filename.reactions out
    isMark = 0 # check if there is marked reactions (with !!!)
    with open (path+chem+'.reactions','w') as f:
        for i in rout: 
            f.write(i+'\n')
            if '!!!' in i: isMark+=1
        f.write('END\n')
        
    if isMark > 0 and out_mode[1] == '0':
        print('DS to_SSH_sets: !!! in output reaction file: ',path+chem+'.reactions',', total ',isMark,' times!!! check!!!')
        return isMark

    # return 0 if only out files for GENOA reduction
    if out_mode[0] == '0': return 0
    
    ### prepare to save other files for SSH-aerosol simulations
    # filename.species species
    fsps=[]
    fsps.append('File for '+chem+' (name and molar mass)\n')
    fsps.append('# gaseous species # aqueous species\n')
    fsps.append('{:d} 0\n'.format(len(species)))
    fsps.append('---Gas-phase----\n')

    # filename.RO2 RO2
    fRO2=[]
    fRO2.append('# RO2 species in {:s} mechanism. Find {:d} RO2 species with RO2 reactions.\n'.format(chem,len(spRO2)))

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
            #if tag_First and i.organic and not i.Radical:
            #    faer.append(i.to_species_list_aer(Type = 'smiles', firstRun = True))
            #el
            if i.condensed:
                faer.append(i.to_species_list_aer(Type = 'smiles'))
                faer_vec.append(i.to_species_list_aer(Type = 'vectors'))
            if tag_fake and i.Radical: # fake i.RO2 
                sname.append('FA'+i.name)
                fsps.append('FA{:s}    {:6.2f}\n'.format(i.name,i.mass))
            #if i.RO2 and i.SMILES != '-' : 
            if i.name in spRO2 and i.SMILES != '-' : # not input PRAM species
                #fRO2.append('{:s}\t{:d}\n'.format(i.name,i.groupID))
                fRO2.append('{:s}\n'.format(i.name))
                #if not i.RO2: print('to_SSH_set: species not RO2 are with RO2-RO2 reactions: ',i.name)
            #elif i.RO2 and i.SMILES != '-' :
            #    print('to_SSH_set: RO2 species without RO2 reaction: ',i.name)
                #species[species.index(i)].RO2 = False
    # add water
    faer.append(faer0[-1])
    faer_vec.append(faer0[-1])

    # change number of species if need
    fsps[2] = '{:d} 0\n'.format(len(sname))

    ### output files    
    # filename.species species
    with open (path+chem+'.species','w') as f:
        for i in fsps: f.write(i)
    # filename.RO2 RO2
    with open (path+chem+'.RO2','w') as f:
        for i in fRO2: f.write(i)

    # filename.aer.vec aerosol species list with vector
    with open (path+chem+'.aer.vec','w') as f:
        for i in faer_vec: f.write(i)

    # other files
    if out_mode[0] == '2':

        # filename.aer
        with open (path+chem+'.aer','w') as f:
            for i in faer: f.write(i)

        # filename.cst
        #with open (path+chem+'.cst','w') as f:
        #    for i in cst: f.write(i+'\n')

        # filename.viz
        # output flow chart available generated with graphviz in viz-js.com
        viz_from_chem(reactions,species,path+chem+'.viz')

        # filename.latex
        # output species and reaction lists as a latex table
        #to_latex_sets(path,chem,reactions,species,soapfile = None)
        
    # printout check if there is unfinished reactions
    tmp = 0
    if kept_species != [] and kept_sps != []: 
        print('DS to_SSH_sets: !!! kept species not found: ',kept_sps)
        tmp += len(kept_sps)
    if isMark:
        print('DS: !!! in output reaction file: ',path+chem+'.reactions',', total ',isMark,' times!!! check!!!')
        tmp += isMark
        
    return tmp

def read_properties_SSH(species,filename):
    """read property information from files"""
    """
        Objectives:

        Inputs:
            
        Outputs:
    """
    with open (filename) as f: info=f.read().splitlines()
    # species name list
    sps=[i.name for i in species]
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
                sname = info[i+1].split('"')[-2]
                #sname = info[i+1].split('\t')[1].replace('"','')
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
            species[-1].fromText(info[ind[0]:ind[1]])
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

    if filename and os.path.exists(filename): 
        print('read file in update_AIOMFACformat: ',filename)
    else:
        print('generate functional groups using UManSysProp in update_AIOMFACformat.')

    # get names
    sps=[i.name for i in species]

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
        tag = 0
        for i in range(len(newsp[0])):
            print('Update species structure from SOAP: ', species[sps.index(newsp[0][i])].name)
            species[sps.index(newsp[0][i])].SOAPStructure=[newsp[1][i],newsp[2][i]]
            tag += 1
        print('In total read Structure from SOAP: ', tag, ' species.')

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
            #species[sps.index(items[0])].dHvap_KJ = dHvap_simpol(newsp)
            print(items[0], newsp, species[sps.index(items[0])].dHvap_KJ)
    else: 
        raise TypeError('DS: update_AIOMFACformat not recognize input type '+Type)

def viz_from_chem(reactions, species, filename = 'chem.viz'):

    # 3 for [products, branch ratio, reactants]
    # shape for products color for reactants # colorscheme = set19
    # colors see: http://www.graphviz.org/doc/info/colors.html
    rccols = {'O3': 'blue', 'OH': 'red', 'NO3': 'green',
              'NO2': 'yellow3', 'NO': 'orange', 'RO2': 'gray65', 'HO2': 'cyan3',
              'CO': 'black', 'SO2': 'black', 'H2O': 'orchid3',
              'J': 'red4', 'CL': 'black'} # darkolivegreen4 # yellowgreen # deeppink2 # antiquewhite3
    rcarrowheads = {'O3': 'vee', 'OH': 'odiamond', 'NO3': 'onormal',
              'NO2': 'box', 'NO': 'dot', 'RO2': 'diamond', 'HO2': 'odot',
              'CO': 'normal', 'SO2': 'normal', 'H2O': 'tee',
              'J': 'obox', 'CL': 'tee'} 
    tag_letter = False # output viz with letter if possible
    
    sps = [s.name for s in species]
    
    with open (filename,'w') as f:

        if tag_letter: rccols_tmp = {}
        # init
        snote = []
        for i in range(len(sps)): snote.append([[],[],[]]) # record links with others

        for rc in reactions:
            if not rc.status: continue
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
        # layout
        f.write('    layout=dot; graph [ranksep=0.8];\n')
        if 'APINENE' in primaryVOCs:
            f.write('APINENE [fontsize="50",style="filled",penwidth=0,shape="ellipse",fillcolor="#ff880022"]\n')
        if 'BPINENE' in primaryVOCs:
            f.write('BPINENE [fontsize="50",style="filled",penwidth=0,shape="ellipse",fillcolor="#0044ff22"]\n')
        if 'LIMONENE' in primaryVOCs:
            f.write('LIMONENE [fontsize="50",style="filled",penwidth=0,shape="ellipse",fillcolor="#00ff0022"]\n')
        f.write('\n')

        iline_layout = 'penwidth=3,len=1,arrowsize=3'
        isps_layout = 'fontsize="50"'
        
        # write node
        tmp = [[],[],[],[],[]] # svoc, lvoc, elvoc, voc, radical
        rcshape=['box','hexagon','octagon','ellipse','plain']
        msps = []

        for i in range(len(sps)):
            if species[i].status and species[i].organic:
                if species[i].name in ['APINENE','BPINENE','LIMONENE']: continue
                if species[i].lump != []: msps.append(species[i].name)
                if species[i].condensed:
                    if species[i].psat_atm <= 1E-13: tmp[2].append(sps[i]) #elvoc
                    elif species[i].psat_atm <= 1E-9: tmp[1].append(sps[i]) #lvoc
                    else: tmp[0].append(sps[i]) # svoc
                elif species[i].Radical: tmp[4].append(sps[i]) # radical
                else: tmp[3].append(sps[i]) # voc

        for i in range(len(tmp)):
            for j in range(len(tmp[i])):
                if j == 0: f.write('\tnode [{:s},shape={:s}];'.format(isps_layout,rcshape[i]))
                if tmp[i][j] in msps: f.write(' "m{:s}";'.format(tmp[i][j]))
                else: f.write(' "{:s}";'.format(tmp[i][j]))
            f.write('\n')

        for i in range(len(sps)):
            if sps[i] in msps: sp = 'm'+sps[i]
            else: sp = sps[i]
            if snote[i][0] != []:
                for j in range(len(snote[i][0])):

                    if snote[i][0][j] in msps: tmp = '\t"{:s}"->"m{:s}"[{:s}'.format(sp,snote[i][0][j],iline_layout)
                    else: tmp = '\t"{:s}"->"{:s}"[{:s}'.format(sp,snote[i][0][j],iline_layout)

                    # format
                    #if snote[i][1][j] != 1.0: tmp +=', label = "{:.1f}"'.format(float(snote[i][1][j])*100.)
                    if snote[i][2][j] is not None: # products
                        psp = snote[i][2][j]
                        tmp +=',color ={0:s},fontcolor={0:s},arrowhead={1:s}'.format(rccols[psp],rcarrowheads[psp])
                        # add letter
                        if tag_letter and snote[i][1][j] == 1.0 and psp not in list(rccols_tmp.keys()): 
                            tmp += ', label = "{:s}"'.format(psp)
                            rccols_tmp[psp] = True
                    f.write(tmp+']\n')
        # add rank
        tmp = '{rank = min'
        for j in primaryVOCs: tmp+=(';'+j)
        f.write(tmp+'}\n') 
          
        # end
        f.write('}\n')
        if tag_letter: print('rccols_tmp',rccols_tmp)

if __name__ == '__main__':
    rc,sp=read_chem_sets('../fromMCM/BCARY/reaction.txt','../fromMCM/BCARY/mcm_subset_mass.txt')
    to_SSH_sets('../toSSH/','BCARY',rc,sp)

