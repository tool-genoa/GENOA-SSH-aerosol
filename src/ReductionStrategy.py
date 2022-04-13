# -*- coding: utf-8 -*-
#================================================================================
#
#     GENOA v1.0: the GENerator of reduced Organic Aerosol mechanism
#
#     Copyright (C) 2022 CEREA (ENPC) - INERIS.
#     GENOA is distributed under GPL v3.
#
#================================================================================

import math
import numpy as np
from copy import deepcopy 

from Parameters import primaryVOCs,Roptions,NokeepSp
from SSHResultProcess import get_conc, ref_conc
from KineticMCMtoSSH import MCMtoSSH_kinetic_rate
from Functions import duplicates,list_all_fit, \
                      isfloat,compare, \
                      isFitScale,array_to_ratio
from ChemRelation import get_species_from_reactions,lifetime, \
                         productions,destructions, \
                         branch_ratio,tree,generation

lump_props = ['New Species','name','lumpRatio','lumps','formula','SMILES','mass',\
              'RO2','psat_atm','dHvap_KJ','rConcs','lifetime']

read_props = ['name','formula','SMILES','mass','RO2','psat_atm','dHvap_KJ'] # read from species

# read options
if 'OrderType'in Roptions.keys() and Roptions['OrderType'] in ['conc', 'reaction']:
    OrderType = Roptions['OrderType']
    #print('RS: read OrderType from Roptions: ', OrderType)
else: OrderType = 'reaction'

# type of candidates lumped species for A
if 'PeerType'in Roptions.keys() and Roptions['PeerType'] in ['reaction', 'all', 'psat']: 
    PeerType = Roptions['PeerType']
    #print('RS: read PeerType from Roptions: ', PeerType)
else:  PeerType = 'all'

# if A and candidates can produce one or the other
if 'CheckLink'in Roptions.keys() and isfloat(Roptions['CheckLink']): 
    CheckLink = Roptions['CheckLink']
    #print('RS: read CheckLink from Roptions: ', CheckLink)
else: CheckLink = False # default

# treatment for RO2 species
if 'RO2Treat'in Roptions.keys(): RO2Treat = Roptions['RO2Treat']
else: RO2Treat = None 
#if RO2Treat: print('RS: read RO2Treat from Roptions: ', RO2Treat)
# (fixed is for remove by Psat_S/NVOC; remove RO2 after lumping)

# if check lifetime, No.generation, saturation vapor pressure, stracture
if 'CheckTau'in Roptions.keys(): CheckTau = Roptions['CheckTau'] # False, value of tolerance scale
else: CheckTau = False
#if CheckTau: print('RS: read CheckTau from Roptions: ', CheckTau)

if 'CheckGen'in Roptions.keys(): CheckGen = Roptions['CheckGen']
else: CheckGen = False
#if CheckGen: print('RS: read CheckGen from Roptions: ', CheckGen)

if 'CheckPsat'in Roptions.keys(): CheckPsat = Roptions['CheckPsat'] # False, value of tolerance scale
else: CheckPsat = False
#if CheckPsat: print('RS: read CheckPsat from Roptions: ', CheckPsat)

if 'CheckStructure'in Roptions.keys() and isinstance(Roptions['CheckStructure'],dict):
    CheckStructure = Roptions['CheckStructure']
    CheckStructureKeys = list(Roptions['CheckStructure'].keys())
    #print('RS: read CheckStructure from Roptions: ', CheckStructure)
else: CheckStructure, CheckStructureKeys = [], []

# freeze compounds
if 'FreezeSpecies' in Roptions.keys() and isinstance(Roptions['FreezeSpecies'], list): 
    NoLumps = Roptions['FreezeSpecies']
    #print('RS: read FreezeSpecies from Roptions. Input frozen species: ',NoLumps)
else: NoLumps = []

def reduction_chemistry_scheme(reactions,species,Type,criteria, concs_ave, concs_min, refconc_paths, RefConcRead):
    """the main function used for reduction:
        reduce reactions and species depending on the input reduction type and criteria"""
    # reduction keywords
    # Bratio, lump
    # lump species
    Type_lump=['kdec','tau','jump','replace']
    # remove trees
    Type_tree=['number','Psat_NVOC','Psat_SVOC','gen','conc']

    if Type in Type_tree :
        # remove species tree and all related reactions
        return reduce_reaction_byTree(reactions,species,Type,criteria)

    elif Type in Type_lump:
        # remove species and set their precursors form their products directly
        lsps,lpros = reduce_reaction_bySpecies(reactions,species,Type,criteria)
        return lsps

    elif Type == 'Bratio':
        # removing by branch ratio
        # RRR: need to rewrite the criteria
        Type = int(criteria.split('_')[1])
        Type = 'simple_'+str(Type)
        criteria = float(criteria.split('_')[0])

        lsps, lrcn = reduce_reaction_byRemove(reactions,species,Type,criteria,concs_min)
        return lsps

    elif Type == 'lump':
        # return lump species and their properties
        lsps,lpros,lchn = reduce_reaction_byLump(reactions,species,criteria,refconc_paths,concs_ave, RefConcRead = RefConcRead)
        return lsps

    else:
        raise TypeError('RS reduction_chemistry_scheme: not recognize type of reduction '+Type)

def reduce_reaction_byLump(reactions,species,criteria, 
                           refconc_paths, concs, 
                           frozen = [], RefConcRead = None,
                           target = [], lumptype = ['tau']):
    """prepare species and reactions depending on the given lump species list"""

    # control the maximumlump times
    if len(criteria.split('_')) == 1: 
        tag_ctl= int(criteria)#0 # no limitation
    else:
        tag_ctl = int(criteria.split('_')[1])
        criteria = criteria.split('_')[0]

    # get species name
    sps=[]
    for i in species:sps.append(i.name)

    # get names in reactions
    rcSps=get_species_from_reactions(reactions)

    # only keep organic species in rcSps and rcLump
    n=0 # index in rcSps
    rcLump=[] #[name, lumped species]
    while n<len(rcSps):
        i = species[sps.index(rcSps[n])]
        # check if i is actived
        if not i.status:
            raise AttributeError('RS: non active species is in the reaction list, species '+i.name)
        if i.organic: 
            n+=1
            # record lump in rcLump
            rcLump.append([i.name,i.lump+i.jump])
        else:
            rcSps.remove(i.name)

    # get the reference concentration: use 'fake' conc for radicals
    rcSps,rcConcs=ref_conc(refconc_paths,species=rcLump,
                           Type='sum',radical='FA',savname=RefConcRead)

    # get properties in the current reaction list
    # number of generation
    gen_sps,gens = [], [] #generation(reactions,'None')

    # lifetime
    tau_sps,taus=lifetime(reactions,concs)

    # production # return species name, [products]
    peers = productions(reactions,Type='pro')

    if CheckLink:
        plist = productions(reactions,Type='rea')
        for i in range(len(plist[0])):
            n = 0
            while n < len(plist[1][i]): # extend
                s = plist[1][i][n]
                if s in plist[0] and species[sps.index(s)].organic:
                    for k in plist[1][plist[0].index(s)]:
                        if k not in plist[1][i] and species[sps.index(k)].organic: plist[1][i].append(k)
                n += 1

        dlist = destructions(reactions,Type='pro')
        for i in range(len(dlist[0])):
            n = 0
            while n < len(dlist[1][i]): # extend
                s = dlist[1][i][n]
                if s in dlist[0] and species[sps.index(s)].organic:
                    for k in dlist[1][dlist[0].index(s)]:
                        if k not in dlist[1][i] and species[sps.index(k)].organic: dlist[1][i].append(k)
                n += 1


    # get lumping order for species with conc
    SpOrderInd = []
    if OrderType == 'reaction':
        # the reversed order of products in the reaction list
        nrs = len(reactions)
        for i in range(nrs):
            # reversed index
            ind = nrs - i - 1
            if not reactions[ind].status: continue
            for s in reactions[ind].products:
                if s in rcSps:
                    j = rcSps.index(s)
                    if j not in SpOrderInd: SpOrderInd.append(j)
    elif OrderType == 'conc':
        # ordered by the refered concentration and generation number
        tmp = sorted(rcConcs)
        tmp1=[]
        # in case species with the same conc
        for i in range(len(tmp)):
            if tmp[i] not in tmp1:
                tmp1.append(tmp[i])
                # get all species with the same conc
                jj = np.where(rcConcs==tmp[i])[0]
                # use generation number to reverse the order
                k = []
                for j in jj: k.append(gens[gen_sps.index(rcSps[j])])
                k=[x for _,x in sorted(zip(k,jj))]
                for j in range(len(k)): SpOrderInd = np.append(SpOrderInd,k[len(k)-j-1])
    # check
    if SpOrderInd == []: print('RS: SpOrderInd is empty')

    # lumping preparation
    lumps=[] # lump info
    lumpItems=[] # store lumped species
    ctl = 0 # lumping times

    # try lumping in the ordered list
    for n in SpOrderInd:

        # species name
        sname=rcSps[int(n)]

        # checking
        if target != [] and sname not in target: continue
        if sname in lumpItems or sname in primaryVOCs: continue
        if sname in NoLumps: continue

        if CheckLink:
            if sname in plist[0]: 
                plist_now = plist[1][plist[0].index(sname)]
            else:
                plist_now = []
            if sname in dlist[0]:
                dlist_now = dlist[1][dlist[0].index(sname)]
            else: 
                dlist_now = []

        # not in species list
        if sname not in sps:
            raise NameError('RS: not in species list '+int(n)+' '+sname)

        # species status
        msp = species[sps.index(sname)]
        if not msp.status:
            print('RS lumping: in species list but not actived',int(n),sname)
            raise AttributeError

        # do not lump ROO and RO fast species:
        if msp.Radical and msp.name[-1] == 'O': continue

        # select possible lumping target
        SpList = []
        # consider all organic species
        if PeerType == 'all' :
            for ips in rcSps:
                if species[sps.index(ips)].status:
                    if CheckLink and (ips in plist_now or ips in dlist_now): continue
                    elif ips not in SpList: SpList.append(ips)
        # consider species in the same reactions
        elif PeerType == 'reaction':
            # only check peer species
            if sname in peers[0]: SpList = peers[1][peers[0].index(sname)]
            # check the peers for lumped species
            if msp.lump:
                for i in msp.lump:
                    if i not in peers[0]: continue
                    for j in peers[1][peers[0].index(i)]:
                        if CheckLink and (j in plist_now or j in dlist_now): continue
                        if j not in SpList: SpList.append(j)
        elif PeerType == 'psat':
            tmp = species[sps.index(sname)].psat_atm # store psat of A
            if tmp != 0.0: tmp = math.log(tmp) # use the log value
            tmp1 = {} # store the psat for SpList
            for ips in rcSps:
                if species[sps.index(ips)].status:
                    if CheckLink and (ips in plist_now or ips in dlist_now): continue
                    if ips not in SpList:
                        # select items with similar past in the front
                        dnow = species[sps.index(ips)].psat_atm
                        if tmp != 0.0:
                            if dnow != 0.0: dnow = abs(tmp-math.log(dnow))
                            else: dnow = 1E99 # try at last

                        for i in SpList:
                            if dnow < tmp1[i]:
                                SpList.insert(SpList.index(i),ips)
                                tmp1[ips] = dnow
                                break

                        if ips not in SpList: 
                            SpList.append(ips)
                            tmp1[ips] = dnow
        elif PeerType == 'svoc': # only check psat for svoc species
            if species[sps.index(sname)].condensed:
                tmp = math.log(species[sps.index(sname)].psat_atm) # store psat of A # use the log value
                tmp1 = {} # store the psat for SpList
                for ips in rcSps:
                    if species[sps.index(ips)].status:
                        if CheckLink and (ips in plist_now or ips in dlist_now): continue
                        if ips not in SpList:
                            # select items with similar past in the front
                            dnow = species[sps.index(ips)].psat_atm
                            if tmp != 0.0:
                                if dnow != 0.0: dnow = abs(tmp-math.log(dnow))
                                else: 
                                    dnow = 1E99 # try at last
                                    print('RS: find condensed species with no psat.',ips) 
                                    raise AttributeError

                            for i in SpList:
                                if dnow < tmp1[i]:
                                    SpList.insert(SpList.index(i),ips)
                                    tmp1[ips] = dnow
                                    break

                            if ips not in SpList: 
                                SpList.append(ips)
                                tmp1[ips] = dnow
            else:
                for ips in rcSps:
                    if species[sps.index(ips)].status:
                        if CheckLink and (ips in plist_now or ips in dlist_now): continue
                        elif ips not in SpList: SpList.append(ips)

            SpList1 = []
            for ips in rcSps:
                if species[sps.index(ips)].status:
                    if CheckLink and (ips in plist_now or ips in dlist_now): continue
                    elif ips not in SpList1: SpList1.append(ips)
        # check
        if frozen is not []:
            # check and remove frozen pairs
            for i in frozen:
                if sname in i: # if sname in the list
                    if len(i) == 1:
                        SpList = []
                    for j in i:
                        if j in SpList: SpList.remove(j)
        for i in NoLumps:
            if i in SpList:
                SpList.remove(i)
        if SpList == []: continue

        # search lumping
        tmplump = SpeciesFitLump(msp,SpList,[sps,species],[rcSps,rcConcs],
                                 [gen_sps,gens],[tau_sps,taus],criteria,lumptype)

        # record lumping
        if tmplump:
            lumps.append(tmplump)
            ctl += 1 # lump time
            for i in tmplump[0]: 
                if i not in lumpItems:
                    lumpItems.append(i) # add name
                else:
                    raise ValueError('RS: repeat lumped species in lumpItems.',lumpItems,i,tmplump)

        # lump times control
        # only check tag_ctl times effective species
        if tag_ctl and ctl == tag_ctl: break

    # check if no lump
    if ctl == 0: return [],[],[]

    # update lumping results
    changes,newsps,lumpinfo = [],[],[]

    # newsps
    for l in lumps:
        if l[0] not in newsps: newsps.append(l[0])
        else: raise ValueError('RS: find lumpsps in newsps!',newsps,l[0],lumps)

    # get original species/reactions
    lp_sps, lp_rcn = {}, {} # 'id': need to be changed
    # species
    for s in lumpItems: 
        i = sps.index(s) 
        lp_sps[i] = species[i]
    # reactions
    for i,rcn in enumerate(reactions):
        tmp = [i for i in rcn.reactants + rcn.products if i in lumpItems]
        if tmp != []: # find species
            lp_rcn[i] = reactions[i]
    # prepare copies for multiple lumptype
    if len(lumptype) > 1:
        changes = [[lp_sps, lp_rcn]]
        for n in range(len(lumptype)-1):
            changes.append([deepcopy(lp_sps),deepcopy(lp_rcn)])

    for n in range(len(lumptype)):
        if n > 0: lp_sps, lp_rcn = changes[n]
        # change species properties
        for l in lumps:
            # order species by ratio
            lumpsps = sorted(l[0], key=dict(zip(l[0],l[1][n])).get,reverse=True)
            lumpratio = sorted(l[1][n], reverse=True)

            # name of the new surrogate with maximum ratio
            newsp = lumpsps[0]
            # name of previously lumped species
            lumpTot = [lumpsps[1]]
            for s in lumpsps:
                lumpTot += species[sps.index(s)].lump

            # update species
            # properties in order: average MWs, Psat, dHvap_KJ, functional groups, henry, gamma
            tmp = [0.0,1,0.0,[0.0]*56, {}, 1., 1.]
            for s in lumpsps:
                ind = sps.index(s) # species index
                i = lp_sps[ind] # species
                rt = lumpratio[lumpsps.index(s)] # ratio

                if s != newsp: i.status = False # status
                i.lump = lumpTot # lumped species
                # calculate new values
                tmp[0] += i.mass*rt # MWs 
                tmp[1] *= pow(i.psat_atm,rt) # psat_atm
                tmp[2] += i.dHvap_KJ * rt # dHvap_KJ
                if i.SOAPStructure is not None:
                    for k in range(56):
                        if k in i.SOAPStructure[0] :
                            tmp[3][k]+=(i.SOAPStructure[1][i.SOAPStructure[0].index(k)])*rt
                for k in i.functionalGroups.keys(): 
                    if k not in tmp[4].keys(): tmp[4][k] = i.functionalGroups[k] * rt
                    else: tmp[4][k] += i.functionalGroups[k] * rt
                tmp[5] *= pow(i.henry,rt) #henry

            # assign new values
            i = lp_sps[sps.index(newsp)]
            i.mass = tmp[0]
            i.psat_atm = tmp[1]
            i.psat_torr = i.psat_atm*760.
            i.dHvap_KJ = tmp[2]
            i.SOAPStructure = [[],[]]
            i.henry = tmp[5]
            i.gamma = tmp[6]
            for k in range(56):
                if tmp[3][k] != 0.0: 
                    i.SOAPStructure[0].append(k)
                    i.SOAPStructure[1].append(tmp[3][k])
            i.functionalGroups = tmp[4]
            if 'C' in i.functionalGroups.keys() and i.functionalGroups['C'] != 0.:
                tmp = i.functionalGroups['C']
                i.DU = 2 + 2 * i.functionalGroups['C'] # DU: 1/2 (2 + 2C + N - H - X)
                # ratio
                i.ratios['OM/OC'] = round(i.mass / (tmp*12.),3)
                if 'H' in i.functionalGroups.keys():
                    i.ratios['H/C'] = round(i.functionalGroups['H'] / (tmp),3)
                    i.DU -= i.functionalGroups['H']
                if 'O' in i.functionalGroups.keys():
                    i.ratios['O/C'] = round(i.functionalGroups['O'] / (tmp),3)
                if 'N' in i.functionalGroups.keys():
                    i.ratios['N/C'] = round(i.functionalGroups['N'] / (tmp),3)
                    i.DU += i.functionalGroups['N']
                i.DU /= 2.

            # lump reactions
            for rcn in lp_rcn.values(): LumpReaction(rcn,newsp,lumpsps,lumpratio)

            # prepare output info for species properties
            lump_prop = deepcopy(lump_props)
            #lump_props = ['lumpRatio','lumps','formula','SMILES','mass',\
            #               'RO2','psat_atm','dHvap_KJ','rConcs','lifetime']
            lump_prop[lump_props.index('New Species')] += '\tm{:s}'.format(newsp)
            # record previous lumped species
            lump_prop[lump_props.index('lumps')] += '\t{:}'.format(lumpTot)
            # read properties
            for i in lumpsps:
                ind = sps.index(i) # species index
                j = lp_sps[ind] # species in lp
                # ratio
                s = lump_props.index('lumpRatio')
                if i == lumpsps[0]: lump_prop[s] +=  '\t{:}'.format(lumptype[n])
                lump_prop[s] += '\t{:}'.format(lumpratio[lumpsps.index(i)])
                # record properties
                for k in read_props:
                    s = lump_props.index(k)
                    lump_prop[s] += '\t{:}'.format(eval('j.'+k))
                # 'rConcs' [rcSps,rcConcs]
                s = lump_props.index('rConcs')
                if i in rcSps: lump_prop[s] += '\t{:}'.format(rcConcs[rcSps.index(i)])
                else: lump_prop[s] += '\t-'
                # 'lifetime' [tau_sps,taus]
                s = lump_props.index('lifetime')
                if i in tau_sps: lump_prop[s] += '\t{:}'.format(taus[tau_sps.index(i)])
                else: lump_prop[s] += '\t-'

            lumpinfo.append(lump_prop)

    return newsps,lumpinfo,changes

def reduce_reaction_bySpecies(reactions,species,Type,
                               criteria,frozen = [],target = []):
    """
        Inputs:
            remove species according certain ceriteria
            reactions: reaction list
            species: species list
            type: what type of reduction
            criteria: how many reudctions are wannted
            forzen: if certain combination is not wantted to be removed
            target: if only want to conduct relavant reduction
        Outputs:
            
     """

    # get species that need to be lumped
    lumpsps = [] # species list
    lumppds = [] #  products, prodauct ratios

    sps = [i.name for i in species]

    # remove certain reactions
    if Type == 'kdec': 
        # mute reaction with the type of kdec and replace the kdec reactants by all the kdec products
        for rcn in reactions:
            
            if not rcn.status : continue

            # check kdec species
            #if rcn.rate.str == 'KDEC':
            if 'KDEC' in rcn.rate.str:
                if rcn.rate.str == 'KDEC': ratio = 1.0
                else: 
                    tmp = rcn.rate.str.replace('KDEC','').replace('*','')
                    if isfloat(tmp): ratio = float(tmp)
                    else:
                        raise ValueError('RS: kdec can not get ratio from: '+rcn.rate.str)

                if len(rcn.reactants) != 1 :
                    print(rcn.reactants)
                    raise NameError('RS kdec: with multiple reactants')

                s,tag = rcn.reactants[0],1
                # check if in other reactions as reactant:
                for rcn1 in reactions:
                    if rcn1.status and s in rcn1.reactants and rcn1 != rcn:
                        print('Can not remove kdec species cuz it is also in other reactions.',rcn1.rate.str,rcn1.reactants,s)
                        tag = 0
                        break
                if tag:

                    if s in lumpsps:
                        ind = lumpsps.index(s)
                        for i in rcn.products: 
                            if i in lumppds[ind][0]: # update ratio
                                lumppds[ind][1][lumppds[ind][0].index(i)] += rcn.ratiosPD[rcn.products.index(i)]*ratio
                            else: 
                                lumppds[ind][0].append(i) # add pd
                                lumppds[ind][1].append(ratio) # ratio, KDEC base
                    else:
                        pdrt = []
                        if ratio == 1.0: pdrt = deepcopy(rcn.ratiosPD)
                        else: 
                            for i in rcn.ratiosPD: pdrt.append(i*ratio)

                        lumppds.append([deepcopy(rcn.products),pdrt])
                        lumpsps.append(s)

                    # mute this reaction in output
                    rcn.status = False
                    #break # check only one kdec species

        for rcn in reactions:
            if rcn.status and set(rcn.reactants)&set(lumpsps):
                print(rcn.rate.str,set(rcn.reactants)&set(lumpsps))
                raise AttributeError('RS: kdec is mixed with other reactions')

    elif Type == 'tau':
        # remove species that react fast (kinetic rate > certain value)

        #TEMP = 244.97 # the lowest temperature in profile_conc_2015.nc
        tag_multiple = False

        sps_num = [[],[],[]] # name and number and index in reaction
        for i in range(len(reactions)):
            rcn = reactions[i]
            if not rcn.status : continue
            for s in rcn.reactants:
                if s in sps_num[0]: sps_num[1][sps_num[0].index(s)].append(i)
                else:
                    sps_num[0].append(s)
                    sps_num[1].append([i])

        for rcn in reactions:
            if not rcn.status : continue

            tag = 0
            rcn.rate.update_value()
            if rcn.rate.pyformat and rcn.rate.pyformat > criteria:
                print(rcn.rate.str, rcn.rate.pyformat, criteria)
                tag = 1

            if tag:
                # check kdec species
                if len(rcn.reactants) != 1 :
                    print(len(rcn.reactants))
                    raise ValueError('RS: tau multiple reactants',criteria)

                # allow multiple or not
                if tag_multiple or len(sps_num[1][sps_num[0].index(rcn.reactants[0])]) == 1:
                    if rcn.rate.str != 'KDEC': print('Find fast reaction!: ',rcn.reactants,rcn.rate.str, criteria)
                    lumppds.append([deepcopy(rcn.products),deepcopy(rcn.ratiosPD)])
                    lumpsps.append(rcn.reactants[0])

                    # mute this reaction in output
                    rcn.status=False

                    #break # check only one kdec species

                # not in other reactions
                else: #elif rcn.reactants[0] in lumpsps: 
                    print('RS: in multiple reactions', rcn.reactants[0])

        if tag_multiple:
            for s in lumpsps:
                for i in sps_num[1][sps_num[0].index(s)]: reactions[i].status = False

    elif Type == 'jump':
        # try replacing A by B if A only forms B
        # distruction list
        dlist = destructions(reactions,Type='pro')
        nrc = len(dlist[0]) # reverse order
        for s in range(nrc): # reactants
            # dlist[0][i],dlist[1][i][0] -> A and B
            i = nrc-s-1
            if len(dlist[1][i]) == 1:
                if dlist[0][i] in primaryVOCs: continue # not jump primary voc
                if dlist[0][i]+'->'+dlist[1][i][0] in frozen: continue
                if target != [] and dlist[0][i] not in target and dlist[1][i][0] not in target: continue
                # index
                isp, jsp = sps.index(dlist[0][i]), sps.index(dlist[1][i][0])
                # keep mass balance
                if abs(species[isp].functionalGroups['C']-species[jsp].functionalGroups['C']) > 3: 
                    continue # mass conservation
                # check Psat, do not jump LVOC !! new
                if species[isp].condensed and species[isp].psat_atm < 1E-6: continue
                # check mass if less than 100, reserve for removing
                #if species[sps.index(dlist[1][i][0])].mass <= 100. or species[sps.index(dlist[0][i])].mass <= 100.: continue

                lumppds.append([[dlist[1][i][0]],[1.0]])
                lumpsps.append(dlist[0][i])
                if len(lumpsps) >= criteria: break

    elif Type == 'replace':
        # reaction list
        # try replacing A by B if A, B are reactants of the same reaction (rtB > rtA)
        nrcn = len(reactions)
        for c in range(nrcn):
            rcn = reactions[nrcn-c-1] # reversed order            
            if not rcn.status or len(rcn.products) <= 1: continue
            # order products by brt
            tmp = []
            for i in [x for _, x in sorted(zip(rcn.ratiosPD, rcn.products))]:
                if i not in NokeepSp:
                    # add limitation on mass to avoid replace too different species
                    if species[sps.index(i)].mass <= 100. : continue 
                    tmp.append(i) #replaced pds
            ntmp = len(tmp)
            if ntmp >= 1:
                for i in range(ntmp):
                    for j in range(ntmp):
                        if i == j : continue
                        isp,jsp = sps.index(tmp[i]),sps.index(tmp[j])

                        # check target
                        if target != [] and tmp[i] not in target and tmp[i] not in target: continue

                        # check type: not replace radicals
                        if species[isp].Radical:
                            if not species[jsp].Radical: continue
                        else:
                            if species[jsp].Radical: continue

                        # check mass
                        if abs(species[isp].functionalGroups['C']-species[jsp].functionalGroups['C']) > 3:
                            continue

                        if tmp[i]+'->'+tmp[j] not in frozen: # replace species A by B [A,B]
                            lumppds.append([[tmp[j]],[1.0]])
                            lumpsps.append(tmp[i])
                            if len(lumpsps) >= criteria: break
                    if len(lumpsps) >= criteria: break
            if len(lumpsps) >= criteria: break
        #print('find the replace sets.',lumpsps,lumppds)

    # modify species list

    for i,j in enumerate(lumpsps):
        tag = 1
        if j in sps:
            js = sps.index(j) # index
            tmp = 'to'
            species[js].status = False

            for k in lumppds[i][0]:
                if k in NokeepSp: continue
                if k in sps : # both in species list
                    jp = sps.index(k)
                    species[jp].jump.append(j)
                    tmp += ' '+k
                else: 
                    print('RS: not find in species list', k, lumppds[i][0])
                    tag = 0
            species[js].jump.append(tmp)
        else: 
            print('RS: not find species for jumping', j, lumppds[i][0])
            tag = 0
        if tag == 0: 
            raise RuntimeError('RS: not found both jumpped and jumpping species.')

    # merge reactions and ratio
    for rcn in reactions:
        if rcn.status:
            n=0 # change products
            while n < len(rcn.products):
                # if find species that need to replace
                if rcn.products[n] in lumpsps:

                    # replace i by pds
                    pds = lumppds[lumpsps.index(rcn.products[n])]

                    # new ratio of pds
                    rts = []
                    if rcn.ratiosPD[n] == 1.0:
                        rts = pds[1]
                    else:
                        for j in range(len(pds[1])): 
                            rts.append(pds[1][j]*rcn.ratiosPD[n])

                    # delete old i
                    rcn.products.pop(n)
                    rcn.ratiosPD.pop(n)

                    # add new elements
                    rcn.products.extend(pds[0])
                    rcn.ratiosPD.extend(rts)

                    # trim 
                    rcn.merge_duplicates()
                else: n+=1

    if Type != 'kdec': # change reactants 1 to 1
        for rcn in reactions:
            if rcn.status:
                for i in range(len(rcn.reactants)):
                    if rcn.reactants[i] in lumpsps: 
                        j = lumpsps.index(rcn.reactants[i])
                        rcn.reactants[i] = lumppds[j][0][0]
                # mute reactions if both reactants and products
                for i in lumppds:
                    if i[0][0] in rcn.reactants and i[0][0] in rcn.products:
                        rcn.status = False
    return lumpsps, lumppds

def reduce_reaction_byRemove(reactions,species,Type,val,concs,frozen=[]):
    """remove all reactions related to the children species if the parent species are not included
        frozen are reaction()"""

    # prepare for the reduction of species
    sps1=get_species_from_reactions(reactions)
    sps=[]
    for i in species : sps.append(i.name)

    # RRR: change certieria remove times and Type
    if '_' in Type:
        nctl = int(Type.split('_')[1])
        outType = Type.split('_')[0]
    else:
        nctl = 0
        outType = 'simple'

    # check frozen
    ifros = []
    if frozen != []:
        for i in frozen:
            if isinstance(i,int):
                #print 'frozen provided by the index of the reactions'
                #ifros = frozen
                ifros.append(i)
            else:
                #print 'search for the index of the frozen reactions...'
                for j in reactions:
                    #if j.status and i.isEqual(j): 
                    if i.isEqual(j):
                        ifros.append(reactions.index(j))
        # check
        if len(ifros) != len(set(ifros)):
            print(ifros)
            raise ValueError('RS: duplicates in ifros')
        if len(ifros) != len(frozen):
            #logger.debug(ifros,frozen) 
            #logger.debug('RS: len(ifros) != len(frozen)')
            pass

    # prepare for removing
    lump = []
    rm_rcn=[]
    nc = len(reactions)

    if val >= 1.0: # try all ratios
        for i in range(nc):
            # get reactions
            rcn = reactions[nc-i-1]
            if rcn.status and (nc-i-1) not in ifros:
                # check primary voc in case all removed
                if set(primaryVOCs) & set(rcn.reactants): #continue
                    isp = [s for s in rcn.reactants if s not in primaryVOCs]
                    tag = 0
                    for j in reactions:
                        if j.status and (set(primaryVOCs)&set(j.reactants)):
                            jsp = [s for s in j.reactants if s not in primaryVOCs]
                            if sorted(isp) == sorted(jsp): tag += 1 # keep initial reactions
                    if tag <= 1: continue # not remove the only valid reaction

                rcn.status = False
                lump.append(rcn.toSSH(Type=outType))
                rm_rcn.append(nc-i-1)

            if nctl and nctl == len(lump):break
    else:
        #get reactants and branch ratios 
        sprs,brts = branch_ratio(reactions, concs) # brts=[ratios,index of reactions, side reactants]
        # reverse sequences: remove reactions with small brts
        for i in range(nc):
            # get reactions
            rcn = reactions[nc-i-1]
            if not rcn.status: continue #and (nc-i-1) not in ifros: continue

            # get valid species
            # check reactant
            for s in rcn.species[0]:
                # PRAM not check brt
                if 'PRAM' in species[sps.index(s)].source and (nc-i-1) not in ifros:
                    rcn.status = False
                    lump.append(rcn.toSSH(Type=outType))
                    rm_rcn.append(nc-i-1)
                elif s in sprs:
                    # get branching ratio of all reactions react with rcn.species[0][0]
                    rts = brts[0][sprs.index(s)]
                    inds = brts[1][sprs.index(s)]
                    peers = brts[2][sprs.index(s)]
                    # check branching ratios of all reactions from rcn.species[0][0]
                    for j in range(len(rts)):
                        if rts[j] <= val and reactions[inds[j]].status and inds[j] not in rm_rcn and inds[j] not in ifros:
                            # remove
                            reactions[inds[j]].status = False
                            lump.append(reactions[inds[j]].toSSH(Type=outType))
                            rm_rcn.append(inds[j])
                            if nctl and nctl == len(lump): break # check lump times
                    
                if nctl and nctl == len(lump):break
            if nctl and nctl == len(lump):break
 
    # trim species list depending on reactions list
    sps2=get_species_from_reactions(reactions)
    for i in sps1:
        if i not in sps2 and species[sps.index(i)].organic : species[sps.index(i)].status=False

    return lump, rm_rcn

def reduce_reaction_byTree(reactions,species,Type,val):
    """remove all reactions related to the children species if the parent species are not included"""

    # get tree
    parents,children = tree(deepcopy(reactions))
    ns = len(parents)
    sps1=get_species_from_reactions(reactions)

    # get species list
    sps=[]
    for i in species : sps.append(i.name)
    nsp=len(species)

    # record parent name that need to reduce
    lump=[]

    # criteria: Type and val
    # remove species if the number of its children species is larger than val
    if Type == 'number':
        for n in range(ns):
            # criteria: number of derivied species <= val
            if len(children[n]) <= val : lump.append(parents[n])
                #addNoDuplicate(children[n],lump)

    # remove species if the minimum Psat of its children species > val
    elif Type == 'Psat_SVOC':
        for n in range(ns):
            # init
            psat = 1E99

            # mark RO2
            mRO2 = False

            # get the minimum psat of a tree
            for i in children[n] : 
                # find RO2
                if not mRO2 and species[sps.index(i)].RO2: mRO2 = sps.index(i)

                if species[sps.index(i)].psat_atm != 0.0: psat=min(psat,species[sps.index(i)].psat_atm)

            # check
            if psat == 1E99: # all the Psat is 0.0 (uncomputated) 
                print('!!! all the children Psat is zeros: ',parents[n], children[n])

            # criteria: the minimum Psat of derivied species > val
            if psat > val :
                # treatment for RO2: if found RO2, skip
                if RO2Treat == 'fixed' and (mRO2 or species[sps.index(parents[n])].RO2): 
                    print('### SVOC with RO2: ', mRO2, parents[n],species[sps.index(parents[n])].RO2,children[n])
                    continue
                lump.append(parents[n])

    # remove if its Psat < val
    elif Type == 'Psat_NVOC':
        for n in range(ns):
            # criteria: the psat of this parent species < val
            if species[sps.index(parents[n])].psat_atm < val :
                # test RO2
                if RO2Treat == 'fixed' and species[sps.index(parents[n])].RO2:
                    print('### NVOC with RO2: ',parents[n],species[sps.index(parents[n])].psat_atm)
                    continue
                lump.append(parents[n])

    # remove if its No.generation > val
    elif Type == 'gen':
        # get generation list
        spsGen,gens=generation(reactions)

        ns = len(spsGen)
        for n in range(ns):
            if gens[n] > val: lump.append(spsGen[n])

    # remove if the ref conc of all derived species is lower than val
    elif Type == 'conc':
        # criteria: val in format example: [['gas',1E-3,'<='],['aero',1E-4,'<=']]
        # val[0]: threshold for gas
        # val[1]: threshold for aerosol

        # process criteria
        criteria={'gas':0,'aero':0,'TM':0}
        for i in val.keys():
            if i in criteria.keys():
                criteria = val[i]
                if i == 'gas': cind = [0]
                elif i == 'aero': cind = [1]
                elif i == 'TM': cind=[0,1]
                print('criteria for ',i,' with value ',criteria)
                break

        # reference
        print('conc ref ',refconc_paths)
        # get all conc. [species name, [gas,aerosol], time array]
        datas = get_conc(refconc_paths)

        lump = []

        # process conc
        for data in datas:
            tmp_lump=[]

            # timestep
            nt = len(data[2])

            # compute total conc. in each timestep
            total=np.zeros(nt)
            for j in range(nsp):
                if species[j].organic and sps[j] in data[0]:
                    s = data[0].index(sps[j])
                    # compute
                    for k in cind: 
                        for i in range(nt):total[i] += data[1][s][k][i]
            # reduction
            for s in range(ns): # parents
                tag=0
                if 0: # compare for each children
                    for i in children[s] : # children
                        if i in NokeepSp or not i in data[0]: continue
                        for j in range(nt):
                            tmp = 0.0
                            for k in cind:
                                tmp += data[1][data[0].index(i)][k][j]
                            if total[j] > 0.0 and tmp/total[j] > criteria:
                                tag=1
                                break
                        if tag: break
                    if tag==0: tmp_lump.append(parents[s])
                else: # sum all conc of children
                    tmp = np.zeros(nt)
                    for i in children[s] : # children
                        if i in NokeepSp or not i in data[0]: continue
                        for j in range(nt):
                            for k in cind:
                                tmp[j] += data[1][data[0].index(i)][k][j]
                    for j in range(nt):
                        if total[j] > 0.0 and tmp[j]/total[j] > criteria:
                            tag=1
                            break
                        if tag: break
                    if tag==0: tmp_lump.append(parents[s])
            # process tmp_lump
            if lump == []: lump = tmp_lump
            else:
                tmp=[]
                for i in lump:
                    if i in tmp_lump: tmp.append(i)
                    else: i,' not fit in other condition'
                lump = tmp
            print(len(tmp_lump),len(lump))

    else:
        raise TypeError('not found type: '+Type)

    # reduce reactions
    for i in reactions:
        for j in i.reactants:
            if j in lump: i.status=False

    if Type == 'gen' or RO2Treat != 'fixed':
        for i in range(len(reactions)):
            rcn = reactions[i]
            # remove the production of RO2 species
            n = 0
            while n < len(rcn.products):
                js = rcn.products[n]
                if js in lump and species[sps.index(js)].RO2:
                    rcn.products.pop(n)
                    rcn.ratiosPD.pop(n)
                else: n += 1

    # trim species list depending on reactions list
    sps2=get_species_from_reactions(reactions)
    for i in sps1:
        if i not in sps2: species[sps.index(i)].status=False

    return lump

def SpeciesFitLump(X,Ys,sp,rc,gen,tau,criteria,lumptype=['plain']):
    """check if input two species can be lumped together
       compulsory criteria: Psat_atm (saturation vapor pressure)
                            tau (chemical lifetime
                                 careful: the input conc. of NokeepSp species (Parameters.py) matters)
       optional criteria: DU (degree of unsaturation)
                        formula (chemical formula)
                        MW (molar mass in tenth digit)
                        FGs (number and type of functional groups)
                        FGT (only the type of functional groups)"""
    # lump species:
    # Structure: formula, MW, number/type of the functional groups, OM/OC, atomic ratios
    # Partitioning: Psat, Hvap
    # Chemistry: lifetime tau, number of generation, degree of unsaturation, degradation rate, tree/reversed tree, reactants
    # Results: concentrations analysis

    # (1) prepare for lumping
    # record species names
    lumpsps=[X.name]

    # record the estimated parameter ratios in lumping
    # currently depending on the conc in the original case
    lumpratio=[rc[1][rc[0].index(X.name)]]

    # the generation number of X
    igen = 0
    # the new surrogate species
    # currently have the same properties as the lumped species with the largest lumpratio
    newsp = X

    # (2) search for species that can be lumped with X
    for s in Ys:

        # checking basic: not the same species X
        if s == X.name: continue

        # non-zero ref conc.
        #if rc[1][rc[0].index(s)] <= 0.0: continue

        # get s in species list
        i = sp[1][sp[0].index(s)]

        # check status
        if not i.status: continue

        # check PRAM species
        if 'PRAM' in i.source:
            if 'PRAM' not in newsp.source: continue
        elif 'PRAM' in newsp.source: continue

        # check species type
        if newsp.Radical:

            if not i.Radical: continue # radical and non-radical
            else: # check radical type: lump RO3 with RO3, ROO with ROO, RO with RO, RO2 with RO2
                if newsp.name[-2:] in ['OO','O3','O2'] and newsp.name[-2:] != i.name[-2:]: continue
        else:

            if i.Radical: continue
            else: # both non-radical
                tag = 0
                for j in ['PAN','OOH','CO3H','CO2H']: # functional group indicated by mcm name
                    if j in CheckStructureKeys:
                        if j in newsp.name and not j in i.name: tag = 1
                        elif not j in newsp.name and j in i.name: tag = 1
                if tag: continue

                # check Psat
                if CheckPsat: # not both NVOC, check the scale of log
                    if abs(math.log10(newsp.psat_atm) - math.log10(i.psat_atm)) > CheckPsat: continue

                if 'C=C' in CheckStructureKeys:
                    if not i.SOAPStructure or not newsp.SOAPStructure:
                        print(i.SOAPStructure, newsp.SOAPStructure)
                        raise AttributeError('RS: SOAPStructure is not in both i and newsp.')
                    tmp0 = 0 # for newsp
                    tmp1 = 0 # for i
                    for j in [16,17,18,19,20]: # C=C
                        if j in i.SOAPStructure[0]: tmp1 += i.SOAPStructure[1][i.SOAPStructure[0].index(j)]
                        if j in newsp.SOAPStructure[0]: tmp0 += newsp.SOAPStructure[1][newsp.SOAPStructure[0].index(j)]
                    if tmp0 + tmp1 != 0 and tmp0 * tmp1 == 0: continue

        # for both radicals and non-radicals
        # check the minimum number of carbon
        if min(newsp.functionalGroups['C'],i.functionalGroups['C']) <= 3 : continue 

        # check mass diff
        if abs(newsp.mass - i.mass) > 100.: continue
        
        if 'C' in CheckStructureKeys:
            if abs(newsp.functionalGroups['C'] - i.functionalGroups['C']) > CheckStructure['C'] : 
                continue

        # check nitrate        
        if 'N/C' in CheckStructureKeys:
            if newsp.ratios['N/C'] * i.ratios['N/C'] == 0.0 and newsp.ratios['N/C'] + i.ratios['N/C'] != 0.0: continue

        # check oxygen  ! only for API    
        if 'O' in CheckStructureKeys:
            if  abs(newsp.functionalGroups['O'] - i.functionalGroups['O']) > CheckStructure['O'] : continue

        # generation
        if CheckGen and not isFitScale(gen[1][gen[0].index(i.name)],igen,CheckGen):
            continue

        # lifetime
        if CheckTau and not ifFitLifetime(X,i,tau,CheckTau): continue

        # structure and partitioning
        tag = 0
        for k in ['DU', 'formula', 'MW']:
            if k in CheckStructureKeys and ifFitStructure(newsp,i,k,CheckStructure[k]): tag = 1
        if tag: continue

        ### fits all criteria, now record i
        newsp = i
        lumpsps.append(i.name)
        lumpratio.append(rc[1][rc[0].index(i.name)])

        break #only merge once

    # check if find at least one species for lumping
    if len(lumpsps) > 1:
        if len(lumpratio) != 2: raise ValueError('len lumpratio != 2!')
        # compute lumpratios
        lumpratios = []
        for ltype in lumptype:
            if ltype == 'plain': rts = lumpratio
            elif ltype == 'tau': # change ratio computed by tau * Conc.
                tmp,rts = 0.0, [0.,0.]
                for i in range(2):
                    tmp += tau[1][tau[0].index(lumpsps[i])]*lumpratio[1-i]
                for i in range(2): 
                    if tmp != 0.0: 
                        rts[i] = tau[1][tau[0].index(lumpsps[1-i])]*lumpratio[i]/tmp
                    else: 
                        rts[i] = 0.0
            elif isfloat(ltype): # ratio of A
                rts = [float(ltype),100.-float(ltype)]
            rts=array_to_ratio(rts)
            lumpratios.append(rts)

        return [lumpsps,lumpratios]
    else: return False

def ifFitStructure(X,Y,criteria,val):
    """types: DU, formula, MW"""
    if criteria == 'DU':
        val = 2.
        if isFitScale(X.DU,Y.DU,val): #if X.DU == Y.DU:
            return True
        else:
            return False

    elif criteria == 'formula':
        if X.formula == Y.formula:
            return True
        else:
            return False

    elif criteria == 'MW':
        val = 100.
        #if isFitScale(X.mass,Y.mass,val):
        if abs(X.mass - Y.mass) <= val:
            return True
        else:
            return False
    else:
        raise TypeError('criteria for lumping is not recognized: '+ criteria)


def ifFitLifetime(X,Y,tau,val = 3.):
    """check if the ratio of two inputs is in the range of val [1/val to val]"""
    if X.name in tau[0] and Y.name in tau[0]:
        t1 = max(tau[1][tau[0].index(X.name)],1.0)
        t2 = max(tau[1][tau[0].index(Y.name)],1.0)

        if isFitScale(t1,t2,val):
            #logger.debug('lifetime tuned and fits', X.name, t1, Y.name, t2,t1/t2)
            return True
        else:
            pass
    return False

def LumpReaction(rcn,newsp,lumpsps,lumpratio,Type = 'lump'):

    if not rcn.status: return False

    # reactants
    # get if in lumping
    tmp = list(set([x for x in rcn.reactants if x in lumpsps]))

    if len(tmp) == 1:
        # compute the estimated parameter ratio
        rt = lumpratio[lumpsps.index(tmp[0])]
        if Type == 'lump':
            if rt == 0.0:
                rcn.status = False
            elif rt < 0.0 or rt > 1.0:
                print(rt)
                raise ValueError('RS: wrong rt.',rt,newsp,lumpsps,lumpratio,rcn.toSSH())
            elif rt != 1.0: # use
                rcn.rate.str += ('*{:6.4E}'.format(rt))
                rcn.rate.update()
        elif Type == 'replace': # replace the smallesrt by the largerst one
            if tmp[0] != newsp:
                # remove reactions whose reactants is lumped species
                rcn.status = False
                print('replace' ,tmp[0],' by ',newsp)
        elif Type == 'mix': # mix removing and replace
            rRT = 0.001
            if rt <= rRT:
                print('rt <= replacing ratio !!!!',rRT,rt,newsp,lumpsps,lumpratio)
                rcn.status = False
            elif rt > rRT : # both ratios > rRT: lumping
                if (1-rt) > rRT:
                    rcn.rate.str += ('*{:6.4E}'.format(rt))
                    rcn.rate.update()
                #else: # only one larger than rRT, replaceing
        else:
            raise TypeError('RS: LumpReaction Type is not recognized, '+Type)

        # change name
        rcn.reactants[rcn.reactants.index(tmp[0])] = newsp

    elif len(tmp) > 1: # find more than one lumped species as reactants in the same reaction
        tmp = list(set([x for x in rcn.reactants if x not in lumpsps]))
        #if tmp == []: rcn.status = False
        raise ValueError('find two lumped species as reactants in the same reaction',lumpsps,rcn.toSSH())

    # products
    # replace name if any products in lumpsps
    for s in list(set([x for x in rcn.products if x in lumpsps])):
        rcn.products[rcn.products.index(s)] = newsp

    # trim
    rcn.merge_duplicates()

def trim_scheme(reactions,species):
    """check the scheme and remove unused ones"""
    dlist = productions(reactions,'rea')

    sps = []
    for i in range(len(dlist[0])):
        n = 0
        while n < len(dlist[1][i]) and not (set(primaryVOCs)&set(dlist[1][i])): # extend
            s = dlist[1][i][n]
            if s in dlist[0]:
                for k in dlist[1][dlist[0].index(s)]:
                    if k not in dlist[1][i]: dlist[1][i].append(k)
            n += 1
        if not (set(primaryVOCs) & set(dlist[1][i])): # find species can not track back primary voc
            #logger.debug('RS trim_scheme: found isolated species, ',dlist[0][i])
            sps.append(dlist[0][i])

    # remove species
    for i in species:
        if i.status and i.organic and i.name not in primaryVOCs:
            if i.name not in dlist[0] and i.name not in sps: sps.append(i.name)  
            if i.name in sps: i.status = False
    if sps != []: 
        #logger.debug('RS trim_scheme: ',sps)
        pass

    # remove reactions
    for r in range(len(reactions)):
        for i in reactions[r].reactants:
            if i in sps:
                reactions[r].status = False

    return reactions, species, sps
        
