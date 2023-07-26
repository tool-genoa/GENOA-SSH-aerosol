# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  ReductionStrategy.py contains functions that reduces reactions/species 
#
#   based on different reduction strategies.
#
# ================================================================================

import math
import numpy as np
from copy import deepcopy 
from itertools import combinations #permutations

import Module as md
from Parameters import primaryVOCs, Roptions, NokeepSp
from SSHResultProcess import get_conc, ref_conc
from KineticMCMtoSSH import get_type_rate
from Functions import duplicates, list_all_fit, \
                      isfloat, compare, get_negative_str, \
                      isFitScale, array_to_ratio
from ChemRelation import get_species_from_reactions, lifetime, \
                         productions, destructions, \
                         branch_ratio, tree, generation

# if add species in the output, add both in lump_props and read_props
# output of lumped species
lump_props = ['New Species','name','lumpRatio','lumps','formula','SMILES','mass',\
              'RO2','condensed','groupID','source','psat_atm','dHvap_KJ','rConcs','lifetime']
# outputs read from the Species class
read_props = ['name','formula','SMILES','mass','RO2','condensed','groupID','source','psat_atm','dHvap_KJ'] # read from species

groupIDs = [['APINENE', 'BPINENE', 'LIMONENE'],['APINENE'],['BPINENE'],['LIMONENE'],['APINENE', 'BPINENE'],['BPINENE', 'LIMONENE'],['APINENE', 'LIMONENE']]
# not use group id
tag_source = False

# read options
print('\nRead settings for lumping...')
if 'OrderType'in Roptions.keys() and Roptions['OrderType'] in ['conc', 'reaction']:
    OrderType = Roptions['OrderType']
    #print('\tRS: read OrderType from Roptions: ', OrderType)
else: OrderType = 'reaction'

# type of candidates lumped species for A
if 'PeerType'in Roptions.keys() and Roptions['PeerType'] in ['reaction', 'all', 'psat']: 
    PeerType = Roptions['PeerType']
    #print('\tRS: read PeerType from Roptions: ', PeerType)
else:  PeerType = 'all'

# if A and candidates can produce one or the other
if 'CheckLink'in Roptions.keys() and isfloat(Roptions['CheckLink']): 
    CheckLink = Roptions['CheckLink']
    #print('\tRS: read CheckLink from Roptions: ', CheckLink)
else: CheckLink = False # default

# treatment for RO2 species
if 'RO2Treat'in Roptions.keys(): RO2Treat = Roptions['RO2Treat']
else: RO2Treat = None 
#if RO2Treat: print('\tRS: read RO2Treat from Roptions: ', RO2Treat)
# (fixed is for remove by Psat_S/NVOC; remove RO2 after lumping)

# if check lifetime, No.generation, saturation vapor pressure, stracture
if 'CheckTau'in Roptions.keys(): CheckTau = Roptions['CheckTau'] # False, value of tolerance scale
else: CheckTau = False
if CheckTau: print('\tRS: read CheckTau from Roptions: ', CheckTau)

if 'CheckGen'in Roptions.keys(): CheckGen = Roptions['CheckGen']
else: CheckGen = False
if CheckGen: print('\tRS: read CheckGen from Roptions: ', CheckGen)

if 'CheckPsat'in Roptions.keys(): CheckPsat = Roptions['CheckPsat'] # False, value of tolerance scale
else: CheckPsat = False
if CheckPsat: print('\tRS: read CheckPsat from Roptions: ', CheckPsat)

if 'CheckStructure'in Roptions.keys() and isinstance(Roptions['CheckStructure'],dict):
    tmp = [i for i in Roptions['CheckStructure'] if Roptions['CheckStructure'][i]]
    CheckTypeKeys = [i for i in ['CO3H','CO2H','PAN','OOH'] if i in tmp]
    CheckStructureKeys = [i for i in ['C=C','C','O','N/C'] if i in tmp] # ['DU', 'formula', 'MW']
    CheckStructure = {}
    for i in CheckStructureKeys: CheckStructure[i] = Roptions['CheckStructure'][i]
    print('\tRS: read CheckStructure from Roptions: ', CheckStructure, ', Check type keys: ', CheckTypeKeys)
    for i in tmp:
        if i not in CheckTypeKeys and i not in CheckStructureKeys:
            raise TypeError('\tRS: CheckStructure type not recognized.',i,Roptions['CheckStructure'])

else: CheckStructure, CheckStructureKeys, CheckTypeKeys = [], [], []

def reduction_chemistry_scheme(reactions,species,Type,criteria, concs_ave, concs_min, refconc_paths, RefConcRead):
    """the main function used for reduction:
        reduce reactions and species depending on the input reduction type and criteria"""
    # reduction keywords
    # bratio, lump
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
        
    elif Type == 'Psat_aero':
        lsps = []
        # pre-reduction, remove aerosols with Psat >= 1E-3 atm
        for i in species:
            if i.status and i.condensed and i.psat_atm >= criteria:
                print('Psat_aero: remove aerosol: ',i.name, i.psat_atm)
                i.condensed = False
                lsps.append(i.name)
        return lsps

    elif Type == 'bratio':
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
                           target = [], lumptype = ['tau'],
                           rcSpsCut = None):
    """prepare species and reactions depending on the given lump species list"""

    # control the maximumlump times
    if len(criteria.split('_')) == 1: 
        tag_ctl= int(criteria)#0 # no limitation
    else:
        tag_ctl = int(criteria.split('_')[1])
        criteria = criteria.split('_')[0]

    # get species name
    sps=[i.name for i in species]

    # get names in reactions
    rcSps = get_species_from_reactions(reactions)
    # cut rcSps into several piecies
    if isinstance(rcSpsCut,list) and len(rcSpsCut)==2: # [i,n] i = 1,2,...,n
        nsps = len(rcSps)
        n = int(nsps/rcSpsCut[1]) # interval
        i = rcSpsCut[0] # index + 1
        if n == 0:
            if i <= nsps: rcSps = [rcSps[i-1]]
            else: rcSps = []
        else:
            if i == rcSpsCut[1]: # last bin
                rcSps = rcSps[n*(i-1):]
            else:
                rcSps = rcSps[n*(i-1):n*i]
        if target != [] and rcSps != []: 
            for i in target:
                if i not in rcSps: rcSps.append(i)

    # remove primary VOCs
    for i in primaryVOCs:
        if i in rcSps: rcSps.remove(i)

    # remove other non-lumped species
    phos = []
    for i in reactions:
        if 'Tab' in i.rate.str or 'KINEITC' in i.rate.str:
            for j in i.reactants:
                if species[sps.index(j)].organic and j not in phos: phos.append(j)
    for i in phos:
        if i in rcSps: rcSps.remove(i)

 
    # check if need stop
    if rcSps == []: return [],[],[]
    
    # only keep organic species in rcSps and rcLump
    rcLump=[] #[name, lumped species]
    for i in rcSps:
        s = species[sps.index(i)]
        # check if i is actived
        if not s.status:
            raise AttributeError('RS: non active species is in the reaction list, species '+s.name)
        if s.organic: # record lump in rcLump
            rcLump.append([s.name,s.lump+s.jump])

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
    #if SpOrderInd == []: print('RS: SpOrderInd is empty')

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
        if sname in frozen: continue
        if sname in phos: continue

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
        #if msp.Radical and msp.name[-1] == 'O': continue

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

            #SpList1 = []
            #for ips in rcSps:
            #    if species[sps.index(ips)].status:
            #        if CheckLink and (ips in plist_now or ips in dlist_now): continue
            #        elif ips not in SpList1: SpList1.append(ips)
        # check
        if frozen != []:
            # check and remove frozen pairs
            for i in frozen:
                if sname in i: # if sname in the list
                    if len(i) == 1:
                        SpList = []
                    for j in i:
                        if j in SpList: SpList.remove(j)

        for i in phos:
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
            # properties in order: average MWs, Psat, dHvap_KJ, functional groups, henry, source # gamma
            tmp = [0.0,1,0.0,[0.0]*60, {}, 1., []]
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
                    for k in range(60):
                        if k in i.SOAPStructure[0] :
                            tmp[3][k]+=(i.SOAPStructure[1][i.SOAPStructure[0].index(k)])*rt
                for k in i.functionalGroups.keys(): 
                    if k not in tmp[4].keys(): tmp[4][k] = i.functionalGroups[k] * rt
                    else: tmp[4][k] += i.functionalGroups[k] * rt
                tmp[5] *= pow(i.henry,rt) #henry
                if tag_source:
                    for k in i.source:
                        if k not in tmp[6]: tmp[6].append(k)

            # assign new values
            i = lp_sps[sps.index(newsp)]
            # copy newsp for output
            oldsp = deepcopy(i)
            i.mass = tmp[0]
            i.psat_atm = tmp[1]
            i.psat_torr = i.psat_atm*760.
            i.dHvap_KJ = tmp[2]
            i.SOAPStructure = [[],[]]
            i.henry = tmp[5]

            if tag_source:
                # update source and group if needs
                tmp[6] = list(sorted(tmp[6]))
                if i.source != tmp[6]:
                    i.source = tmp[6] # update
                    if i.source in groupIDs:
                        i.groupID = groupIDs.index(i.source) + 1
                    elif i.source != []:
                        print('i.source not in groupIDs!', i.source, i.name)

            #i.gamma = tmp[6]
            for k in range(60):
                if tmp[3][k] != 0.0: 
                    i.SOAPStructure[0].append(k)
                    i.SOAPStructure[1].append(tmp[3][k])
            i.functionalGroups = tmp[4]
            for k in ['C','H','N','O']:
                if k not in i.functionalGroups.keys(): 
                    print('??? lp!!! functionalGroups miss key!!!',k,i.name.lump)
                    i.functionalGroups[k] = 0.
            if i.functionalGroups['C'] != 0.:
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

            # read properties of all species
            for i in range(len(lumpsps)+1):
                if i == 0: 
                    j = lp_sps[sps.index(newsp)] # new sps
                    isps = newsp
                else: # ! change index to index - 1
                    isps = lumpsps[i-1] # name
                    if isps == newsp: j = oldsp
                    else: j = lp_sps[sps.index(isps)] # species in lp

                # ratio
                s = lump_props.index('lumpRatio')
                if i == 0: lump_prop[s] +=  '\t{:}\t{:}\t'.format(lumptype[n],sum(lumpratio))
                else: lump_prop[s] += '\t{:}'.format(lumpratio[i-1])
                # record properties
                for k in read_props:
                    s = lump_props.index(k)
                    lump_prop[s] += '\t{:}'.format(eval('j.'+k))
                # 'rConcs' [rcSps,rcConcs]
                s = lump_props.index('rConcs')
                if i == 0 or isps not in rcSps:
                    lump_prop[s] += '\t-'
                else: 
                    lump_prop[s] += '\t{:}'.format(rcConcs[rcSps.index(isps)])

                # 'lifetime' [tau_sps,taus]
                s = lump_props.index('lifetime')
                if i == 0 or isps not in tau_sps:
                    lump_prop[s] += '\t-'
                else: 
                    lump_prop[s] += '\t{:}'.format(taus[tau_sps.index(isps)])

            lumpinfo.append(lump_prop)

    return newsps,lumpinfo,changes

def reduce_reaction_bySpecies( reactions,species,Type,
                               criteria,frozen = [],target = [],
                               cut_order = None):
    """
        Inputs:
            remove species according certain ceriteria
            reactions: reaction list
            species: species list
            type: what type of reduction
            criteria: how many reudctions are wannted
            frozen: if certain combination/species is not wanted to be removed
            target: if only want to conduct relavant reduction
            cur_order: cut possibility for multiprocessing. if 1: from large to small; if 2: from small to large
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
                else: print('remove kdec: ',rcn.reactants)

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

        sps_num = [[],[],[]] # reactant, index in reaction
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

            tag = 0 # check if it is a fast reaction
            rcn.rate.update_value()
            if rcn.rate.pyformat and rcn.rate.pyformat > criteria:
                #print(rcn.rate.str, rcn.rate.pyformat, criteria)
                tag = 1

            if tag:
                # check kdec species
                if len(rcn.reactants) != 1 : continue
                    #print(len(rcn.reactants))
                    #raise ValueError('RS: tau multiple reactants',criteria)

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
        #dlist = destructions(reactions,Type='pro')
        #plist = productions(reactions,Type='pro')
        dlist, plist = [[],[],[],[]] ,[[],[],[],[]] # sps, prd, rates, rts
        for i in reactions:
            if not i.status: continue
            rate, rt = get_type_rate(i.rate.SSH)
            for j in i.reactants:
                if j not in NokeepSp:
                    if j not in dlist[0]: # check
                        dlist[0].append(j)
                        dlist[1].append([deepcopy(i.products)])
                        dlist[2].append([rate])
                        dlist[3].append([[rt]*len(i.products)])
                    else: # add new species
                        s = dlist[0].index(j)
                        # check rate
                        if rate in dlist[2][s]:
                            k = dlist[2][s].index(rate)
                            for l in i.products:
                                if l in dlist[1][s][k]: # find pd, update rt
                                    dlist[3][s][k][dlist[1][s][k].index(l)]+=rt
                                else: # add new pd and ratio
                                    dlist[1][s][k].append(l)
                                    dlist[3][s][k].append(rt)
                        else: # add new rate, pds, rts 
                            dlist[2][s].append(rate)
                            dlist[1][s].append(deepcopy(i.products))
                            dlist[3][s].append([rt]*len(i.products))

            for j in i.products:
                if j not in NokeepSp:
                    if j not in plist[0]: 
                        plist[0].append(j)
                        plist[1].append([deepcopy(i.reactants)])
                        plist[2].append([rate])
                        plist[3].append([[rt]*len(i.reactants)])
                    else: # find in list, check rate
                        s = plist[0].index(j)
                        if rate in plist[2][s]:
                            k = plist[2][s].index(rate)
                            for l in i.reactants:
                                if l in plist[1][s][k]: # find pd, update rt
                                    plist[3][s][k][plist[1][s][k].index(l)]+=rt
                                else: # add new pd and ratio
                                    plist[1][s][k].append(l)
                                    plist[3][s][k].append(rt)
                        else: # add new rate, pds, rts 
                            plist[2][s].append(rate)
                            plist[1][s].append(deepcopy(i.reactants))
                            plist[3][s].append([rt]*len(i.reactants))


        # order species and rts
        for i in range(len(dlist)):
            for j in range(len(dlist[1][i])): # order by ratios
                if len(dlist[3][i][j]) > 1:
                    dlist[3][i][j],dlist[1][i][j] = zip(*sorted(zip(dlist[3][i][j],dlist[1][i][j])))
        for i in range(len(plist)):
            for j in range(len(plist[1][i])): # order by ratios
                if len(plist[3][i][j]) > 1:
                    plist[3][i][j],plist[1][i][j] = zip(*sorted(zip(plist[3][i][j],plist[1][i][j])))

        n = 0 # jumping times
        for i in reversed(range(len(sps))):
            isp = sps[i] # name
            if not species[i].organic: continue
            elif isp in primaryVOCs: continue # not jump primary voc
            elif target != [] and isp not in target: continue
            
            # find in the destruction list
            if isp in dlist[0]: idl = dlist[1][dlist[0].index(isp)]
            else: idl = []
            # find in the production list
            if isp in plist[0]: ipl = plist[1][plist[0].index(isp)]
            else: ipl = []

            # check len
            if len(idl) == 1 or len(ipl) == 1:
                isps = []
                
                if len(idl) == 1: # A1,A2,A3 -> B or A-> B
                    if cut_order == 1: 
                        #print('Reject idl cuz cut_order',cut_order,idl,frozen)
                        pass
                    else:
                        isps,dplist,jsps,jrts = [isp], dlist, idl[0],dlist[3][dlist[0].index(isp)][0]
                        #if len(jsps) != len(jrts): raise ValueError('Accepted idl',cut_order,idl,frozen,isps,jsps,jrts)
                        
                if isps == [] and len(ipl) == 1: # not one to more / more to one
                    if cut_order == 2:
                        #print('Reject ipl cuz cut_order',cut_order,ipl,frozen)
                        pass
                    else:
                        #print('Accepted ipl',cut_order,idl,frozen)
                        isps, dplist,jsps,jrts = [isp], plist, [],[]  # A -> B1,B2,B3
                        # build jsps, jrts
                        for j in ipl[0]: # remove NokeepSp
                            if j not in NokeepSp:
                                jsps.append(j)
                                jrts.append(dplist[3][dplist[0].index(isp)][0][ipl[0].index(j)])

                if isps == []: continue

                if 0: # check all
                    # may fit jumping, check other properties
                    tag = 0 # tag = 1 not jump
                    rate = dplist[2][dplist[0].index(isp)][0] # B name
                    # find other A
                    for j in range(len(dplist[0])):
                        if jsps in dplist[1][j]: # find other A
                            if len(dplist[1][j]) == 1 and dplist[2][j][0] == rate:
                                if dplist[0][j] not in isps: 
                                    isps.append(dplist[0][j])
                                    jrts.append(dplist[3][j][dplist[1][j].index(jsps)])
                            else: # find multi B, not jump
                                tag += 1
                    if tag: # jumping species is involved in other reactions
                        print(tag, 'involved in other reactions', jsps, isps) 
                        continue
                    else: isps = list(sorted(isps)) # sort
                    jrts = array_to_ratio(np.average(jrts,axis=0))
                else: 
                    jrts = array_to_ratio(jrts)

                # check frozen species & relationship in jumping
                for j in frozen:
                    # check individual species
                    if '->' not in j:
                        if j in isps+jsps: continue
                # check relationship
                itmp,jtmp = '',''
                for j in isps: itmp+=(j+',')
                for j in jsps: jtmp+=(j+',')
                itmp,jtmp = itmp[:-1],jtmp[:-1]
                if len(idl) == 1 : # 'A1 A2 A3 -> B' or 'A1 -> B' 
                    tmp = '{:s}->{:s}'.format(itmp,jtmp)
                elif len(ipl) == 1: # 'A -> B1 B2 B3'
                    tmp = '{:s}->{:s}'.format(jtmp,itmp)
                if tmp in frozen: continue
                #else: print('find new jp', tmp, frozen)
                
                # check
                if len(isps) < 1 or len(jsps) < 1:
                    print('jump: isps or jsps len < 1. continue',j,isps,jsps,jrts,tmp)
                    continue
                elif set(isps) & set(jsps):
                    raise ValueError('jump: found the same species in both isps and jsps.',j,isps,jsps,jrts,tmp)
                elif len(jsps) != len(jrts):
                    raise ValueError('jump: len not match.',j,isps,jsps,jrts,tmp)

                # check properties of species for jumping
                imass = np.average([species[sps.index(j)].mass for j in isps])
                if len(idl) == 1 : jmass = np.sum([species[sps.index(jsps[j])].mass*jrts[j] for j in range(len(jsps))])
                elif len(ipl) == 1: jmass = np.average([species[sps.index(j)].mass for j in jsps])
                if np.isnan(imass) or np.isnan(jmass):
                    print([species[sps.index(j)].mass for j in isps],[species[sps.index(jsps[j])].mass*jrts[j] for j in range(len(jsps))],[species[sps.index(j)].mass for j in jsps])
                    raise ValueError('jump: sum of mass is nan.',j,isps,jsps,jrts,tmp,imass,jmass)
                elif abs(imass-jmass)/(imass+jmass)*2 > 0.5:
                    print('jump refuse: mass diff',imass,jmass,abs(imass-jmass)/(imass+jmass)*2)
                    continue

                # record all possibilities
                if len(idl) == 1 :
                    for j in isps:
                        lumpsps.append(j)
                        lumppds.append([jsps,jrts]) # species, ratios
                        #print('jp. idl == 1 ',j,[jsps,jrts])
                elif len(ipl) == 1:
                    for j in jsps:
                        lumpsps.append(j)
                        lumppds.append([isps,[1.0]*len(isps)])
                        #print('jp. ipl == 1 ',j,[isps,[1.0]*len(isps)])
                n += 1
                # check out
                #if not (isinstance(cut,list) and len(cut)) == 2:
                #    if n >= criteria: break
                if n >= criteria: break

    elif Type == 'replace':

        # reaction list
        # try replacing A by B if A, B are reactants of the same reaction (rtB > rtA)
        plist = [[],[],[],[]] # rate, reactants, products, ratios
        for i in range(len(reactions)):
            if not reactions[i].status: continue
            rcn = reactions[i]
            # reactants
            rea = list(sorted(rcn.reactants))
            # rate
            rate, rt = get_type_rate(rcn.rate.SSH)
            if rea in plist[1]: 
                # check reactants, get index
                j = plist[1].index(rea)
                # check rate
                if rate in plist[0][j]:
                    k = plist[0][j].index(rate)
                    for l in rcn.products:
                        if l in plist[2][j][k]:
                            # update ratio
                            plist[3][j][k][plist[2][j][k].index(l)]+=rt
                        else: # add new species and ratio
                            plist[2][j][k].append(l)
                            plist[3][j][k].append(rt)
                else: # save new rate
                    plist[0][j].append(rate)
                    plist[2][j].append([k for k in rcn.products if k not in NokeepSp])
                    plist[3][j].append([rt]*len(plist[2][j][-1]))

            else: # new
                plist[0].append([rate]) # add rate
                plist[1].append(rea) # add reactants
                plist[2].append([[k for k in rcn.products if k not in NokeepSp]])
                plist[3].append([[rt]*len(plist[2][-1][0])]) # add ratios for pbs
                
        for l in reversed(range(len(plist[1]))):
            for i in range(len(plist[2][l])):
                s = plist[2][l][i]
                if len(s) <= 1: continue # at least two products
                # check target
                elif target != [] and not (set(target) & set(s)): continue
                # order by mass
                tmp = [x for _,x in sorted(zip([species[sps.index(j)].mass for j in s],s))]
                tmp1 = [] # used species

                for j in list(combinations(tmp,2)): # all possibilities
                    if (set(target) & set(tmp1)): continue
                    elif target != [] and not (set(target) & set(j)): continue
                    isp,jsp = sps.index(j[0]),sps.index(j[1]) # get A -> B
                    
                    # check species properties for replacing
                    # need to be same types
                    if not species[isp].organic or not species[jsp].organic: continue
                    elif species[isp].Radical != species[jsp].Radical: continue
                    elif  species[isp].condensed != species[jsp].condensed: continue

                    # check mass
                    k = abs(species[isp].mass-species[jsp].mass)/(species[isp].mass+species[jsp].mass)*2.
                    if k > 0.5:
                        if cut_order != 2: print('replace: mass ratio > 0.5',j,species[isp].mass,species[jsp].mass,k)
                        continue

                    # carbon diff can not be huge
                    if abs(species[isp].functionalGroups['C']-species[jsp].functionalGroups['C']) > 3:
                        if cut_order != 2: print('replace refused: C number diff.',j,species[isp].functionalGroups['C']-species[jsp].functionalGroups['C'],species[isp].functionalGroups['C'],species[jsp].functionalGroups['C'])
                        continue
                    #elif abs(species[isp].mass-species[jsp].mass) >= 100.:continue

                    # last step: check frozen + set isp,jsp by species name
                    tag = 0
                    # check species in replacement
                    for k in frozen:
                        if '->' not in k:
                            # check species
                            if k in j: continue
                            
                    # first check
                    if sps[jsp]+'->'+sps[isp] not in frozen:
                        if cut_order == 1: 
                            #print('Reject j->i cuz cut_order',cut_order,sps[jsp],sps[isp],frozen)
                            pass
                        else:
                            #print('Accepted j->i',cut_order,sps[jsp],sps[isp],frozen)
                            isp,jsp,tag = sps[jsp],sps[isp],1
                    # second check
                    if tag != 1 and sps[isp]+'->'+sps[jsp] not in frozen:
                        if cut_order == 2: 
                            #print('Reject i->j cuz cut_order',cut_order,sps[jsp],sps[isp],frozen)
                            pass
                        else:
                            #print('Accepted j->i',cut_order,sps[jsp],sps[isp],frozen)
                            isp,jsp,tag = sps[isp],sps[jsp],1

                    if tag != 1: continue # check next set

                    # replace species A by B [A,B]
                    lumpsps.append(isp)
                    lumppds.append([[jsp],[1.0]])

                    # save replaced species to prevent using twice
                    for k in isp,jsp: 
                        if k not in tmp1: tmp1.append(k)

                # check out condition
                    if len(lumpsps) >= criteria: break
                if len(lumpsps) >= criteria: break
            if len(lumpsps) >= criteria: break
        #print('find the replace sets.',lumpsps,lumppds)


    # modify species list
    for i,j in enumerate(lumpsps):
        tag = 1
        if j in sps:
            js = sps.index(j) # index
            tmp = '->'
            species[js].status = False

            # get all jumped species
            tmpjp = []
            for k in range(len(lumpsps)):
                if lumppds[k]==lumppds[i]: tmpjp.append(lumpsps[k])
            if tag_source:
                srcs = []
                for k in tmpjp:
                    if k in sps:
                        for l in species[sps.index(j)].source:
                            if l not in srcs: srcs.append(l)
                srcs = list(sorted(srcs))
                        
            for k in lumppds[i][0]:
                if k in NokeepSp: continue
                if k in sps : # both in species list
                    jp = sps.index(k) # get index
                    if not species[jp].organic: continue # check organic
                    species[jp].jump += tmpjp
                    tmp += (k+',')
                    
                    if tag_source:
                    # update source
                        if srcs != [] and species[jp].source != srcs:
                            isrc = 0
                            for l in srcs: 
                                if l not in species[jp].source: 
                                    species[jp].source.append(l)
                                    isrc += 1
                            # update group ID  
                            if isrc:
                                species[jp].source = list(sorted(species[jp].source))
                                if species[jp].source in groupIDs:
                                    species[jp].groupID = groupIDs.index(species[jp].source) + 1
                                #else:
                                #    print('species[jp].source not in groupIDs!', species[jp].source, species[jp].name)
                else: 
                    print('RS: not find in species list', k, lumppds[i][0])
                    tag = 0
            # add info for jumping
            if tmp[-1] == ',': tmp = tmp[:-1]
            species[js].jump.append(tmp)
        else: 
            print('RS: not find species for jumping', j, lumppds[i][0])
            tag = 0
        if tag == 0: 
            raise RuntimeError('RS: not found both jumpped and jumpping species.')

    # merge reactions and ratio
    if lumpsps != []:
        for rcn in reactions:
            if rcn.status:
                n=0 # change products
                while n < len(rcn.products):
                    # if find species that need to replace
                    i = rcn.products[n]
                    if i in lumpsps:

                        # replace i by pds
                        pds = lumppds[lumpsps.index(i)]
                        #if len(pds[0]) != len(pds[1]): raise ValueError(i,pds,Type)
                        # new ratio of pds
                        rts = [j*rcn.ratiosPD[n] for j in pds[1]]

                        # delete old i
                        rcn.products.pop(n)
                        rcn.ratiosPD.pop(n)

                        # add new elements
                        rcn.products.extend(pds[0])
                        rcn.ratiosPD.extend(rts)

                        # trim 
                        rcn.merge_duplicates()
                    else: n+=1

        # change reactants by products
        n = 0 # index
        while n < len(reactions):
            rcn = reactions[n]
            if rcn.status:
                # get reactants and products
                irt, ipd = rcn.reactants, rcn.products
                # get sp need to be lumped
                tmp = list(set(irt) & set(lumpsps))
                if tmp == []: # not find lump sps 
                    n += 1
                else:
                    if len(tmp) > 1: 
                        raise ValueError('RS: one rcn should not have two reactants to be replaced.',Type,tmp, irt,ipd)
                    # get index in reactants
                    j = irt.index(tmp[0])

                    # get the replacement
                    s = lumppds[lumpsps.index(tmp[0])][0]
                    tmp1 = [k for k in s if species[sps.index(k)].organic]
                    
                    i = 0 # number of active reaction
                    for l,k in enumerate(tmp1):
                        if i > 0: # need insert rcn
                            reactions.insert(n+i,deepcopy(rcn))
                            #print('RS: insert new rcn. ',Type,i,n+i,j,k,ipd,irt,s,lumpsps,lumppds)
                        reactions[n+i].reactants[j] = k
                        i += 1

                    if i == 0:  # mute cuz not replacement
                        reactions[n].status = False
                        n += 1
                    else: # count insert number
                        n += i 
            else: n += 1
            
        # check reaction list
        for rcn in reactions:
            if rcn.status:
                # mute reactions if not organic reactant
                tmp = [i for i in rcn.reactants if species[sps.index(i)].organic]
                if tmp == []: 
                    print('RS: mute rcn cuz no organic reactants.',Type, i, rcn.reactants, rcn.products, rcn.rate.str)
                    rcn.status = False
                # change reactions if sps is in both reactants and products
                tmp1 = [i for i in tmp if i in rcn.products]
                if tmp1 != []:
                    if len(tmp1) > 1: 
                        raise ValueError('Found more than one sps repeat in reactants and products.',Type, tmp, tmp1, rcn.reactants, rcn.products, rcn.rate.str)
                    else: 
                        # get product index
                        i = rcn.products.index(tmp1[0])
                        # get new ratio
                        rt = 1. - rcn.ratiosPD[i]
                        #print('RS: Find sps in both reactants and products.',Type, i, tmp1[0], rcn.reactants, rcn.products, rcn.rate.str)
                        if rt <= 0. or rt > 1.: 
                            #print('-- mute rcn cuz rt not fits.',rt)
                            rcn.status = False
                        else:  
                            # new rate & ratios & species
                            for j in range(len(rcn.ratiosPD)):
                                rcn.ratiosPD[j] /= rt
                            rcn.rate.str += ('*{:6.4E}'.format(rt))
                            rcn.rate.update()
                            rcn.products.pop(i)
                            rcn.ratiosPD.pop(i)
                            #print('RS: Find sps in both reactants and products.',Type, i, tmp1[0])
                            #print('-- trim rcn with rt.',rt, rcn.reactants, rcn.products, rcn.rate.str)
                        
    return lumpsps, lumppds

def reduce_reaction_byRemove(reactions,species,Type,val,concs,frozen=[]):
    """remove all reactions related to the children species if the parent species are not included
        frozen are reaction()"""

    # prepare for the reduction of species
    sps1=get_species_from_reactions(reactions)
    sps=[i.name for i in species]

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
    sps=[i.name for i in species]
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

def species_type(X):
    """return the type of species"""
    
    if not X.organic: return 'inorg'
    s = X.name
    # basic type: radical, gas, aero
    if X.Radical: 
        itype = 'R' # ROO, ROOX, RO3, RO2, R
        if s[-2:] in ['OO','O3','O2']: return itype+s[-2:]
        elif s[-1] == 'O': return itype+'O'
        elif s[-3:-1] == 'OO': return itype+'OOX'
        elif  s[-1] == 'C': return itype+'C'
        else:
            if X.RO2: itype+'O2'
            else: itype
            #print('species_type: X radical type not found.', s)
    else:
        if X.condensed: itype = 'A' # aerosol
        else: itype = 'G' # gas 
        for i in CheckTypeKeys:
            if i in s: 
                return itype+i
    return itype

def ssh_reaction_type(rcn):
    """read the input reaction and return the type in SSH"""
    if '!!!' in rcn.rate.SSH: 
        print('find !!! in rcn.rate.SSH.')
        return None
    elif 'PHO' in rcn.rate.SSH:  return 'PHO'
    else:
        itype = [i for i in rcn.rate.SSH.split(' ') if (not isfloat(i) and  i not in ['KINETIC','']) ]
        if itype == []: raise NameError('No keywords (e.g. SPEC, ARR1, ARR2, ARR3, MCM, TB) is found.',rcn.rate.SSH)
        return itype
        
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
    

    # record species names
    Xname = X.name
    lumpsps=[Xname]
    # record the estimated parameter ratios in lumping
    # currently depending on the conc in the original case
    lumpratio=[rc[1][rc[0].index(Xname)]]

    # the generation number of X
    igen = 0
    # X type. if X is radical: itype in ['RO3','RO2','ROO','RO']
    # if X is non radical: 
    Xtype = species_type(X)

    # (2) search for species that can be lumped with X
    for s in Ys:

        # checking basic: not the same species X
        if s == Xname: continue

        # non-zero ref conc.
        #if rc[1][rc[0].index(s)] <= 0.0: continue

        # get s in species list
        i = sp[1][sp[0].index(s)]

        # check status
        if not i.status: continue

        # check PRAM species
        #if 'PRAM' in i.source:
        #    if 'PRAM' not in X.source: continue
        #elif 'PRAM' in X.source: continue

        # check species type
        if Xtype != species_type(i): continue
        
        # check Psat for aerosols
        if Xtype[0] == 'A' and CheckPsat:
            if abs(math.log10(X.psat_atm) - math.log10(i.psat_atm)) > CheckPsat: continue

        # for both radicals and non-radicals
        # check the minimum number of carbon
        #if min(X.functionalGroups['C'],i.functionalGroups['C']) <= 3 : continue 

        # check mass diff, min max
        #if abs(X.mass - i.mass) > 100. or min(X.mass,i.mass) <= 100. : continue
        if abs(X.mass-i.mass)/(X.mass+i.mass)*2 > 0.5: 
            #print('lp refuse: mass',Xname,s,X.mass,i.mass,abs(X.mass-i.mass)/(X.mass+i.mass)*2 )
            continue
            
        # check if has C=C structure: not for Radical
        if Xtype[0] =='A': #!= 'R':
            if 'C=C' in CheckStructureKeys:
                if not i.SOAPStructure or not X.SOAPStructure:
                    print(Xtype,Xname,s,i.SOAPStructure, X.SOAPStructure)
                    raise AttributeError('RS: SOAPStructure is not in both i and X.',i.name,X.name)
                tmp0 = 0 # for X
                tmp1 = 0 # for i
                for j in [16,17,18,19,20]: # C=C
                    if j in i.SOAPStructure[0]: tmp1 += i.SOAPStructure[1][i.SOAPStructure[0].index(j)]
                    if j in X.SOAPStructure[0]: tmp0 += X.SOAPStructure[1][X.SOAPStructure[0].index(j)]
                if tmp0 + tmp1 != 0 and tmp0 * tmp1 == 0: continue

        if 'C' in CheckStructureKeys:
            if abs(X.functionalGroups['C'] - i.functionalGroups['C']) > CheckStructure['C'] : 
                continue

        # check nitrate        
        if 'N/C' in CheckStructureKeys:
            if X.ratios['N/C'] * i.ratios['N/C'] == 0.0 and X.ratios['N/C'] + i.ratios['N/C'] != 0.0: continue

        # check oxygen  ! only for API    
        if 'O' in CheckStructureKeys:
            if  abs(X.functionalGroups['O'] - i.functionalGroups['O']) > CheckStructure['O'] : continue

        # generation
        if CheckGen and not isFitScale(gen[1][gen[0].index(i.name)],igen,CheckGen):
            continue

        # lifetime
        if CheckTau and not ifFitLifetime(X,i,tau,CheckTau): continue

        # structure and partitioning
        tag = 0
        for k in ['DU', 'formula', 'MW']:
            if k in CheckStructureKeys and ifFitStructure(X,i,k,CheckStructure[k]): tag = 1
        if tag: continue

        ### fits all criteria, now record i
        #print('lump: find fit: ',Xname,s, flush = True)
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

def trim_scheme(reactions,species): #,Type = None):
    """check the scheme and remove unused ones"""
    dlist = productions(reactions,'rea')
    sps0 = [i.name for i in species]

    sps = []
    for i in range(len(dlist[0])):
        n = 0
        while n < len(dlist[1][i]): #and not (set(primaryVOCs)&set(dlist[1][i])): # extend
            #if Type is None and (set(primaryVOCs)&set(dlist[1][i])): break
            if set(primaryVOCs)&set(dlist[1][i]): break
            s = dlist[1][i][n]
            if s in dlist[0]:
                for k in dlist[1][dlist[0].index(s)]:
                    if k not in dlist[1][i]: dlist[1][i].append(k)
            n += 1

        srcs = list(sorted(set(primaryVOCs) & set(dlist[1][i])))
        if srcs == []: # find species can not track back primary voc
            #logger.debug('RS trim_scheme: found isolated species, ',dlist[0][i])
            sps.append(dlist[0][i])
        #elif Type != None and srcs != species[sps0.index(dlist[0][i])]: # update source and group ID if need
        #    species[sps0.index(dlist[0][i])].source = srcs
        #    if srcs in groupIDs:
        #        species[sps0.index(dlist[0][i])].groupID = groupIDs.index(srcs) + 1
            #else:
            #    print('trim_scheme: i.source not in groupIDs!', srcs, dlist[0][i])

    # remove species
    for i in species:
        if i.name not in dlist[0]:
            if i.status and i.organic and i.name not in primaryVOCs:
                if i.name not in sps: sps.append(i.name)  
                if i.name in sps: i.status = False
    if sps != []: 
        #logger.debug('RS trim_scheme: ',sps)
        # remove reactions
        for r in range(len(reactions)):
            for i in reactions[r].reactants:
                if i in sps: reactions[r].status = False

    return reactions, species, sps

def reaction_merge(reactions,species):
    """Merge reactions with the same type in the reaction list"""
    
    reas = [[],[],[]] # store reactants, types, index
    inds = [] # index in reas[2] [a,b]
    # for each reaction, inds record the index in reas[0], reas[1], and reas[2]
    # sort
    for i in range(len(reactions)):
        rcn = reactions[i]
        if not rcn.status: continue
        tmp = sorted(rcn.reactants)
        if tmp in reas[0]: # find the same reactants
            s = reas[0].index(tmp) # get index in reas
            rate = ssh_reaction_type(rcn)
            if rate in reas[1][s]: # same type
                p = reas[1][s].index(rate)
                inds.append([s,p,len(reas[2][s][p])])
                reas[2][s][p].append(i)
            else:
                reas[1][s].append(rate)
                inds.append([s,len(reas[2][s]),0])
                reas[2][s].append([i])
        else: # not find, build new
            reas[0].append(tmp)
            reas[1].append([ssh_reaction_type(rcn)])
            inds.append([len(reas[2]),0,0])
            reas[2].append([[i]])
        #print(i, len(reas[0]), reas[0][-1], reas[1][-1],reas[2][-1],rcn.reactants)

    # merge
    new_rcns,nomers,merged = [],[],[]
    for i,j in enumerate(inds):
        s,p,q = j 
        # s is species index: reas[0][s]
        # p is rate type index: reas[1][s][p]
        # q is rcn index in reactions reas[2][s][p][q]
        # ircn is the index of targeted reaction in reaction list
        ircn = reas[2][s][p][q]
        
        # check reaction index
        if ircn in merged: continue
        
        # check number of similar reactions
        k = len(reas[2][s][p])
        
        if k == 1 or ircn in nomers: # one reaction, no merge
            new_rcns.append(reactions[ircn])
            
        elif k >= 2 : # need to merge
            tmp = [reactions[l] for l in reas[2][s][p]] # merged reactions
            new_rcn = merge_reactions(tmp, reas[1][s][p]) # reactions and the type
            
            # check merged results
            if new_rcn: 
                new_rcns.append(new_rcn)
                # pervent from merge the same reactions
                for l in reas[2][s][p]:
                    if l not in merged: merged.append(l)
            else: # can not merge
                new_rcns.append(reactions[ircn])
                # update nomers
                for l in reas[2][s][p]:
                    if l not in nomers: nomers.append(l)
        else:
            raise ValueError('len(reas[2][s]) <= 1.',k,i,j,reas[0][s],reas[1][s],reas[2][s])

    return  new_rcns
    
def reaction_seperate(reactions,species,Type = 'nokdec'):
    """rewrite reactions to those with seperqte products"""
    # init
    new_rcns = []

    # get species list
    sps=[i.name for i in species]

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
        elif 'SPEC' in i.rate.SSH and 'TB' in i.rate.SSH: # not change if SPEC with TB
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
    """merge input reactions if they have the same reactants"""

    if len(rcns) < 2 or not Type or Type == 'PHO': # check number of input reactions
        #print('merge_reactions input number of reactions is less than 2! or type is not recognized.',len(rcns),Type)
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
