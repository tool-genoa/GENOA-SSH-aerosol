# -*- coding: utf-8 -*-
#================================================================================
#
#     GENOA v1.0: the GENerator of reduced Organic Aerosol mechanism
#
#     Copyright (C) 2022 CEREA (ENPC) - INERIS.
#     GENOA is distributed under GPL v3.
#
#================================================================================

# this python script contains functions that processes a list of reactions and returns species properties related to the reaction list

import numpy as np

from Parameters import NokeepSp,primaryVOCs
from Functions import compare,isfloat,isContain,array_to_ratio

def get_species_from_reactions(reactions):
    """return a list of species name that appear in the reactions"""
    # species list
    sps=[]
    for i in reactions:
        if not i.status: continue
        for j in i.reactants,i.products:
            for k in j:
                if k == '': continue
                elif k not in sps: sps.append(k)

    return sps

def destructions(reactions,Type=None):
    """return a list of info related to the destruction of species"""
    sps=[[],[]] # name, [[reactants],[products],rate]
    sps0 = [] # all species that occurs in the products 
    for i in reactions:
        if not i.status: continue
        # check species name
        #k = reactions.index(i)
        for j in i.products:
            if j not in NokeepSp and j not in sps0: sps0.append(j)
        for j in i.reactants:
            if j not in NokeepSp:
                # add in the list if not exists
                if j not in sps[0]: 
                    sps[0].append(j)
                    sps[1].append([])
                # record produced species and ratio, kinetic rate
                if Type is None:
                    sps[1][sps[0].index(j)].append([list(zip(i.ratiosRC,i.reactants)),list(zip(i.ratiosPD,i.products)),i.rate.pyformat])
                elif Type =='reactant':sps[1][sps[0].index(j)].append(list(zip(i.ratiosRC,i.reactants)))
                elif Type == 'rea': sps[1][sps[0].index(j)].append(i.reactants)
                elif Type =='product': sps[1][sps[0].index(j)].append(list(zip(i.ratiosPD,i.products)))
                elif Type =='pro': 
                    for k in i.products: 
                        if k not in NokeepSp and k not in sps[1][sps[0].index(j)]: sps[1][sps[0].index(j)].append(k)
                elif Type =='rate':sps[1][sps[0].index(j)].append(i.rate.pyformat)
                else:
                    raise TypeError('CR destruction: given type unknown: '+Type)
    for j in sps0:
        if j not in sps[0]: # not in the destruction list
            sps[0].append(j) # add the species that not in the destruction list
            sps[1].append([])
    return sps

def productions(reactions,Type=None):
    """return a list of info related to the production of species"""
    sps=[[],[]] # name, [[reactants],[products],rate]
    for i in reactions:
        if not i.status: continue
        # check species name
        for j in i.products:
            if j not in NokeepSp:
                # add in the list if not exists
                if j not in sps[0]: 
                    sps[0].append(j)
                    sps[1].append([]) 
                # record produced species and ratio, kinetic rate
                if Type is None:
                    sps[1][sps[0].index(j)].append([list(zip(i.ratiosRC,i.reactants)),list(zip(i.ratiosPD,i.products)),i.rate.pyformat])
                elif Type =='reactant':sps[1][sps[0].index(j)].append(list(zip(i.ratiosRC,i.reactants)))
                elif Type =='product':sps[1][sps[0].index(j)].append(list(zip(i.ratiosPD,i.products)))
                elif Type =='rate':sps[1][sps[0].index(j)].append(i.rate.pyformat)
                elif Type == 'rea':
                    for k in i.reactants: 
                        if k not in NokeepSp and k not in sps[1][sps[0].index(j)]: sps[1][sps[0].index(j)].append(k)
                elif Type == 'pro':
                    for k in i.products: 
                        if k not in NokeepSp and k not in sps[1][sps[0].index(j)]: sps[1][sps[0].index(j)].append(k)
    return sps

def reaction_chains(reactions, species):
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
    return plist,dlist

def branch_ratio(reactions, concs):
    "compute 1/k[side reactants]"
    reactants = [] # record all reactants
    ratios=[]
    reaction_inds=[]
    other_reactants = []

    for i in range(len(reactions)):
        if not reactions[i].status: continue

        rcn = reactions[i]

        rcn.check_type() # update species
        rcn.rate.update_value() # update kinetic rate

        # check number of species
        for sps in rcn.species[0]:
            #sps = rcn.species[0][0]

            # branch ratios
            tauD = rcn.rate.pyformat

            # species list
            tmp = [] # side reactants
            if len(rcn.reactants) > 1:
                for j in rcn.reactants:
                    if j not in rcn.species[0]: tmp.append(j)

            if sps in reactants:
                # update
                ind = reactants.index(sps)
                # add ratios,reaction_inds,other_reactants
                other_reactants[ind].append(tmp)
                # add and reaction index
                ratios[ind].append(tauD)
                reaction_inds[ind].append(i)
            else:
                # build new
                reactants.append(sps) 
                ratios.append([tauD])
                reaction_inds.append([i])
                other_reactants.append([tmp])

    # get ratios
    nc = len(concs[list(concs.keys())[0]]) # number of case
    nt = len(concs[list(concs.keys())[0]][0]) # number of time
    nr = 0
    for i in range(len(ratios)):
        tag = 0
        # check i need concs
        for j in other_reactants[i]:
            if j != []:
                tag = 1
                break

        if tag:
            rt = []
            for t in range(nt):
                for c in range(nc):
                    tmp = [] # current ratio
                    for j in range(len(other_reactants[i])):
                        tmp.append(ratios[i][j])
                        for k in other_reactants[i][j]: 
                            if k in concs.keys(): tmp[-1] *= concs[k][c][t] # concs unit
                            else: tmp[-1] = 0.0 # concs of this side reactant is 0.0
                    tmp = array_to_ratio(tmp)
                    if rt == []: rt = tmp
                    else:
                        for k in range(len(tmp)):
                            if tmp[k] > rt[k]: rt[k] = tmp[k] # only keep the maximum in all case and all time
            ratios[i] = rt 
        else: ratios[i] = array_to_ratio(ratios[i]) # no side reactant

    return reactants,[ratios,reaction_inds,other_reactants]

def tree(reactions):
    """ input a list of reactions
        extract all the species (except the common species) in a tree format
        reactants -> products as parents -> children species
        return lists of parent and children species"""

    parents=[]
    children=[]

    # traverse all reactions
    for i in reactions:
        # check status
        if not i.status: continue
        for j in i.reactants:
            if j in NokeepSp : continue
            # add new children and/or parent
            if j not in parents:
                parents.append(j)
                children.append([])
                ind = -1
            else: ind = parents.index(j)
            for k in i.products: 
                if k not in children[ind]: children[ind].append(k)
            
    # traverse all species
    tag=1
    while tag:
        tag=0
        for l in range(len(children)):
            for j in children[l]:
                # if children in parent list, check if their children is recorded
                if j in parents:
                    ind=parents.index(j)
                    for s in children[ind]:
                        if s not in children[l]: 
                            children[l].append(s)
                            tag=1

    return parents,children

def reversed_tree(reactions):
    """ input a list of reactions
        extract all the species (except the common species) in a reversed tree format
        children species -> all reactants before
        return lists of children and parent species"""

    parents=[] + primaryVOCs
    children=[[]]
    # traverse all reactions
    for i in reactions:
        # check status
        if not i.status: continue
        for j in i.products:
            if j in NokeepSp : continue
            # add new parents
            if j not in parents:
                parents.append(j)
                children.append([])
                ind = -1
            else: ind = parents.index(j)
            # add i.reactants
            for k in i.reactants:
                if k not in NokeepSp and k not in children[ind]: children[ind].append(k)

    # traverse all species
    tag=1
    while tag:
        tag=0
        for l in range(len(children)):
            for j in children[l]:
                if j in parents:
                    ind=parents.index(j)
                    for s in children[ind]:
                        if s not in children[l]: 
                            children[l].append(s)
                            tag=1

    return parents,children

def lifetime(reactions, concs):
    """return species name and their lifetime depending on the given reaction list
        lifetime: tau=1/(sum(k[Conc. of other reactants]))"""

    # update reaction
    for i in reactions: i.rate.update_value()

    # get destruction info
    destructionList = destructions(reactions) # name, [[reactants],[products],rate]
    #productionList = productions(reactions)

    # set species name
    speciesName=list(set(destructionList[0]))

    #speciesName=set(productionList[0]+destructionList[0])

    # storage
    lftime=[]

    #for i in concs.keys(): concs[i] = np.average(concs[i])

    for i in speciesName:

        # destruction
        tauD=0. # for one species
        if i in destructionList[0]:
            ind = destructionList[0].index(i)
            for j in destructionList[1][ind]:
                subtau = j[2] # rate
                for k in j[0]:
                    # reactant == species
                    if k[1] == i: tmp=1./k[0]
                    # reactant != species, given in concs.
                    else: 
                        if k[1] in concs.keys():
                            tmp=k[0]*concs[k[1]]
                        else:
                            tmp=k[0]*concs['RO2'] # RO2 + RO2 reactions dimer
                            #raise NameError('lifetime: NokeepSp species conc. not given, '+k[1])
                    subtau *= tmp # [conc]*k
                tauD+=subtau # sum([conc]*k)
        if tauD == 0.: lftime.append(1E9) # not destruction
        else: lftime.append(1./tauD)

    #print('lifetime is computed for {:d} species.'.format(len(speciesName)))
    return speciesName,lftime

def generation(reactions, Type = None):
    """return species and their number of generation depending on the given reaction list"""

    pros = destructions(reactions,Type='pro')
    reas = destructions(reactions,Type='rea')
    gens = {}
    for i in primaryVOCs: gens[i] = 1
    used = []

    while len(list(gens.keys())) < len(pros[0]) or len(used) == len(pros[0]):
        for i in list(gens.keys()):
            if i not in pros[0]:
                if i in primaryVOCs:
                    print(i, 'primary VOC is not in pros list')
                    raise AttributeError('CR generation: Check reaction list.')
                else: print(i) 
                continue
            n = pros[0].index(i) # index
            if n in used: continue
            tag = 0
            for j in reas[1][n]: # check oxidants in the atmosphere
                for k in j:
                    if k in ['OH','O3','NO3']:
                        tag = 1
                        break
                if tag: break

            for j in pros[1][n]: # products
                if j not in list(gens.keys()): gens[j] = []
                for k in gens[i]:
                    if i == j:
                        print(i,j,k)
                        #raise AttributeError('CR generation: species form itself.')
                        continue
                    if tag: gens[j].append(k+1) # gen +1
                    else: gens[j].append(k) # gen no change
            used.append(n)
    
    for i in list(gens.keys()):
        gens[i] = list(set(gens[i]))
        gens[i].sort()

    if Type == 'all': return gens

    speciesName = get_species_from_reactions(reactions)
    sps_generation = []
    nmax = 0

    for i in speciesName: 
        if i in list(gens.keys()):
            sps_generation.append(min(gens[i]))
            nmax = max(nmax, max(gens[i]))
        else:
            sps_generation.append(0)
    return speciesName,sps_generation
