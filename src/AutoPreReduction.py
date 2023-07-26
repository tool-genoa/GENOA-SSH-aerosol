# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  Python script to run prereduction depending on the given reduction
#   parameters and options. 
#
# ================================================================================

from Parameters import RunSets, primaryVOCs
from ReductionStrategy import reduction_chemistry_scheme, trim_scheme, \
                              reaction_seperate, reaction_merge
from KineticMCMtoPython import update_kinetic_rate

# order to conduct prereduction
RunSetsOrder = ['Psat_aero','Psat_NVOC','Psat_SVOC','gen','conc','kdec','tau','lump','bratio']
        
def auto_prereduction(reactions,species, setups = RunSets[0]):

    fout=[] # init for output reductions
    
    print('\nActive pre-reduction ...')
    n0, s0 = len([i for i in reactions if i.status]), len([s for s in species if s.status])
    # print numbers before reduction to check
    print('Number of reactions and species before prereduction: ', n0, s0)

    # Remove the gas-particle partitioning of primaryVOCs
    # Those won't be removed in currrent reduction
    for s in species:
        if s.name in primaryVOCs and s.condensed:
            print('Remove the condensability of primary VOC: ', i.name)
            s.condensed = False

    # Update branching ratios if prereduction with lumping or removing reactions
    if setups['lump'] or setups['bratio']:
        concs_ave, concs_min, refconc_paths, RefConcRead = update_kinetic_rate(locs)
    else: 
        concs_ave, concs_min, refconc_paths, RefConcRead = [],[],[], None
    
    # Check options by the order in RunSetsOrder
    for key in RunSetsOrder:

        if setups[key]:

            print('\nPrepare prereduction by the keyword: ',key,' value: ', setups[key])
            
            rdcs = reduction_chemistry_scheme(reactions,species,key,
                                              setups[key], 
                                              concs_ave, concs_min,
                                              refconc_paths, RefConcRead)

            print('criteria: ',key,', reduced species: ',rdcs, ', No. ',len(rdcs))

            # write output in detail
            fout.append('Reduction by keyword {:s} with value {:s}, number: {:d}\n'.format(key,str(setups[key]),len(rdcs)))
            
            # record reductions
            if key == 'lump': 
                for i in rdcs: 
                    for j in i: fout.append(j+'\t')
                    fout.append(i+'\n')
            else:
                for i in rdcs: fout.append(i+'\t')
                fout.append('\n\n')

    # print numbers before reduction to check
    n1, s1 = len([i for i in reactions if i.status]), len([s for s in species if s.status])
    if n0 != n1 or s0 != s1:
        print('Number of reactions and species after prereduction phase I: ', n1, s1)

    # Phase II prereduction
    # Trim mechanism: reactions & species that can not be formed from primary VOCs
    print('\nTrimming useless reactions & species that can not be formed from primary VOCs: ', primaryVOCs)
    reactions,species,trim_sp = trim_scheme(reactions,species)
    n2, s2 = len([i for i in reactions if i.status]), len([s for s in species if s.status])
    if trim_sp != [] or n1 != n2 or s1 != s2: 
        tmp = 'Current No.rcn: {:d} & No.sps: {:d}.\nRemove useless species No.{:d}: '.format(n2,s2,len(trim_sp))
        
        # ouput to check
        print(tmp, trim_sp)
        fout.append(tmp)
        for s in trim_sp: fout[-1] += ('\t'+s)
        fout[-1] += '\n\n'
    
    # Combine reactions
    print('\nMerging reactions ...')
    reactions = reaction_merge(reactions,species)
    n3, s3 = len([i for i in reactions if i.status]), len([s for s in species if s.status])
    # check numbers of reactions
    if n2 != n3: 
        print('Number of reactions is changed from {:d} to {:d}.\n'.format(n2,n3))

    # save output report
    fout.append('Before prereduction: No.rcn: {:d} \t No.sps: {:d}\n'.format(n0,s0))
    fout.append('After prereduction phase I: No.rcn: {:d} \t No.sps: {:d}\n'.format(n1,s1))
    fout.append('After prereduction phase II: No.rcn: {:d} \t No.sps: {:d}\n'.format(n3,s3))
        
    return reactions, species, fout
