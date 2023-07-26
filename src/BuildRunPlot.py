# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  BuildRunPlot.py converts the knietic rate contants from MCM format
#
#   to a format can be executed by SSH-aerosol. 
#
#  The converted format is used to output SSH-aerosol aerosol species list.
#
# ================================================================================

import os

def build_run_plot():

    from Parameters import reactionfile,speciesfile, \
                           speciestype,reactiontype, \
                           locs,RunSets, \
                           pathNewChem, IDchem, prefix, \
                           pathNewRes, pathSSH,SSH_mode, \
                           Ttot,DeltaT, \
                           Tnow,initfile,pathInitFiles, tail, \
                           primaryVOCs, ResultFolder

    # read mechanisms in other formats & prereduction
    if RunSets[0]['NewChem']: # test tag for part A

        from DataStream import to_SSH_sets,read_chem_sets
        from AutoPreReduction import auto_prereduction

        print('Generate New SOA Mechanism Files...\n')

        # get lists for reactions and species
        rc,sp=read_chem_sets(reactionfile,speciesfile,speciestype,reactiontype)

        # if reduction options are avtived
        if RunSets[0]['prereduction']:
            rc,sp,fout = auto_prereduction(rc,sp,RunSets[0])
        else: fout = []

        # output as SSH files
        to_SSH_sets(pathNewChem,IDchem,rc,sp,'20')
        if fout != []:
            with open (pathNewChem+'/'+IDchem+'/Reduction.sav','w+') as f:
                for i in fout: f.write(i)
                
    # embed chemistry in SSH and init/run simulations
    if max(RunSets[1].values()):
        from SSHSetUp import set_SSH,get_namelist,run_SSHsimulation
        from Functions import move_results

        # check if chem file exist
        if IDchem not in os.listdir(pathNewChem):
            print('Chem file not exists: ',IDchem,pathNewChem, os.listdir(pathNewChem))
            raise NameError('BRP: not found IDchem '+IDchem)

        if RunSets[1]['compile']:
            # clean and compile SSH-aerosol with the new chem
            print('Compile new SSH... in ',pathSSH)
            set_SSH(IDchem,pathNewChem,pathSSH,mode = SSH_mode)

        if RunSets[1]['simulate']:
            # write the namelist with the repository of the new chem
            print('Write new SSH namelist... ')

            Namelist = '{:s}/namelist_init.ssh'.format(pathSSH)
            ssh_namelist, ssh_inds = get_namelist()
            
            # Run simulations with settings provided in Parameters.py
            print('Run simulation... results in ',ResultFolder)
            results,errs  = run_SSHsimulation(tail,ResultFolder,pathSSH,
                                              ssh_namelist, ssh_inds,
                                              Ttot,DeltaT,[[],locs],
                                              Tnow,initfile = initfile,
                                              pathInit = pathInitFiles)
        else: results = []

        if 'clean' in RunSets[1].keys() and RunSets[1]['clean']:
            # save results: clean and mv
            path_old = '{:s}/{:s}'.format(pathSSH,ResultFolder)
            path_new = '{:s}/{:s}'.format(pathNewRes,ResultFolder)
            print('Move SOA simulation results from {:s} to {:s}...'.format(path_old,path_new))
            if results != []: # get secondary path, lc
                items = [i.split('/')[1] for i in results]
            else: items = [] # move everything in the folder
            move_results(path_old, path_new, items)

    # analysis simulation results and reduce the chemical scheme again if need 
    if RunSets[2]['display']:
        from SSHResultProcess import get_organic_conc, display_total_conc

        pathSSH = pathNewRes

        # get tails for ref and current cases
        if '_tmp' in ResultFolder: 
            case_tail = ResultFolder.split('_')[-2].replace(prefix,'')+'_tmp'
        else:
            case_tail = ResultFolder.split('_')[-1].replace(prefix,'')

        ref_tail = RunSets[2]['refPath'].split('_')[-1].replace(prefix,'')

        # prepare for ploting
        if RunSets[2]['cmp']: 
            plot_paths = RunSets[2]['cmp']
            savpathplot = '{:s}'.format(RunSets[2]['savpath'])
        else: 
            plot_paths = [RunSets[2]['refPath']] + [pathSSH+'/'+ResultFolder]
            savpathplot = '{:s}/{:s}'.format(RunSets[2]['savpath'],ResultFolder)

        np = len(plot_paths)

        # tail
        # plot_tails = [ref_tail,case_tail]
        plot_tails = []
        for i in plot_paths: plot_tails.append(i.split('_')[-1].replace(prefix,''))

        # labels
        if RunSets[2]['cmpLabel']: plot_labels = RunSets[2]['cmpLabel']
        else: plot_labels = plot_tails

        # styles
        if RunSets[2]['cmpstyle']: plot_styles = RunSets[2]['cmpstyle']
        else: plot_styles = ['-','--']

        # color
        if RunSets[2]['cmpcolor']: plot_colors = RunSets[2]['cmpcolor']
        else: plot_colors = ['b','r','g','k','c','orange']

        print('tails: ',plot_tails,', labels: ',plot_labels, 'styles: ', plot_styles, 'colors: ', plot_colors)

        os.makedirs(savpathplot, exist_ok=True)
        fdtl = open ('{:s}/details'.format(savpathplot),'w+')

        fdtl.write('path\tave SOA\tmax SOA\tfinal SOA\n')
        # plotting
        for inow in Tnow:
            for idx in range(len(Ttot)): # for different time settings
                # time settings
                iT=Ttot[idx] #total time
                iD=DeltaT[idx] #time step

                # locations
                for il in locs:
                    y = il[0]
                    x = il[1]

                    if initfile in ['month','storage']: lc = 'm{:d}y{:d}x{:d}'.format(il[2],il[0],il[1])
                    else: lc = 'y{:d}x{:d}'.format(il[0],il[1])

                    paths = []
                    for p in range(np):
                        # get path for folder
                        tmp ='{:s}/{:s}/c_{:d}_{:d}_{:d}h_{:s}'.format(plot_paths[p],lc,iT,iD,inow,plot_tails[p])

                        # check if toto file exist
                        tag = 0
                        for p1 in os.listdir(tmp):
                            if 'toto' in p1: tag = 1

                        # add path
                        if tag:
                            paths.append(tmp)
                        else:
                            raise NameError('BRP: no toto file found in the path: '+tmp)

                    displayPath = '{:s}/{:s}/SOAs_{:d}h'.format(savpathplot,lc,inow)

                    # get total organic conc
                    data = get_organic_conc(paths)

                    # plot total organics
                    rmse = display_total_conc(data,plot_labels,savpath=displayPath,
                               Type='plain',lst=plot_styles,col = plot_colors,
                               errType = RunSets[2]['err_type'])

                    print('path: ', '{:s}/{:s}/c_{:d}_{:d}_{:d}h_{:s}'.format(plot_paths[1],lc,iT,iD,inow,plot_tails[1]))

                    for k in range(len(paths)-1):
                        tmp = [paths[k+1]]
                        tmp.append('{:6.4f}'.format(sum(data[1][k+1])/len(data[1][k+1]))) # average SOA
                        tmp.append('{:6.4f}'.format(max(data[1][k+1]))) # max SOA
                        tmp.append('{:6.4f}'.format(min(data[1][k+1][1:]))) # min SOA
                        tmp.append('{:6.4f}'.format(data[1][k+1][-1])) # final SOA
                        print('path, average, max, min, final SOAs:',tmp)

                        fdtl.write(tmp[0])
                        for sub in tmp[1:]: fdtl.write('\t'+sub)
                        fdtl.write('\n')
        fdtl.close()

    # check up information of the new chemistry
    if not RunSets[0]['NewChem']: 
        from DataStream import read_chem_sets
        rc,sp=read_chem_sets(pathNewChem+'/'+IDchem+'/'+IDchem+'.reactions',speciesfile,'SSH','SSH')

    # original print out messages # complete mode
    # get number of reactions/species
    from Functions import get_info
    tmp = pathNewChem+'/'+IDchem+'/'+IDchem
    rea, gas, aer, nro2, nph = get_info(rc,sp,'com')

    print('Number of gas species: ', gas)
    print('Number of aerosol species: ', aer)
    print('Number of reaction: ', rea)
    print('Number of RO2: ', nro2)
    print('Number of Photolysis: ', nph) 
