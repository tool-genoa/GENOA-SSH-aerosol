# -*- coding: utf-8 -*-
#================================================================================
#
#     GENOA v1.0: the GENerator of reduced Organic Aerosol mechanism
#
#     Copyright (C) 2022 CEREA (ENPC) - INERIS.
#     GENOA is distributed under GPL v3.
#
#
# This is a file to clean up the output info from GENOA
#================================================================================

import os
import shutil

import Parameters

#### NEED TO CHECK HERE !!!
# put 1 for the folder need to be cleaned
options = {'SSHs':1,'smls':1,'chems':1}
# all the generated chemical mechanism with ID mentioned in the cfg file will be cleaned
# except for the ones kept in the list and the final reduced mechanism (if testing file exists):
keep_chems = ['BCARYorg', 'BCARYRa']
####
# execute command:
# cd ../GENOA/src # go to the src code
# ptyhon3 ../examples/clean.py [cfg_name]
# [cfg_name] example: ../inputs_bcary/rdc_cfg.ini
########

## clean procedure
# clean ssh-aerosol
if options['SSHs']:
    print('Cleaning SSHs... ')
    for i in pathSSH_rdc, pathSSH_ref:
        shutil.rmtree(i, ignore_errors=True)
        print('Folder ',i,' is removed.')

if options['smls'] or options['chems']:
    # get paths
    path_sav_rec = '{0:s}/{1:s}_recs'.format(pathNewChem,tail)
    path_sav_chem = '{0:s}/{1:s}_chems'.format(pathNewChem,tail)
    path_sav_res = '{0:s}/{1:s}'.format(pathNewRes,tail)

    # get the final case name from Testing file
    chem_files = os.listdir(path_sav_chem)
    for i in chem_files:
        if 'Testing_' in i: # get name of the final reduced mechanism
            ichem = i.split('_')
            if len(ichem) == 3 and prefix in ichem[1]: # found
                if ichem[1] not in keep_chems: # add into keep_chems if not in
                    keep_chems.append(ichem[1])

    if options['smls']:
        print('\nCleaning smls... ')
        files = os.listdir(path_sav_res)
        for i in files:
            if 'Results_' in i:
                if i.replace('Results_','') not in keep_chems:
                    shutil.rmtree(path_sav_res+'/'+i, ignore_errors=True)
                else:
                    print(i,' is kepted.')
        print('Folder ', path_sav_res,' is cleaned.')

    if options['chems']:
        print('cleaning chems... ')
        shutil.rmtree(path_sav_rec, ignore_errors=True)
        print('Folder ',path_sav_rec,' is removed.')
        for i in chem_files:
            tag = 0
            if 's_' in i: tag = 1
            elif prefix in i and i[0:len(prefix)] == prefix:
                if i not in keep_chems: tag = 1
            if tag: shutil.rmtree(path_sav_chem+'/'+i, ignore_errors=True)
            else: print(i,' is kepted.')
        print('Folder ', path_sav_chem,' is cleaned.')

    # print info 
    print('Keep the chemical files and SOA simulation results related to the reduced mechanisms: ', keep_chems)
