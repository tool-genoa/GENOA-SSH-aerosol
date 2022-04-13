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
import sys
import json
import numpy as np
import configparser

# output chem name VOC+[tail] and simulation path in ssh.chem.mcm/
IDchem = None

# mode for ssh-aerosol
## 'fast': only output namelist, total SOA conc.
## 'complete': output everything
# mode for generating ref. concs. from fake chem case
# takes long time
## 'compile': compile ssh-aerosol
SSH_mode = 'fast_compile'

# basic namelist, default in the path: src/ssh-genoa/
namelist_pre = 'namelist_ssh'
initfile = 'storage'
data_type = 'aero' # type of analysed data
Tnow = [0,12] # under the condition at which hour

# aerosol output settings
RunSets=[
         # 0: new chem sets
         {'NewChem':0},
         # 1: simulation
         {'compile':0,'simulate':0,'clean':0},
         # 2: postprocessing
         {'refPath': None,# reference case
          'display':0,
          'err_type':'fe-2',
           # 0: default, campare with the ref case
          'cmp':[],
          'cmpcolor':[],
          'cmpLabel': ['MCM','Rdc.'],
          'cmpstyle': ['-','--',':'],
          'items': [],#['ratio','fun_8,72','Kp_8,72'],
          'savpath':'../results/graphs'
          },

         # 3: Setups for Auto-Processing.py
         {
          'from_ref': False,
          # pre case
          'IDchemPre': None, # pre chem
          'preChemPath': None,
          'prePath': None,
          # ref case
          'IDchemRef': None, # ref chem
          'refChemPath': None,
          'refPath': None,
          # fake case
          'IDchemFake': None, # fake
          'fakePath': None,

          # options
          'tag_quick_mode': False,
          'tag_pre_reduction': False,
          'tag_redo': False, # redo if error exceed pre-testing results
          'tag_hybrid': False, # use parallel in series reduction

          # reduction strategies
          'auto_types':['rm','rp','lp','mg','rs','da'], # rm: remove reactions; rs: remove not used species, da: remove aerosol species, rr:remove radicals (not used now); lp: autolumping
          'BranchRatio': [5E-2, 1E-1, 5E-1], # only for AutoRemoving, ratio_remove times, but need to set for each auto_types
          # error analysis
          'err_ref':[0.01,0.02,0.02,0.03,0.03,0.03,0.04,0.04,0.06,0.06,0.08,0.08,0.10,0.10], # error tolerance compared to ref case
          'err_pre':[0.01,0.01,0.02,0.01,0.02,0.03,0.02,0.04,0.04,0.06,0.04,0.08,0.08,0.10], # error tolerance compared to pre case
          'tag_check_RO2': False,
          'tag_check_average': True,
          'tag_check_median': True,

          # with testing after each reduction or not
          'tag_testing':True,

          # setting for restriction on pre-testing
          'try_at_err': 99.,  # where to start if err_ref >= try_at_err
          'try_ave_ref': 99., # limitation on the err_ave_max for pre-testing
          'try_max_ref': 99., # limitation on the err_max for pre-testing

          # number of pre-testing
          'nPreTest': 150,
          'naer_late-stage': 20 # number of aerosol that indicates late stage of reduction
         },

        # 4: Setups for Auto_Testing.py
        {
         # reduction set
         'Test_file': None,
         # reference
         'IDchemRef':None, # ref chem
         'refChemPath':None,
         'orgpath': 'ref', # if save organics_1.txt for multiple testing
         'Trytimes': False, # test locs number e.g., 100 # not used if ind is provided.
         # error
         'err_max': 1., # error tolerance. Program exits if finds error >= err_max * 3.
         'err_out': 1.
        }]

# order of reduction
pathNewChem = '../toSSH/'
pathNewRes = '../results/'

primaryVOCs = [None]
prefix = None

# Reduction Options
Roptions = {
            'OrderType' : 'reaction', # 'conc' 'reaction'
            'PeerType' : 'svoc', # 'reaction' 'all' 'svoc'
            'CheckLink' : True,  # 'True' 'False'
            'RO2Treat' : None,   # 'fixed'
            'CheckTau' : 10,     # value of tolerance, float
            'CheckGen' : False,  # generation
            'CheckPsat' : False, # Psat_atm for SVOC
            'CheckStructure' : {'PAN': True,'N/C': True,
                                'OOH': True,'CO3H': True,
                                'C=C': True, # check if exists
                                #'O': 3, # diff O <= 3
                                'C': 2 # diff C <= 2
                                },
            'FreezeSpecies':[] # frozen species
            }

# aerosol output criteria # only used when building single reduction case
AeroDict = {
            #'min','max','ave','avelog10','evap','simpol', 'ave_log10_3'
            'vpType':'VP1BP2',
            'Psat_aer':0.0,
            # reduction options
            'Psat_NVOC': 1E-13,
            'Psat_SVOC': 1E-4,
            'conc':{'aero':1E-5,'gas':0.0,'TM':0.0},
            # lifetime
            'tau':0,
            'gen':10,
            'kdec':True,
            'Bratio':'1E-3_1',
            'lump':'DU_1'}

# output species besides RO2, Organics
outputs = {'aero': [], # 'PH2O'
           'gas': []} # primaryVOCs[0]

# repositories of MCM files and SSH model
speciestype = 'SSH'
speciesfile = None
reactiontype = 'SSH'
reactionfile= None

pathSSH = '../SSHs/'
pathInitFiles = None

# location # [y,x,m] [lat,lon,month] # lat: 32-70 lon: -17-39.8 month 0-11
# index of locations in ncfile
locations = []
locs = []

print('Running python file: ',sys.argv[0])

# read config file if provided
if len(sys.argv) > 1:
    cfgfile = sys.argv[1]
    print('Reading parameters from config file: ', cfgfile)

    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(cfgfile)

    # assign parameters directly in certain sections
    for s in ['chemistry_id','conditions','input','output']:
        for item in cfg.items(s):
            i, j = item # get name and value
            try:
                tmp = eval(i) # check if i exist
                exec('{:s} = {:s}'.format(i, j)) # assign
            except NameError:
                print('not able find parameter: ',i)
                raise
            tmp1 = eval(i)

    # assign conditional parameters
    inds_cfg = { # index in RunSets
                'training' : 3, 'testing' : 4,
                'simulation': 1, 'postprocessing': 2}

    for s in cfg['action']:
        if s not in list(inds_cfg.keys()):
            print('Not recognize action: ',s)
            continue
        elif cfg['action'][s] in ['True', '1']:
            k = inds_cfg[s]
            if s == 'postprocessing': RunSets[k]['display'] = 1
            for item in cfg.items(s):
                i, j = item # get name and value
                if i in list(RunSets[k].keys()):
                    tmp = RunSets[k][i]
                    try:
                        exec('RunSets[k][i] = {:s}'.format(j)) # assign
                    except:
                        RunSets[k][i] = j # assign
                    # print to check update
                else:
                    print('not able find parameter: ',i, 'in RunSets[{:d}]'.format(k), cfgfile)
                    raise NameError('cfg: RunSets')

if locs == []:
    print('The input locs is empty.', cfgfile)
    raise ValueError('Check locs in the config file')
elif isinstance(locs,str):
    print('Read locs from the entire list', locs)
    with open (locs,'r') as f: locs = json.loads(f.read())
    print('read: ',len(locs),' locs')
elif  isinstance(locs,list):
    if isinstance(locs[0],str):
        n = int(locs[1])
        print('Read the first ',n,' locs from the list', locs[0])
        with open (locs[0],'r') as f: locs = json.loads(f.read())[:n]
    else:
        print('Read locs list: ',locs)
else:
    raise TypeError("locs type not recognized.", locs)
##
# settings for running SSH simulations
tail = IDchem.replace(prefix,'')

if pathSSH[-1] == '/': pathSSH = pathSSH[:-1]

pathSSH_rdc = '{:s}_rdc'.format(pathSSH)
pathSSH_ref = '{:s}_ref'.format(pathSSH)
pathInitFiles = os.path.abspath(pathInitFiles)

# other default settings
#TINYM = 1E-20
TINYM=1E-99

# simulation time
Ttot = [432000]
DeltaT = [3600]

### species list and relative properties
# used in datastream: SSHSpeciesInit and SSHSpeciesInitMW: species name & their MWs
# NOT in reactions file but need to put in SSH species lists
SSHSpeciesInit={'SULF':98.,'NH3':17.,'HNO3':63.,'HCL':36.5,
                'RO2':250.,'SOAlP':392.,'airm':29.,'BiMT':136.,
                'OH':17.,'NO3':62.,'O3':48.,'CO':28.,
                'NO2':46.,'NO':30.,'HO2':33.,'H2O2':34.,
                'H2':2., 'SO2':64., 'SO3':80.}

# species soutput as constants in SSH-aerosol
cst=['OH','NO3','O3','NO','NO2','HO2','CO']

# species kept as constants in the simulation
# conc. # unit: ug/m3 # the conc. of those species are required in the computation of lifetime
concsSpecies= {'OH': 17.,'O3': 48. ,'NO3': 62.,'CO': 28.,'SO2': 64.,'NO': 30.,'NO2': 46.,'HO2': 33.,'RO2': 250.}

# species that does not participate in the redution
NokeepSp = ['SULF', 'NH3', 'HNO3', 'HCL', 'HONO', 'N2O5', 'RO2',
             'BiMT', 'SOAlP', 'NO2', 'NO3', 'H2O2', 'H2', 'NO',
             'CO', 'SO3', 'SO2', 'O3', 'OH', 'HO2','MD','BC','H2O']

lats = np.load('{:s}/latitudes.npy'.format(pathInitFiles))
lons = np.load('{:s}/longitudes.npy'.format(pathInitFiles))

pathSSH_sav=['files/ssh-aerosol-genoa/', \
             'files/ssh-files/'] # path to SSH and used files for different simulation modes

