# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  Parameters.py contains reduction parameters and options used for the reduction.
#
#   This script also read the input configuration file to update user-defined 
#
#   parameters. 
#
# ================================================================================

import os
import sys
import json
import numpy as np

# get the number of processor
from multiprocessing import cpu_count

# for configuration file
if len(sys.argv) > 1: import configparser


from ast import literal_eval

                #try:
                #    value = literal_eval(value)
                #except (ValueError, SyntaxError):
                #    print(f"Error: Invalid value '{value}' for parameter '{param_name}'.")
                #    print("Please provide a valid float or boolean value.")
                #    exit(1)
                
# print to check
print('Running python file: ',sys.argv[0])

#### Mechanism related setting

# output mechanism name [prefix] + [IDchem]
IDchem = 'None'
# SOA precursors
primaryVOCs = ['APINENE','BPINENE','LIMONENE']
# prefix for SOA precursotr set
prefix = 'MT'

# Output a viz file if number of aerosols >= nviz
nviz = 0

# Default species and reaction files
speciesfile = None
reactionfile = None


### Directories

# Directory to SSH-aerosols
# GENOA will generate two folders [pathSSH]+'_rdc' [pathSSH]+'_ref' during training
pathSSH = '/../SSHs/test'
# Directory to read/save chemical mechanisms and record files
pathNewChem = '/../toSSH/'
# Directory to save SSH-aerosol simulation results
pathNewRes = '/../results/'
# repositories of MCM files and SSH model
speciestype = 'SSH'
reactiontype = 'SSH'


#### Reduction basic settings

# Number of processors used for reduction
# default use the max number of processors
ncpu = cpu_count()

# Reduction general settings
RunSets=[
         # 0: new chem sets
         {'NewChem':0,
           # redcution options
          'prereduction':0, # 1 to active prereduction
           # prereduction options
           # 0: inactive. Numbers are the attempt times
          'kdec':0,   # remove kedc reactions; 0/1
          'tau':0,    # remove reactions with lifetime >= tau (s); 0/1E0
          'bratio':0, # remove reactions with branching ratios >= bratio (< 1); 1E-3_1
          'gen':0,    # remove reactions with No.generation >= gen; 13
          'Psat_NVOC':0, # set 'Psat_SVOC':0, 1E-13
          'Psat_SVOC':0, # set 'Psat_SVOC':0, 1E-4
          'Psat_aero':0,
          'conc':0,      # {'aero':1E-5,'gas':0.0,'TM':0.0},
          'lump': 0},    # DU_1
            
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
          'items': ['fun_72'],#['ratio','fun_8,72','Kp_8,72'],
          'savpath':'/../results/graphs'
          },

         # 3: Setups for Auto-Training.py
         {
          # reuse chem, results files
          'Reuse_chem': False,
          
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

          # reduction strategies
          'strategy_types':['rm','jp','lp','rp','rs','da'],

          # species not do removed/lumped during training
          'frozenspecies': [], # related reaction/species not tested
          'keptspecies':[],    # keep in the reduced mechanisms
          
          # pre-testing locations: numer read from testing file or a separated file
          'nPreTest': 100, # need to give testing file/ or a list of locations

          # error tolerances
          'training_parameter_table': None, # files to input parameters for training
          # priority: input_file > config file > Parameter.py
          'err_ref':[0.01,0.02,0.03,0.04,0.04,0.06,0.06,0.08,0.08,0.10,0.10], # error tolerance compared to ref case
          'err_pre':[0.01,0.01,0.01,0.03,0.02,0.04,0.04,0.06,0.04,0.08,0.10], # error tolerance compared to pre case
          'delta_err': 10, # ratios >= (1) or real value (< 1)
          'err_pav':2,
          'err_rav':2,
          'err_limit':1, # limitation of error tolerances, input values larger than err_limit is considered as a ratio
          'try_ave_ref': 99., # tolerance on err_rav to stop the reduction cycle
          'try_max_ref': 99., # tolerance on err_ref to stop the reduction cycle
          
          # options specific for series reduction
          #'tag_pre_reduction': False,
          #'tag_hybrid': False, # use multi-processing in series reduction
          #'BranchRatio': [5E-2, 1E-1, 5E-1], 
          # options for both seires and parallel reductions
          #'tag_quick_mode': False,
          
          # options for parallel reductions
          # parallel or series reduction
          'tag_para': True,       # default
          'train_key':'species',  # reaction / species used to order reduction cycle
          'train_order':None,     # generate by train_key or read from file

          # options specific for parallel reductions
          'nrcn_rm1':   0, # number of reactions to active elementory-like treatment - only use if tag_rm1 = 2
          'nrcn_limit': 0, # limitation for number of candidate reductions via removing reactions
          'err_pre_testing': 0.02, # err limitation cncurrent errors - used when tag_pre_testing == 2
          'err_stage_change': 99., # err to change stage II to IV - used when tag_stage != 0

          # options can be read from table, for all: 0 - inactive, 1 - active
          'tag_rm1': 0,         # 2 check nrcn_rm1, if n <= nrcn_rm1, active
          'tag_stage': 0,       # 2 first active aerosol-orient treatment and then inactive  
          'tag_efficient': 0,   # num >= 1: nval 
          'tag_pre_testing': 2, # 0 - training, 1 - pre-testing, 2 - check errors to decide
          
          ## not tested for GENOA v2.0
          'restart': False,  # name of the restat file
          'ntree': 0, # number of reduction tree to traverse
          # Time limitation for a redcution - reduction stop when time >= tlim (s)
          'time_limit':0,    # in seconds - output a restart file
          # if use tag_reach for pre-testing dataset
          'tag_reach':0,
          'err_rav_reach':0.0,
         },

        # 4: Setups for Auto_Testing.py
        {
         # reduction set
         'Test_file': None,
         # reference
         'IDchemRef':None, # ref chem
         'refChemPath':None,
         'orgpath': 'ref', # if save organics_1.txt for multiple testing
         'loc_num': False, # test locs number e.g., 100 # not used if ind is provided.
         # error
         'err_max': 1., # error tolerance. Program exits if finds error >= err_max * 3.
         'err_out': 1.,
         # name suffix to save err file
         'savname': None
        }
        ]

# Lumping options
Roptions = {
            # order species X
            'OrderType' : 'reaction', # 'conc'; 'reaction': their positions in the reaction list
            # order lumpable species for X
            'PeerType' : 'svoc', # 'reaction'; 'all'; 'svoc': check psat for SVOC
            # check relationship - 1: no lumping if A can form B 
            'CheckLink' : True, # 'True' 'False'
            # specific treatment for RO2
            'RO2Treat' : None,  # 'fixed'
            # check lifetime diff: tau1/tau2
            'CheckTau' : 10,    # value of tolerance, float
            # check No.generation diff: abs(gen1-gen2)
            'CheckGen' : False, # generation
            # check psat diff: psat1/psat2
            'CheckPsat' : 2,    # Psat_atm for SVOC
            # check types of species
            # if ture, species with different types won't be lumped
            'CheckStructure' : {
                                # checked by name (valid for MCM species)
                                'PAN': True,'N/C': True,
                                'OOH': True,'CO3H': True, 'CO2H': True,
                                # checked by strucutre
                                'C=C': True, # with C=C
                                # checked by formula
                                'O': 3, # diff in numbers of oxygen O <= 3
                                'C': 2  # diff in numbers of carbon C <= 2
                                }
            }

### Aerosol criteria
AeroDict = {
            # Saturation vapor pressure computation with UManSysProp, Topping et al., 2016
            # VP for vapor ponit : 0 - Nannoolal et al., 2008, 1 - Myrdal and Yalkowsky, 1997
            # BP for boilling point: 0 - Nannoolal et al., 2004, 1 - Stein and Brown, 1994, 2 - Joback and Reid, 1987
            # 'evap': EVAPORATION of Compernolle et al. (2011)
            # 'simpol': SIMPOL.1 of Pankow and Asher (2008)
            'vpType': 'VP0BP0',

            'soapfile':None, # soapfile to generate functional group decomposition
            
            'Psat_aero':0.0, # non-condensable if Psat > Psat_aero or Psat = 0. (not computed)
            'Psat_non_volatile': 1E-99, # aerosols with Psat <= Psat_NVOC: non-volatile
            }


#### SSH-aerosol related settings

# SSH-aerosol compilation mode
# 'fast_compile' - default - only output total SOA concentrations
# 'complete_compile' - output concentrations for all species - take time
SSH_mode = 'fast_compile'

# SSH-aerosol default namelist
namelist_pre = './files/ssh-aerosol-files/namelist_sav'
# Default aerosol species list used in SSH-aerosol
species_list_aer_init = './files/ssh-aerosol-files/species-list-aer-genoa.dat'

# Simulation starting times (h) - can set from 0 to 23 
# default 0h and 12h
Tnow = [0,12]
# Simulation time lists
# total time - 5 days
Ttot = [432000]
# time step - 1 hour
DeltaT = [3600]
# Number of initial conditions - only used for MT reduction: ipvocs = 3
ipvocs = 0

# output species besides RO2, Organics
outputs = {'aero': [], # 'PH2O'
           'gas': []}  # primaryVOCs[0]


### initial conditions
# Tag to read/save conditions in the format: m[x]/y[x]/x[x] 
initfile = 'storage'
# Directory to initial concentrations. Folders ordered by m[X_m]/y[X_y]/x[X_x]/
pathInitFiles = '/../Initial/MT/5ppb'
# numpy files to the index of latitude and longitude
# please download
pathLats = '../examples/2015_lat_round.npy'
pathLons = '../examples/2015_lon_round.npy'

# Training dataset location
# in format [y,x,m] [lat,lon,month] # lat: 32-70 lon: -17-39.8 month 0-11
# the index is in the format [X_y,X_x,X_m], where X is the index of the identifiers
# X_y for latitude (relationship in latitudes.npy), X_x for longitude (longitudes.npy), and X_m for month (X_m = real month-1)
# resolution of the simulation domain: 32 N to 70 N with a step of 0.25 N, 17 W to 39.8 E with a step of 0.4 E
locs = [[0,0,0]]


#### update options/parameters from the config file if provided
tag_copy = 0 # if copy config file
if len(sys.argv) > 1:

    cfgfile = sys.argv[1]
    print('\n\nReading parameters from config file: ', cfgfile)

    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(cfgfile)

    # assign parameters directly in certain sections
    cfgkeys = list(cfg.keys())

    # assign action section
    inds_cfg = { # index in RunSets
                'training' : 3, 'testing' : 4,
                'simulation': 1, 'postprocessing': 2,
                'build-mechanism' : 0, 'prereduction': 0}

    # check keys
    for s in cfgkeys:
        if s == 'action': continue
        if s in list(inds_cfg.keys()): continue
        
        for item in cfg.items(s):
            i, j = item # get name and value
            try:
                exec('{:s} = {:s}'.format(i, j)) # assign
                print('\tRead variable: {:s} = {:s}'.format(i, j))
            except:
                raise NameError('Not able find parameter: {:s} in section [{:s}]'.format(i,s))


    for s in cfg['action']:
    
        if s not in list(inds_cfg.keys()):
            print('Not recognize action option: ',s, 'not in ', list(inds_cfg.keys()))
            continue

        elif cfg['action'][s] in ['True', '1', 'true']:
        
            if s not in cfgkeys:
                raise NameError('Section [{:s}] not found.'.format(s))
                
            k = inds_cfg[s]

            # active postprocessing
            if s == 'postprocessing': RunSets[k]['display'] = 1
            #print('update info in ',s,k, RunSets[k])

            # copy and save the config file to check
            elif s == 'training': tag_copy = 1

            # build new mechanism or conduct prereduction
            elif s in ['build-mechanism','prereduction']:
                RunSets[k]['NewChem'] = 1
                if s == 'prereduction': RunSets[k]['prereduction'] = 1

            # read parameters
            for item in cfg.items(s):
                i, j = item # get name and value
                        
                if s != 'build-mechanism' and i in list(RunSets[k].keys()):
                    try:
                        exec('RunSets[k][i] = {:s}'.format(j)) # assign
                        print('\tRead options: ',s,i,j)
                    except:
                        RunSets[k][i] = j # assign
                        
                # assign aerosol properties
                elif i in list(AeroDict.keys()):
                    try:
                        exec('AeroDict[i] = {:s}'.format(j)) # assign
                        print('\tRead aerosol properties: ',s,i,j)
                    except:
                        AeroDict[i] = j # assign
                        
                else:
                    try:
                        exec('{:s} = {:s}'.format(i, j)) # assign
                        print('\tRead variables {:s} = {:s}'.format(i, j))
                    except:
                        print('not able find parameter: ',i, 'in RunSets[{:d}]'.format(k), cfgfile)
                        raise NameError('cfg: RunSets')
else:
    cfgfile = None
    print('No input config file!')


### Checking inputs
# check the training conditions
if locs == []:
    # no input
    if cfgfile: print('The input locs is empty.', cfgfile)
    raise ValueError('Check locs in the config file')
# input as a file
elif isinstance(locs,str):
    print('Read locs from the entire list', locs)
    with open (locs,'r') as f: locs = json.loads(f.read())
    print('read: ',len(locs),' locs')
# input as a list
elif  isinstance(locs,list):
    if isinstance(locs[0],str):
        n = int(locs[1])
        print('Read the first ',n,' locs from the list', locs[0])
        with open (locs[0],'r') as f: locs = json.loads(f.read())[:n]
    else:
        print('Read locs list: ',locs)
# input type not recognized
else:
    raise TypeError("locs type not recognized.", locs)

# Check the number of assigned processors
if ncpu < 1: 
    raise ValueError('ncpu cannot less than 1.', ncpu)
print('The number of processors for multiprocessing :',ncpu)

# if parallel reduction, ncpu is better to be larger than No. of strategies
if RunSets[3]['tag_para'] and ncpu <= len(RunSets[3]['strategy_types']):
    raise ValueError('ncpu is too small for an efficient run. ncpu should > No.strategy', ncpu, len(RunSets['strategy_types']))


### Update parameters 
# settings for running SSH simulations
tail = IDchem.replace(prefix,'')
# add prefix
if prefix not in IDchem: IDchem = prefix + IDchem
ResultFolder = 'Results_{:s}'.format(IDchem)

if pathSSH[-1] == '/': pathSSH = pathSSH[:-1]

pathSSH_rdc = '{:s}_rdc'.format(pathSSH)
pathSSH_ref = '{:s}_ref'.format(pathSSH)
pathInitFiles = os.path.abspath(pathInitFiles)

if speciesfile == None: 
    speciesfile = '{0:s}/{1:s}/{1:s}.mol'.format(pathNewChem,prefix+IDchem)
if reactionfile == None:
    reactionfile = '{0:s}/{1:s}/{1:s}.reactions'.format(pathNewChem,prefix+IDchem)


#### Other default settings - can not be read from configuration file

# defalut tiny value
TINYM = 1E-99

### Added species list and their molar weight
# all species in this list will output in the species list of the generated SOA mechanism
SSHSpeciesInit={'SULF':98.,'NH3':17.,'HNO3':63.,
                'HCL':36.5,
                'RO2pool':250.,'SOAlP':392.,'BiMT':136.,
                'OH':17.,'NO3':62.,'O3':48.,'CO':28.,
                'NO2':46.,'NO':30.,'HO2':33.,'H2O2':34.,
                'H2':2., 'SO2':64., 'SO3':80.}

# Molar mass of the given chemical element
element_mass={'H':1.,'C':12.,'O':16.,'N':14.,'S':32.}

# basic SMART structures
keyels={'C':'[#6]', # carbon
        'H':'[H]', # hydrogen
        'O':'[O]', # oxygen 6
        'N':'[N]', # nitrogen 7
        #'S':'[S]'
       }

# Species use constant files in SSH-aerosol simulations
cst=['OH','NO3','O3','NO','NO2','HO2','CO']
# conc. # unit: ug/m3 # the conc. of those species are required in the computation of lifetime
concsSpecies= {'OH': 17.,'O3': 48. ,'NO3': 62.,'CO': 28.,'SO2': 64.,'NO': 30.,'NO2': 46.,'HO2': 33.,'RO2pool': 250.}

# Species that does not participate in the redution
NokeepSp = ['SULF', 'NH3', 'HNO3', 'HCL', 'HONO', 'N2O5', 'RO2pool',
             'BiMT', 'SOAlP', 'NO2', 'NO3', 'H2O2', 'H2', 'NO',
             'CO', 'SO3', 'SO2', 'O3', 'OH', 'HO2','MD','BC','H2O']

# list of latitude
lats = np.load(pathLats)
# list of longitude
lons = np.load(pathLons)
# path to a simplified version of SSH-aerosol v1.3, and files for different compilation modes
pathSSH_sav=['./files/ssh-aerosol-genoa', \
             './files/ssh-aerosol-files']
             
# copy the current configure
if tag_copy:
    # copy config file in [ID]_recs
    path_sav_rec = '{0:s}/{1:s}_recs'.format(pathNewChem,tail)
    os.system('mkdir -p {:s}'.format(path_sav_rec))
    os.system('cp {:s} {:s}/config_used.ini'.format(cfgfile, path_sav_rec))
    print('config file is copied to ', path_sav_rec)
