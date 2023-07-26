# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  GENOA.py is the main script for executing the GENOA program.
#
#  To run the program using the command: python3 GENOA.py [configuration_file]
#
# ================================================================================

import os
import sys
import configparser

if __name__ == '__main__':

    # print program info
    tmp =  '='*75+'\n\n'
    tmp += ' '*5+'GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism\n\n'
    tmp += ' '*5+'Copyright (C) 2023 CEREA (ENPC) - INERIS.\n'
    tmp += ' '*5+'GENOA is distributed under GPLv3.'+'\n\n'
    tmp += '='*75+'\n\n'
    print(tmp)

    # get job id
    n = os.getpid()
    print('Current jobs id:',n, flush=True)

    if len(sys.argv) <= 1:
        print('No configuration file is provied. Please enter: python3 GENOA.py [config_file]')
        raise FileNotFoundError
    else:
        cfg_file = sys.argv[1]
        if not os.path.exists(cfg_file):
            print('Not find the configuration file.', cfg_file)
            raise FileNotFoundError
    try:
        cfg = configparser.ConfigParser()
        cfg.optionxform = str
        cfg.read(cfg_file)
    except IndexError:
        print('Error: check the provied configuration.', cfg_file)
        raise

    try:
        action = cfg['action']
    except KeyError:
        print('Error: not found action section in the configuration file.', cfg_file)
        raise

    tag_action = 0 # 1: training; 2: testing; 3: reduction
    tags = [i for i in list(action.keys()) if action[i] in ['True', '1', 'true']] 

    if 'training' in tags:
    
        if 'testing' in tags:
            print('Run reduction ...')
        else:
            print('Run training ...')
        
        from AutoReduction import auto_reduction
        auto_reduction()
        
    elif 'testing' in tags:
    
        print('Run testing ...')
        
        from AutoTesting import auto_testing
        auto_testing(sav_soapath = True)
        
    else:
        from BuildRunPlot import build_run_plot
        build_run_plot()

