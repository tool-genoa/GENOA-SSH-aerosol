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
import configparser

# add the path to 'src/'
#if 'src' in os.listdir('.'): sys.path.append('src')

if __name__ == '__main__':

    # print program info
    tmp =  '='*75+'\n\n'
    tmp += ' '*5+'GENOA v1.0: the GENerator of reduced Organic Aerosol mechanism\n\n'
    tmp += ' '*5+'Copyright (C) 2022 CEREA (ENPC) - INERIS.\n'
    tmp += ' '*5+'GENOA is distributed under GPLv3.'+'\n\n'
    tmp += '='*75+'\n\n'
    print(tmp)

    if len(sys.argv) <= 1:
        print('No configuration file is provied. Please enter: python3 __init__.py [config_file]')
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

    tag_action, tag_run = 0, 0 # 1: training; 2: testing; 3: reduction
    tags = list(action.keys())   

    if 'training' in tags and action['training'] in ['True','1']:
        tag_action += 1
    if 'testing' in tags and action['testing'] in ['True','1']:
        tag_action += 2

    if tag_action == 3:
        print('Run reduction...')
        from AutoReduction import auto_reduction
        auto_reduction()
    elif tag_action == 2:
        print('Run testing...')
        from AutoTesting import auto_testing
        auto_testing()
    elif tag_action == 1:
        print('Run training...')
        from AutoTrainingSeries import auto_training_srs
        auto_training_srs()
    else:
        print('No actived action.')
