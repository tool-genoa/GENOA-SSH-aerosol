# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  KineticMCMtoPython.py converts the knietic rate contants from MCM format
#
#   to a format can be executed by Python. The converted format is used in 
#
#   reductions via lumping.
#
# ================================================================================

import math
import numpy as np

from Functions import cosx_max
from Parameters import concsSpecies,initfile, \
                       Tnow, prefix, NokeepSp
from KineticMCMtoSSH import mcm_photolysis
from SSHResultProcess import ref_conc

# FORMAT CONVERT
#MCMformat=['EXP','exp','@']#'D-','D+'
#Pyformat=['np.exp','np.exp','**'] # 'E-','E+'

def update_kinetic_rate(locs, IDchem_fake, path_fake, pathInitFiles, RefConcRead):
    #IDchem_fake, path_fake = pathNewRes, pathInitFiles = pathInitFiles, RefConcRead = pathNewRes+'/'+IDchem_fake):
    """get reference concentrations used for reduction from the conditions given by locs"""

    # default CONCENTRATIONS for TBs ['M','O2','N2','H2O','H2','RO2'] # Unit: molecule/cm3
    # M=2.55E+19, O2=0.2095*M, N2=0.7809*M, H2=0.0

    global COSX, SECX, RO2, H2O, TEMP, O2, M

    # get files in lists
    meteo_path = []
    gas_path = []
    refconc_paths = []

    for il in locs:

        if initfile in ['month','storage']: lc = 'm{:d}y{:d}x{:d}'.format(il[2],il[0],il[1])
        else: lc = 'y{:d}x{:d}'.format(il[0],il[1])

        if initfile == 'storage':
            lc = 'm{:d}/y{:d}/x{:d}'.format(il[2],il[0],il[1])
            meteo_path.append('{:s}/{:s}/meteo.dat'.format(pathInitFiles,lc))
            gas_path.append('{:s}/{:s}/gas.cst'.format(pathInitFiles,lc))
        else:
            meteo_path.append('{:s}/meteo/meteo_{:s}.dat'.format(pathInitFiles,lc))
            gas_path.append('{:s}/gas/gas_{:s}.cst'.format(pathInitFiles,lc))


    for il in locs:
        if initfile in ['month','storage']: lc = 'm{:d}y{:d}x{:d}'.format(il[2],il[0],il[1])
        else: lc = 'y{:d}x{:d}'.format(il[0],il[1])

        for j in Tnow: 
            refconc_paths.append('{:s}/{:s}/c_432000_3600_{:d}h_{:s}'.format(path_fake, lc,j, IDchem_fake.replace(prefix,'')))

    # get ref conc and save
    if RefConcRead:
        refconc_file = RefConcRead + str(len(locs))
        ref_conc(refconc_paths,species='all',Type='sum',radical='FA',savname=refconc_file,write=True)#, slist = speciesfile.replace('.mol','.species'))
    else:
        refconc_file = None
        
    # get concs
    concs_ave, concs_min = get_gas_cst(gas_path, concsSpecies)
    # get average TEMP, pressure, RH
    meteo_ave = get_daily_meteo(meteo_path, ['Temperature','Pressure','RelativeHumidity'])

    # useful variables
    TEMP = meteo_ave['Temperature']
    RH = meteo_ave['RelativeHumidity']

    # water
    H2O = 6.1078*math.exp(-1.0E0*(597.3-0.57*(TEMP-273.16))*18./1.986*(1./TEMP-1./273.16))*10./(1.38E-16*TEMP)*RH
    M = meteo_ave['Pressure'] * 7.243E16 / TEMP # air mass
    O2=0.2095*M
    if 'RO2' in list(concs_ave.keys()):
        RO2 = concs_ave['RO2']
    else:
        RO2 = 0.0

    # get solar zenith angle
    COSX,SECX = cosx_max(locs)

    # print check
    if len(locs) <= 10: print('--- KMP: locs:', locs)
    else: print('--- KMP: number of locs:', len(locs))

    print(meteo_ave, 'H2O, M, O2: ', H2O, M, O2)
    print('COSX,SECX: ', COSX,SECX)

    MCM_kinetic_rate(TEMP, M, O2, H2O)

    return concs_ave, concs_min, refconc_paths, refconc_file #, meteo_ave

def MCM_kinetic_rate(TEMP, M, O2, H2O):#(TEMP, M = 2.6E19, O2 = 5.5E18, H2O = 2E15):

    global KRO2NO3, KRO2NO, KRO2HO2, KAPHO2, KAPNO, KNO3AL, KDEC, \
           KROPRIM, KROSEC, KCH3O2, K298CH3O2, KBPAN, KFPAN, \
           KMT01, KMT02, KMT03, KMT04, KMT05, KMT06, KMT07, KMT08,\
           KMT09, KMT10, KMT11, KMT12, KMT13, KMT14, KMT15, KMT16,\
           KMT17, KMT18, KPPN0, KPPNI, KRPPN, FCPPN, NCPPN, FPPN, KBPPN,\
           K14ISOM1

    # simple kinetic rate
    KRO2NO3 = 2.3E-12
    KRO2NO = 2.7E-12*math.exp(360/TEMP)
    KRO2HO2 = 2.91E-13*math.exp(1300/TEMP)
    KAPHO2 = 5.2E-13*math.exp(980/TEMP)
    KAPNO = 7.5E-12*math.exp(290/TEMP)
    KNO3AL = 1.4E-12*math.exp(-1860/TEMP)
    KDEC = 1.00E+06
    KROPRIM = 2.50E-14*math.exp(-300/TEMP)
    KROSEC = 2.50E-14*math.exp(-300/TEMP)
    KCH3O2 = 1.03E-13*math.exp(365/TEMP)
    K298CH3O2 = 3.5E-13
    K14ISOM1 = 3.00E7*math.exp(-5300/TEMP) 

    # complex kinetic rate
    KD0 = 4.90E-3*math.exp(-12100/TEMP)*M
    KDI = 5.4E+16*math.exp(-13830/TEMP)
    KRD = KD0/KDI
    FCD = 0.30
    NCD = 0.75-1.27*(math.log10(FCD))
    FD = 10**(math.log10(FCD)/(1+(math.log10(KRD)/NCD)**2))
    KBPAN = (KD0*KDI)*FD/(KD0+KDI)

    KC0 = 2.7E-28*M*(TEMP/300)**-7.1
    KCI = 1.2E-11*(TEMP/300)**-0.9
    KRC = KC0/KCI
    FCC = 0.30
    NC = 0.75-1.27*(math.log10(FCC))
    FC = 10**(math.log10(FCC)/(1+(math.log10(KRC)/NC)**2))
    KFPAN = (KC0*KCI)*FC/(KC0+KCI)

    K10 = 1.0E-31*M*(TEMP/300)**-1.6
    K1I = 3.00E-11*(TEMP/300)**0.3
    KR1 = K10/K1I
    FC1 = 0.85
    NC1 = 0.75-1.27*(math.log10(FC1))
    F1 = 10**(math.log10(FC1)/(1+(math.log10(KR1)/NC1)**2))
    KMT01 = (K10*K1I)*F1/(K10+K1I)

    K20 = 1.3E-31*M*(TEMP/300)**-1.5
    K2I = 2.3E-11*(TEMP/300)**0.24
    KR2 = K20/K2I
    FC2 = 0.6
    NC2 = 0.75-1.27*(math.log10(FC2))
    F2 = 10**(math.log10(FC2)/(1+(math.log10(KR2)/NC2)**2))
    KMT02 = (K20*K2I)*F2/(K20+K2I)

    K30 = 3.6E-30*M*(TEMP/300)**-4.1
    K3I = 1.9E-12*(TEMP/300)**0.2
    KR3 = K30/K3I
    FC3 = 0.35
    NC3 = 0.75-1.27*(math.log10(FC3))
    F3 = 10**(math.log10(FC3)/(1+(math.log10(KR3)/NC3)**2))
    KMT03 = (K30*K3I)*F3/(K30+K3I)

    K40 = 1.3E-3*M*(TEMP/300)**-3.5*math.exp(-11000/TEMP)
    K4I = 9.7E+14*(TEMP/300)**0.1*math.exp(-11080/TEMP)
    KR4 = K40/K4I
    FC4 = 0.35
    NC4 = 0.75-1.27*(math.log10(FC4))
    F4 = 10**(math.log10(FC4)/(1+(math.log10(KR4)/NC4)**2))
    KMT04 = (K40*K4I)*F4/(K40+K4I)
    KMT05 = 1.44E-13*(1+(M/4.2E+19))
    KMT06 = 1 + (1.40E-21*math.exp(2200/TEMP)*H2O)

    K70 = 7.4E-31*M*(TEMP/300)**-2.4
    K7I = 3.3E-11*(TEMP/300)**-0.3
    KR7 = K70/K7I
    FC7 = math.exp(-TEMP/1420)
    NC7 = 0.75-1.27*(math.log10(FC7))
    F7 = 10**(math.log10(FC7)/(1+(math.log10(KR7)/NC7)**2))
    KMT07 = (K70*K7I)*F7/(K70+K7I)

    K80 = 3.3E-30*M*(TEMP/300)**-3.0
    K8I = 4.1E-11
    KR8 = K80/K8I
    FC8 = 0.4
    NC8 = 0.75-1.27*(math.log10(FC8))
    F8 = 10**(math.log10(FC8)/(1+(math.log10(KR8)/NC8)**2))
    KMT08 = (K80*K8I)*F8/(K80+K8I)

    K90 = 1.8E-31*M*(TEMP/300)**-3.2
    K9I = 4.7E-12
    KR9 = K90/K9I
    FC9 = 0.6
    NC9 = 0.75-1.27*(math.log10(FC9))
    F9 = 10**(math.log10(FC9)/(1+(math.log10(KR9)/NC9)**2))
    KMT09 = (K90*K9I)*F9/(K90+K9I)

    K100 = 4.10E-05*M*math.exp(-10650/TEMP)
    K10I = 4.8E+15*math.exp(-11170/TEMP)
    KR10 = K100/K10I
    FC10 = 0.6
    NC10 = 0.75-1.27*(math.log10(FC10))
    F10 = 10**(math.log10(FC10)/(1+(math.log10(KR10)/NC10)**2))
    KMT10 = (K100*K10I)*F10/(K100+K10I)

    K1 = 2.40E-14*math.exp(460/TEMP)
    K3 = 6.50E-34*math.exp(1335/TEMP)
    K4 = 2.70E-17*math.exp(2199/TEMP)
    K2 = (K3*M)/(1+(K3*M/K4))
    KMT11 = K1 + K2

    K120 = 4.5E-31*M*(TEMP/300)**-3.9
    K12I = 1.3E-12*(TEMP/300)**-0.7
    KR12 = K120/K12I
    FC12 = 0.525
    NC12 = 0.75-1.27*(math.log10(FC12))
    F12 = 10**(math.log10(FC12)/(1.0+(math.log10(KR12)/NC12)**2))
    KMT12 = (K120*K12I*F12)/(K120+K12I)

    K130 = 2.5E-30*M*(TEMP/300)**-5.5
    K13I = 1.8E-11
    KR13 = K130/K13I
    FC13 = 0.36
    NC13 = 0.75-1.27*(math.log10(FC13))
    F13 = 10**(math.log10(FC13)/(1+(math.log10(KR13)/NC13)**2))
    KMT13 = (K130*K13I)*F13/(K130+K13I)

    K140 = 9.0E-5*math.exp(-9690/TEMP)*M
    K14I = 1.1E+16*math.exp(-10560/TEMP)
    KR14 = K140/K14I
    FC14 = 0.4
    NC14 = 0.75-1.27*(math.log10(FC14))
    F14 = 10**(math.log10(FC14)/(1+(math.log10(KR14)/NC14)**2))
    KMT14 = (K140*K14I)*F14/(K140+K14I)

    K150 = 8.6E-29*M*(TEMP/300)**-3.1
    K15I = 9.0E-12*(TEMP/300)**-0.85
    KR15 = K150/K15I
    FC15 = 0.48
    NC15 = 0.75-1.27*(math.log10(FC15))
    F15 = 10**(math.log10(FC15)/(1+(math.log10(KR15)/NC15)**2))
    KMT15 = (K150*K15I)*F15/(K150+K15I)

    K160 = 8E-27*M*(TEMP/300)**-3.5
    K16I = 3.0E-11*(TEMP/300)**-1
    KR16 = K160/K16I
    FC16 = 0.5
    NC16 = 0.75-1.27*(math.log10(FC16))
    F16 = 10**(math.log10(FC16)/(1+(math.log10(KR16)/NC16)**2))
    KMT16 = (K160*K16I)*F16/(K160+K16I)

    K170 = 5.0E-30*M*(TEMP/300)**-1.5
    K17I = 1.0E-12
    KR17 = K170/K17I
    FC17 = 0.17*math.exp(-51/TEMP)+math.exp(-TEMP/204)
    NC17 = 0.75-1.27*(math.log10(FC17))
    F17 = 10**(math.log10(FC17)/(1.0+(math.log10(KR17)/NC17)**2))
    KMT17 = (K170*K17I*F17)/(K170+K17I)
    KMT18 = 9.5E-39*O2*math.exp(5270/TEMP)/(1+7.5E-29*O2*math.exp(5610/TEMP))

    KPPN0 = 1.7E-03*math.exp(-11280/TEMP)*M
    KPPNI = 8.3E+16*math.exp(-13940/TEMP)
    KRPPN = KPPN0/KPPNI
    FCPPN = 0.36
    NCPPN = 0.75-1.27*(math.log10(FCPPN))
    FPPN = 10**(math.log10(FCPPN)/(1+(math.log10(KRPPN)/NCPPN)**2))
    KBPPN = (KPPN0*KPPNI)*FCPPN/(KPPN0+KPPNI)

    return None

def photolysis(num,COSX):

    # IN FORMAT J = L*COSX** M *EXP(-N*SECX) 
    #COSX,SECX=ZENITH(latitude,TIME)
    #print('COSX,SECX',COSX,SECX)
    # zenith angle is at 90 degree COSX = 1
    if COSX <= 0.0: return 0.0
    else: SECX = 1./COSX
    tmp = mcm_photolysis(num)
    if len(tmp) == 3: return tmp[0]*(COSX**(tmp[1]))*math.exp(-1*tmp[2]*SECX)
    else: return None

def kinetic_MCM_to_value(rateStr):
    """return a float value of kineitc rate"""

    TEMP = 245.

    # For TOL
    if 'Tab' in rateStr or 'TROE' in rateStr: return 0

    if 'EXP' in rateStr: rateStr = rateStr.replace('EXP','math.exp')
    elif 'exp' in rateStr: rateStr = rateStr.replace('exp','math.exp')

    if '@' in rateStr: rateStr = rateStr.replace('@','**')

    # check photolysis
    if 'J<' in rateStr:
        num = int(rateStr.split('J<')[1].split('>')[0])
        pVal = '{:5.3E}'.format(photolysis(num,1.))
        rateStr = rateStr.replace('J<{:d}>'.format(num),pVal)

    # check other element
    try:
        tmp = float(eval(rateStr))
        if tmp < 0: print('kinetic rate < 0',tmp,rateStr)
        return tmp
    except NameError as ValueError: # check
        print('KMP: kinetic_MCM_to_value not able to transfer MCM format to a float.',rateStr)

def get_gas_cst(paths, sps_MWs):
    """return the average value (in mole/cm3) of given species (with MWs) from given paths
       paths = ['.../gas_{:s}.cst']
    """
    concs = {}

    for i in paths:
        with open (i,'r') as f: # read cst files
            for j in f.read().splitlines():
                if j.split('\t')[0] in list(sps_MWs.keys()):# species name
                    sname = j.split('\t')[0]
                    s = sps_MWs[sname] # index in sps_MWs
                    val = j.split('\t')[1:]
                    ns = len(val)
                    tmp = np.zeros(ns)
                    for k in range(ns): tmp[k] = float(val[k]) * 6.022E23 * 1E-12/s # record conc.
                    if paths.index(i) != 0: 
                        concs[sname].append(tmp)
                    else:
                        concs[sname] = [tmp]

    concs_ave = {} # lumping , compute lifetime
    for i in NokeepSp:
        if i in list(concs.keys()): concs_ave[i] = np.average(concs[i])
        else: concs_ave[i] = 0.0

    # print to check
    for i in list(sps_MWs.keys()): 
        if i in list(concs_ave.keys()): 
            print('ave conc:',i,"{:.15f}".format(concs_ave[i]),concs_ave[i]/6.022E23/1E-12*sps_MWs[i])

    return concs_ave, concs

def get_daily_meteo(paths, items):
    """return the average value of given items from the given paths
    """
    out = {}
    tmp = [0.0] * len(items)
    for i in paths:
        with open (i,'r') as f: # read cst files
            info = f.read().splitlines()
            titles = info[0].split(',') # Time,Temperature,Pressure,RelativeHumidity
            for j in info[1:25]: # get the first day
                for k in range(len(items)):
                    s = titles.index(items[k]) # index
                    tmp[k] += float(j.split(',')[s])
    # average
    for i in range(len(items)):
        out[items[i]] = (tmp[i] / (24. * len(paths)))
        print('get_daily_meteo', items[i], out[items[i]])
    return out
