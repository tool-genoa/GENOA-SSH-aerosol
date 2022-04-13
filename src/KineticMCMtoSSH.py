# -*- coding: utf-8 -*-
#================================================================================
#
#     GENOA v1.0: the GENerator of reduced Organic Aerosol mechanism
#
#     Copyright (C) 2022 CEREA (ENPC) - INERIS.
#     GENOA is distributed under GPL v3.
#
#================================================================================

import math
from numpy import pi

from Functions import isfloat,isint,multiplyList

#all the simple kinetic rates proposed by MCM
KRO2NO3=['2.3E-12']
KDEC=['1.00E+06']
K298CH3O2=['3.5E-13']
KRO2NO=['2.7E-12','360']    #2.7D-12*EXP(360/TEMP) 
KRO2HO2=['2.91E-13','1300'] #2.91D-13*EXP(1300/TEMP)
KAPHO2=['5.2E-13','980']    #5.2D-13*EXP(980/TEMP) 
KAPNO=['7.5E-12','290']     #7.5D-12*EXP(290/TEMP) 
KNO3AL=['1.4E-12','-1860']  #1.44D-12*EXP(-1862/TEMP) 
KROPRIM=['2.50E-14','-300'] #2.50D-14*EXP(-300/TEMP)
KROSEC=['2.50E-14','-300']  #2.50D-14*EXP(-300/TEMP) 
KCH3O2=['1.03E-13','365']   #1.03D-13*EXP(365/TEMP)
K14ISOM1=['3.00E7','-5300'] #3.00D7*EXP(-5300/TEMP) 
SKR_list=['KRO2NO3','KRO2NO','KRO2HO2','KAPHO2','KAPNO','KNO3AL','KROPRIM','KROSEC','KCH3O2','K298CH3O2','KDEC','K14ISOM1']

# complex kinetic rates proposed by MCM, need to further extend to all CKRs !!! 
# CKR=[KMT13,KFPAN,KMT14,KBPAN] # BCARY
# http://mcm.leeds.ac.uk/MCMv3.3.1/parameters/complex.htt
KFPAN = ['MCM1',3.28E-28,300,6.87, 1.125E-11,300,1.105, 0.30] # TEMP/k
KBPAN = ['MCM2',1.10E-5,-10100,1, 1.9E17,-14100,1, 0.30] # EXP(k/TEMP)
KBPPN = ['MCM2',1.7E-3,-11280, 1, 9.3E16, -13940, 1.0, 0.36] # EXP(k/TEMP)
KMT13 = ['MCM1',2.5E-30, 300,-5.5, 1.8E-11, 0., 0., 0.36] # TEMP/k
KMT14 = ['MCM2',9.0E-5, -9690, 0.0, 1.1E16, -10560, 0., 0.36] # EXP(k/TEMP)
KMT15 = ['MCM1',8.6E-29, 300, -3.1, 9.0E-12, 300, -0.85, 0.48] # TEMP/k

CKR_list=['KMT13', 'KFPAN', 'KMT14', 'KBPAN', 'KBPPN', 'KMT15']

manual_lists={
            '2*(K298CH3O2*8.0E-12)@0.5*RO2':'KINETIC TB RO2 ARR1 3.34664e-12',
            '2*(K298CH3O2*2.9E-12*EXP(500/TEMP))@0.5*RO2':'KINETIC TB RO2 ARR2 2.01494e-12 -250',
            '2*(KCH3O2*7.8E-14*EXP(1000/TEMP))@0.5*RO2':'KINETIC TB RO2 ARR2 2.01494e-12 -250',
            '(KCH3O2*7.8E-14*EXP(1000/TEMP))@0.5' :'KINETIC TB RO2 ARR2 8.963e-14 -682.5',
            '2*KCH3O2*RO2*7.18*EXP(-885/TEMP)':'KINETIC TB RO2 ARR2 1.4791e-12 520',
            # NC12H26
            '2*(KCH3O2*6.4E-14)@0.5*RO2':'KINETIC TB RO2 ARR2 1.6238e-13 -182.5', #2*(1.03E-13*6.4E-14*math.exp(365/TEMP))**0.5 => 1.6238e-13 182.5
            '2*(K298CH3O2*3E-13)@0.5*RO2':'KINETIC TB RO2 ARR1 6.4807e-13',# 2*(3.5E-13*3E-13)**0.5
            '2*(KCH3O2*1.6E-12*EXP(-2200/TEMP))@0.5*RO2':'KINETIC TB RO2 ARR2 1.6238e-12 917.5', # 2 * (1.03E-13*1.6E-12)**0.5 1.6238e-12 (365-2200)*0.5 = -917.5
                }

# Third party species
TBs=['M','O2','N2','H2O','H2','RO2']

def kconstant(val):
    return eval(val)

def dfappend(df,ind,val):
    #['TB','C1','EXP','CKR','unknown']
    df[1][ind].append(val)

def isTB(kin):
    """ third body used in SPACK"""
    if kin in TBs: return True
    else: return False

def tryEXP(kin):
    """return C2 if fits format EXP(C2/TEMP)"""
    if 'EXP' in kin: val=kin.replace('EXP(','').replace('/TEMP)','')
    elif 'exp' in kin: val=kin.replace('exp(','').replace('/TEMP)','')
    else: return False

    if isfloat(val): return val
    else: return False

def trySKR(kin):
    #SKR=[KRO2NO3,KRO2NO,KRO2HO2,KAPHO2,KAPNO,KNO3AL,KDEC,KROPRIM,KROSEC,KCH3O2,K298CH3O2,K14ISOM1]
    if kin in SKR_list: return eval(kin)
    else: return False

def tryCKR(kin):
    """complex kinetic rates proposed by MCM, need to further extend to all CKRs in MCM""" 
    #CKR=[KMT13,KFPAN,KMT14,KBPAN]
    if kin in CKR_list: return eval(kin)
    else: return False

def manual(kin):
    """return complete kinetic used in spack"""
    #The Reactions of RO2 with HO2
    if kin=='3.8E-13*EXP(780/TEMP)*(1-1/(1+498*EXP(-1160/TEMP)))':
        #CH3O2 + HO2 -> CH3OOH
        return 'KINETIC SPEC -10'
    elif kin=='3.8E-13*EXP(780/TEMP)*(1/(1+498*EXP(-1160/TEMP)))':
        #CH3O2 + HO2 -> HCHO
        return 'KINETIC SPEC -11'
    elif kin=='2*KCH3O2*RO2*(1-7.18*EXP(-885/TEMP))':
        #CH3O2 + RO2 -> CH3OH + HCHO
        #1.03E-13*math.exp(365/TEMP)*(1-7.18*EXP(-885/TEMP))
        return 'KINETIC TB RO2 SPEC -12'

    # for seperate
    elif kin=='2*KCH3O2*RO2*0.5*(1-7.18*EXP(-885/TEMP))':
        #CH3O2 + RO2 -> CH3OH + HCHO
        #1.03E-13*math.exp(365/TEMP)*(1-7.18*EXP(-885/TEMP))
        return 'KINETIC TB RO2 SPEC -15'

    # in case that lump CH3COCH3
    #elif kin=='8.8E-12*EXP(-1320/TEMP) + 1.7E-14*EXP(423/TEMP)':
    elif '8.8E-12*EXP(-1320/TEMP)+1.7E-14*EXP(423/TEMP)' in kin:
        if kin=='8.8E-12*EXP(-1320/TEMP)+1.7E-14*EXP(423/TEMP)':
            #CH3COCH3 + OH -> CH3COCH2O2
            return 'KINETIC SPEC -13'
        else:
            print('KtS: add ratio as TB in SSH: ',kin)
            a = eval(kin.replace('8.8E-12*EXP(-1320/TEMP) + 1.7E-14*EXP(423/TEMP)*',''))
            a = 'KINETIC TB {:6.3E} SPEC -13'.format(a)
            print(a)
            return a
    elif kin=='5.00E-12*O2*3.2*(1-EXP(-550/TEMP))':
        #HCOCO -> 1.0000 CO + 1.0000 OH
        #5.00E-12*O2*3.2*(1-EXP(-550/TEMP))
        return 'KINETIC TB O2 SPEC -14'

    # add for LIMONENE species: INDO
    elif kin=='1.80E+13*(TEMP/298)@1.7*EXP(-4079/TEMP)':
        #return 'KINETIC SPEC -16'
        return 'KINETIC ARR3 1.119708E+09 1.7 4079.'
    elif kin=='1.80E+13*(TEMP/298)@1.7*EXP(-4733/TEMP)':
        #return 'KINETIC SPEC -17'
        return 'KINETIC ARR3 1.119708E+09 1.7 4733.'
    elif kin=='2.20E+10*EXP(-8174/TEMP)*EXP(1.00E+8/TEMP@3)':
        return 'KINETIC SPEC -16'
    elif kin=='8.14E+9*EXP(-8591/TEMP)*EXP(1.00E+8/TEMP@3)':
        return 'KINETIC SPEC -17'

    else:
        for key,val in manual_lists.items(): # process ratio if exists 
            if key in kin:
                if key == kin: return val
                elif (key+'*') in kin:
                    #print(kin)
                    try:
                        tmp  = eval(kin.replace(key+'*',''))
                        if isfloat(tmp):
                            if 'ARR1' in val: return val.replace(val.split(' ')[-1],str(tmp*float(val.split(' ')[-1])))
                            elif 'ARR2' in val: return val.replace(val.split(' ')[-2],str(tmp*float(val.split(' ')[-2])))
                            else: return 0
                    except:
                        return 0
                else: return 0
    return 0

def mcm_photolysis(n):
    #J = l*np.cos(X)**m*EXP(-n*(1/np.cos(X))) 
    if n==1:   return [6.073E-05,1.743,0.474]
    elif n==2: return [4.775E-04,0.298,0.080]
    elif n==3: return [1.041E-05,0.723,0.279]
    elif n==4: return [1.165E-02,0.244,0.267]
    elif n==5: return [2.485E-02,0.168,0.108]
    elif n==6: return [1.747E-01,0.155,0.125]
    elif n==7: return [2.644E-03,0.261,0.288]
    elif n==8: return [9.312E-07,1.230,0.307]
    elif n==11:return [4.642E-05,0.762,0.353]
    elif n==12:return [6.853E-05,0.477,0.323]
    elif n==13:return [7.344E-06,1.202,0.417]
    elif n==14:return [2.879E-05,1.067,0.358]
    elif n==15:return [2.792E-05,0.805,0.338]
    elif n==16:return [1.675E-05,0.805,0.338]
    elif n==17:return [7.914E-05,0.764,0.364]
    elif n==18:return [1.140E-05,0.396,0.298]
    elif n==19:return [1.140E-05,0.396,0.298]
    elif n==20:return [7.600E-04,0.396,0.298]
    elif n==21:return [7.992E-07,1.578,0.271]
    elif n==22:return [5.804E-06,1.092,0.377]
    elif n==23:return [1.836E-05,0.395,0.296]
    elif n==24:return [1.836E-05,0.395,0.296]
    elif n==31:return [6.845E-05,0.130,0.201]
    elif n==32:return [1.032E-05,0.130,0.201]
    elif n==33:return [3.802E-05,0.644,0.312]
    elif n==34:return [1.537E-04,0.170,0.208]
    elif n==35:return [3.326E-04,0.148,0.215]
    elif n==41:return [7.649E-06,0.682,0.279]
    elif n==51:return [1.588E-06,1.154,0.318]
    elif n==52:return [1.907E-06,1.244,0.335]
    elif n==53:return [2.485E-06,1.196,0.328]
    elif n==54:return [4.095E-06,1.111,0.316]
    elif n==55:return [1.135E-05,0.974,0.309]
    elif n==56:return [7.549E-06,1.015,0.324]
    elif n==57:return [3.363E-06,1.296,0.322]
    elif n==61:return [7.537E-04,0.499,0.266]
    else:
        print('non recognize input mcm photolysis index: '+n)
        return []

def MCMtoSSH_photolysis(Jin,m=1):
    """ input: photolysis kinetic rate in MCM
        output: photolysis kineitc rate can be read by SPACK"""
    #photolysis index
    n = Jin.split('*')[0].replace('J<','').replace('>','')
    if isint(n): tmp = mcm_photolysis(int(n))
    else: tmp = []
    if len(tmp) == 3:
        #init output
        #outspack='KINETIC PHOTOLYSIS '
        #for SZA in tab_SZA: 
        #    outspack+=(' '+str(m*(a*(math.cos(SZA*pi/180.))**b*math.exp(-c*(1./math.cos(SZA*pi/180.))))))
    #SET TABULATION 11 DEGREES 0. 10. 20. 30. 40. 50. 60. 70. 78. 86. 90.
    #tab_SZA=[0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,78.0,86.0,90.0]
        return 'KINETIC MCM3 {:6.3E} {:5.3f} {:5.3f}'.format(tmp[0]*m,tmp[1],tmp[2])
    else:
        return 'KINETIC photolysis {:s} !!!'.format(Jin)

def MCMtoSSH_kinetic_rate(kin):
    """ output kinetic in string"""
    if ' ' in kin: kin = kin.replace(' ','')
    # photolysis
    if kin[0]=='J':
        # check factor
        if len(kin.split('*')) > 1: # example: J<22>*3
            try:
                m = eval(kin.split('*',1)[1])
                return MCMtoSSH_photolysis(kin.split('*')[0],m)
            except SyntaxError:
                print("non recognize input photolysis after '*': "+Jin)
                return 'KINETIC photolysis non recognize after * !!!'
        else:   return MCMtoSSH_photolysis(kin.split('*')[0])
    # manual
    elif manual(kin):
        return manual(kin)
    else:
        # build a df to store info
        df=[['TB','C1','EXP','CKR','unknown'],[[],[],[],[],[]]]

        # seperate ARR3
        if 'TEMP**' in kin: 
            tag_ARR3 = kin.split('TEMP**')[1].split('*')[0]
            kin = kin.replace('*TEMP**{:s}'.format(tag_ARR3),'')
        else: tag_ARR3 = False
        a_kin=[i.replace(' ','') for i in kin.split('*')]
        for i in a_kin:
            # float
            if isfloat(i):dfappend(df,1,float(i))
            elif isTB(i): dfappend(df,0,i)
            else:
                # exp(c2/T)
                tmp=tryEXP(i)
                if tmp:dfappend(df,2,float(tmp))
                else:
                    # SKR
                    tmp=trySKR(i)
                    if tmp:
                        if len(tmp)==1:dfappend(df,1,float(tmp[0]))
                        else: 
                            dfappend(df,1,float(tmp[0]))
                            dfappend(df,2,float(tmp[1]))
                    else:
                        tmp=tryCKR(i)
                        if tmp: dfappend(df,3,tmp)
                        else:
                            try:
                                dfappend(df,1,float(eval(i)))
                            except (NameError,SyntaxError):
                                dfappend(df,4,i)
        # has unknown part
        if df[1][4]:
            print(df,a_kin)
            print('KtS: ----- unkonwn: '+kin)
            return 'KINETIC !!!'
        #print(df)

        line='KINETIC'
        # TB
        if df[1][0]!=[]:
            if len(df[1][0])>1:
                print('KMS: More than one TB is found in ',kin)
                raise ValueError('number of TB should be <= 1')
            else: line+=' TB '+df[1][0][0]
        # C1
        if df[1][1]!=[]:
            tmp=multiplyList(df[1][1])
            if df[1][2] and df[1][3]: 
                print('KMS: Find both exp and ckr: '+kin)
                raise ValueError
            elif not df[1][2] and not df[1][3]: line+=' ARR1 {:6.3E}'.format(tmp)
            elif df[1][2] and not df[1][3]: 
                if tag_ARR3: line+=' ARR3 {:6.3E} {:s} {:8.2f}'.format(tmp,tag_ARR3,-1.*sum(df[1][2])) # negative value
                else: line+=' ARR2 {:6.3E} {:8.2f}'.format(tmp,-1.*sum(df[1][2])) # negative value
            elif not df[1][2] and df[1][3]: 
                if len(df[1][3])==1:
                    line += (' '+df[1][3][0][0])
                    for l in df[1][3][0][1:]: line+=(' {:6.3E}'.format(l))
                    line += ' {:6.4E}'.format(tmp)
                #raise ValueError('both c1 and ckr: '+kin)
        # CKR
        elif len(df[1][3])==1:
            line += (' '+df[1][3][0][0])
            for l in df[1][3][0][1:]: line+=(' {:6.3E}'.format(l))
            line += ' 1.0'
        else: line+=' !!!'

        return line

def KineticSSHtoStr(SSHstr):
    """update str from SSH format"""

    if 'KINETIC' not in SSHstr:
        print('MD.SSHtoStr: KINETIC is not found in the SSH string.', SSHstr)
        return None

    # special treatments
    if 'PHOTOLYSIS' in SSHstr: 
        return 'J<Tab>'

    # get elements
    tmp = [i for i in SSHstr.split(' ') if i != '']
    n = len(tmp)
    if n > 3: 
        if tmp[-1][0] == '-': val = tmp[-1][1:] # prepare for ARR2 and ARR3
        else: val = '-' + tmp[-1]

    if 'ARR1' in SSHstr:
        return tmp[-1]
    elif 'ARR2' in SSHstr:
        if n <= 3:
            raise ValueError('Check ARR2 type of reaction. '+SSHstr)
        else:
            return '{:s}*EXP({:s}/TEMP)'.format(tmp[-2],val)
    elif 'ARR3' in SSHstr:
        if n <= 3:
            raise ValueError('Check ARR3 type of reaction. '+SSHstr)
        else:
            return '{:s}*TEMP**{:s}*EXP({:s}/TEMP)'.format(tmp[-3],tmp[-2],val)
    else:
        print('MD.SSHtoStr: ',SSHstr)
        return SSHstr
