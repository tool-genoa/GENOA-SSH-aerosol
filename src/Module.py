# -*- coding: utf-8 -*-
#================================================================================
#
#     GENOA v1.0: the GENerator of reduced Organic Aerosol mechanism
#
#     Copyright (C) 2022 CEREA (ENPC) - INERIS.
#     GENOA is distributed under GPL v3.
#
#================================================================================

import json

from Functions import duplicates,trim_items, \
                      isfloat, compare
from Parameters import NokeepSp,AeroDict
from KineticMCMtoSSH import MCMtoSSH_kinetic_rate
from KineticMCMtoPython import kinetic_MCM_to_value

# import MolProperty as mp # not used

sav_species=['name','nameUsed','status',
             'lump','jump',
             'organic','RO2','Radical','condensed',
             'SMILES','InChI','mass','formula','source',
             'psat_torr','psat_atm','dHvap_KJ','gamma','henry','Kp',
             'temp','DU','ratios','functionalGroups','SOAPStructure','generation']

class Species:
    """record properties of a given species"""
    def __init__(self,name=None):
        self.name=name # Name
        self.nameUsed=None # if set, name changed in lumping

        # properties
        self.organic=False
        self.RO2=False 
        self.Radical=False
        self.condensed=False

        self.SMILES='-' # smiles
        self.InChI=None
        self.mass=0. #MW
        self.formula=''

        #reduction
        self.lump=[] # if set, not output
        self.jump=[]
        self.remove=[]

        # for SSH output
        self.type=4 #Type
        self.groupID=2 #group ID
        self.coll_fac=687.
        self.mole_diam=8.39
        self.surf_tens=30.e-03
        self.accomod=1.0
        self.mass_dens=1.30e-06
        self.non_volatile=0

        self.psat_torr=0.0 #Saturation vapor pressure at Tref (torr)
        self.psat_atm=0.0 #Saturation vapor pressure at Tref (atm)
        self.dHvap_KJ=0.0 #Enthalpy of vaporization (kJ/mol)
        self.gamma=None #Activity coefficient at infinite dilution in water
        self.henry=0.0 #Henry's Law constant (M/atm, mol/atm*l)
        self.Kp=None

        self.temp= 298. #Temperature of reference (K)

        self.functionalGroups={}
        self.ratios=None
        self.DU=None
        self.pymol=None
        self.status=True
        self.generation = -1

        self.SOAPStructure = None

        self.source = [] # default

    def add_info(self,key,value,check=False):
        if check:
            if not hasattr(self, key):
                print(self.name, key,value) 
                raise AttributeError('Species add properties: do not have this attr')
        setattr(self, key, value)
            
    def toText(self,filename=None):
        sav = []
        for i in sav_species:
            ctn = json.dumps(eval('self.'+i))
            sav.append(('{:s}\t{}\n').format(i,ctn))

        if filename == None: return sav
        else:
            f = open (filename,'a+')
            f.write('\n')
            for i in sav: f.write(i)
            f.write('\n')
            f.close()
            return None


    def fromTest(self,text):
        for i in text:
            tmp = i.split('\t')
            if len(tmp) == 2: 
                tmp[1] = json.loads(tmp[1])
                if isinstance(tmp[1],str): tmp[1] = str(tmp[1])
                self.add_info(tmp[0],tmp[1],True)

    def to_species_list_aer(self, Type = 'smiles', firstRun = False, henry = 0.0):
        if self.organic:
            spr='\t'
            line='P{:s}'.format(self.name)

            if Type == 'smiles': tmp = self.SMILES
            elif Type == 'vectors': tmp = self.VectorinText()
            else:
                raise NameError('MD: unknown type: '+Type)

            content=[self.type,self.groupID, self.mass,self.name,self.coll_fac,
                    self.mole_diam,self.surf_tens,self.accomod,self.mass_dens,
                    self.non_volatile,tmp,self.psat_torr,self.dHvap_KJ, henry]
            # remove precursor if it is the 1st run
            if firstRun: content[3] = '--'

            style=['{:d}','{:d}','{:6.2f}','{:s}','{:5.2E}',
                    '{:5.2E}','{:5.2E}','{:5.2E}','{:5.2E}',
                    '{:d}','{:s}','{:5.2E}','{:5.2E}', '{:5.2E}']
            for i in range(len(content)):
                tmp=spr+style[i].format(content[i])
                if style[i]=='E': tmp=tmp.replace('E','D')
                line+=tmp
            return line+'\n'

        
    def update(self,AeroCriteria):

        # check
        if not self.status: return None

        if self.SMILES == '-': 
            raise ValueError('MD Species Update: no smiles of the given species: '+self.name)
        elif '/' in self.SMILES or '\\' in self.SMILES:
            print('!!!!/,\\!!!!',self.name,self.SMILES)
            for i in '/','\\': self.SMILES = self.SMILES.replace(i,'')


        # pymol
        self.pymol = mp.get_pybelmol(self.SMILES)

        # formula
        self.formula = self.pymol.formula

        # mass
        self.mass=self.pymol.molwt

        # psat
        self.psat_atm = mp.saturation_vapor_pressure(self.pymol,AeroCriteria['vpType'],self.temp)
        self.psat_torr=self.psat_atm*760. # unit

        # non-volatile
        if self.psat_atm <= AeroCriteria['Psat_NVOC'] : 
            print('species ',self.name,' is non_volatile, with P_sat = ',self.psat_atm,' < ',AeroCriteria['Psat_NVOC'])
            #self.non_volatile = 1

        # Hvap
        self.dHvap_KJ = mp.enthalpy_vaporization(self.pymol,AeroCriteria['vpType'])/1000. # unit in KJ

        # functional group
        self.functionalGroups = mp.functional_groups(self.pymol)

        # is RO2? Radicals? Condensed species? 
        self.is_type(AeroCriteria['Psat_aer'])

        # SOAPstructure ?
        #self.SOAPStructure = mp.SOAPStructure(self.pymol)

    def update_advanced(self,AeroCriteria):

        # check
        if not self.status: return None

        if self.SMILES == '-': 
            raise ValueError('Species Update: no smiles of the given species: '+self.name)
        elif '/' in self.SMILES or '\\' in self.SMILES:
            print('!!!!/,\\!!!!',self.name,self.SMILES)
            for i in '/','\\': self.SMILES = self.SMILES.replace(i,'')

        # pymol
        self.pymol = mp.get_pybelmol(self.SMILES)

        # formula
        self.formula = self.pymol.formula

        # mass
        self.mass=self.pymol.molwt

        # psat
        self.psat_atm = mp.saturation_vapor_pressure(self.pymol,AeroCriteria['vpType'],self.temp)
        self.psat_torr=self.psat_atm*760. # unit

        # non-volatile
        if self.psat_atm <= AeroCriteria['Psat_NVOC'] :
            print('species ',self.name,' is non_volatile, with P_sat = ',self.psat_atm,' < ',AeroCriteria['Psat_NVOC'])
            #self.non_volatile = 1

        # Hvap
        self.dHvap_KJ = mp.enthalpy_vaporization(self.pymol,AeroCriteria['vpType'])/1000. # unit in KJ

        # functional group
        self.functionalGroups = mp.functional_groups(self.pymol)

        # is RO2? Radicals? Condensed species? 
        self.is_type(AeroCriteria['Psat_aer'])

        # SOAPstructure ?
        #self.SOAPStructure = mp.SOAPStructure(self.pymol)

        # activity coefficient
        self.gamma=mp.activity_coefficient_inf(self.pymol,self.temp)

        # henry's law constant (M/atm)
        # H = 1000*760/(18.0*GAMMAinf*Psat)
        self.henry = 1000/(18*self.gamma*self.psat_atm)

        # partitioning coefficient
        self.Kp=mp.partitioning_coefficient(self,self.temp)

        # atomic ratios and degree of unsaturation
        self.ratios,self.DU=mp.organic_ratios_and_du(self.pymol)

    def update_Psat(self,AeroCriteria):

        # check
        if not self.status or self.pymol is None: return None

        # psat
        self.psat_atm = mp.saturation_vapor_pressure(self.pymol,AeroCriteria['vpType'],self.temp)
        self.psat_torr=self.psat_atm*760. # unit

        # Hvap
        self.dHvap_KJ = mp.enthalpy_vaporization(self.pymol,AeroCriteria['vpType'])/1000. # unit in KJ

        # henry's law constant (M/atm)
        # H = 1000*760/(18.0*GAMMAinf*Psat)
        self.henry = 1000/(18*self.gamma*self.psat_atm)

    def update_mass(self):

        if self.pymol is not None: self.mass = self.pymol.molwt
        elif self.SMILES != '-': mol=self.SMILES
        else: mol=self.name

        # mass
        self.mass=mp.molar_mass(mol)

    def is_type(self,Psat):
        
        if self.pymol.spin>1 or '[O+]' in self.SMILES or '[O]' in self.SMILES:
            #print('spin==',self.pymol.spin,self.name,self.SMILES)
            self.Radical=True
            if self.functionalGroups == {} : self.functionalGroups = mp.functional_groups(self.pymol)
            #if 'RO2' in self.functionalGroups.keys(): 
            if self.name[-2:] == 'O2': self.RO2=True
            if self.name[-2:] not in ['OO','O2','O3']:
                if self.name[-1:] != 'O':
                    print('MD: radicals not in O, OO, O2, O3: ',self.name)

        if not self.Radical:
            if Psat ==0.0 or self.psat_atm <= Psat : 
                #print(self.name,self.psat_atm,Psat,self.psat_atm <= Psat)
                self.condensed=True

    def toSOAP(self,info):
        if not self.status or self.SOAPStructure is None: return False

        iStr = self.SOAPStructure
        tmp = ("  double group_tmp_"+self.name+" [] = {")
        tmp0 = [0.0] * 56
        for k in range(len(iStr[0])):
            tmp0[iStr[0][k]] = iStr[1][k]
        for k in range(len(tmp0)):
            if k == len(tmp0)-1: tmp+= (str(tmp0[k])+'};\n\n')
            elif k%6: tmp+= (str(tmp0[k])+', ')
            else: tmp+=(str(tmp0[k])+',\n          ')

        tmp_all = ''
        tag=0
        for j in info:
            if '***' in j:
                tag = info.index(j) 
                break
            else:
                k = j.replace('BiA2D',self.name)
                if 'Psat_ref' in k: tmp_all+=('  {:s}.{:s}={:5.2E};\n'.format(self.name,'Psat_ref',self.psat_torr))
                elif 'deltaH' in k: tmp_all+=('  {:s}.{:s}={:5.2E};\n'.format(self.name,'deltaH',self.dHvap_KJ))
                elif 'aq_type' in k:
                    if 37 in iStr[0]:
                    #if '-COOH' in self.functionalGroups.keys():
                        ind = iStr[1][iStr[0].index(37)]
                        if ind >= 2.0:
                        #if self.functionalGroups['-COOH'] >= 2:
                            tmp_all+=('  {:s}.aq_type="diacid";\n'.format(self.name))
                            tmp_all+=('  {:s}.Kacidity1=3.95e-4;\n'.format(self.name))
                            tmp_all+=('  {:s}.Kacidity2=7.70e-6;\n'.format(self.name))
                        elif ind <= 1.0:
                        #elif self.functionalGroups['-COOH'] >=1 :
                            tmp_all+=('  {:s}.aq_type="monoacid";\n'.format(self.name))
                            tmp_all+=('  {:s}.Kacidity1=6.52e-4;\n'.format(self.name))
                        else:
                            tmp_all+=('  {:s}.aq_type="none";\n'.format(self.name))
                    else: 
                        tmp_all+=('  {:s}.aq_type="none";\n'.format(self.name))
                elif 'nonvolatile' in k:
                    if self.psat_atm <= 1E-14: #AeroCriteria['Psat_NVOC']: 
                        tmp_all+=('  {:s}.nonvolatile=true;\n'.format(self.name))
                    else:
                        tmp_all+=('  {:s}.nonvolatile=false;\n'.format(self.name))
                else:tmp_all+=(k+'\n')

        tmp_all+=tmp
        if not tag:
            raise AttributeError('MD toSOAP: tag = 0')

        for j in info[tag:]:
            tmp_all+=(j.replace('BiA2D',self.name)+'\n')

        return tmp_all

    def VectorinText(self):
        tmp = ""        
        if self.SOAPStructure is not None:
            for i in range(56):
                if i in self.SOAPStructure[0]: tmp += '&{:5.2E}'.format(self.SOAPStructure[1][self.SOAPStructure[0].index(i)])
                else: tmp += '&{:5.2E}'.format(0.0)
        #else: raise AttributeError('species '+self.name+' has SOAPstructure None')
        if tmp != "": return tmp
        else: return "-"

class Reaction:
    """record a chemical reaction"""
    def __init__(self):
        self.reactants=[]
        self.ratiosRC=[]
        self.products=[] 
        self.ratiosPD=[]       
        self.rate=Kinetic()
        self.status=True
        self.type=None
        self.species=None

    def record_product(self,pd,rt):
        """record a single product and its ratio"""
        if pd in self.products:
            self.ratiosPD[self.products.index(pd)]+=rt
        else:
            self.products.append(pd)
            self.ratiosPD.append(rt)

    def record_reactant(self,rc,rt):
        """record a single reactant and its ratio"""
        if rc in self.reactants:
            self.ratiosRC[self.reactants.index(rc)]+=rt
        else:
            self.reactants.append(rc)
            self.ratiosRC.append(rt)

    def check_type(self):
        # get type from the number of reactants
        self.type = len(self.reactants)
        # trim species
        tmp = self.trim_reactants()
        tmp1 = self.trim_products()

        # check products
        if len(tmp1)<1 or len(tmp) == 0:
            # very basic reactions, not consider in the reduction
            self.type=-1
            self.species = [self.reactants,self.products]
        else:
            self.species = [tmp,tmp1]

    def trim_zero(self,val):
        """remove products if their ratios = 0.0"""
        sps= [self.reactants,self.products]
        rts= [self.ratiosRC,self.ratiosPD]

        for n in range(2):
            tag=1
            while tag:
                tag=0
                for i in range(len(rts[n])):
                    if rts[n][i] <= val:
                        del sps[n][i]
                        del rts[n][i]
                        tag=1
                        break

        # check if empty
        for i in sps:
            if i==[]: self.status=False

    def trim_reactants(self):
        return trim_items(self.reactants,['OH','O3','NO3','HO2','SO2','NO2','NO','CO','HCL','CL'])

    def trim_products(self):
        return trim_items(self.products,['HO2','SO3','SO2','NO3','NO2','NO','OH','H2O2','CO'])

    def merge_duplicates(self):
        """if repeat products or repeat reactants, remove one of them"""
        if not self.status: return False

        sps= [self.reactants,self.products]
        rts= [self.ratiosRC,self.ratiosPD]

        # inside check
        for n in range(2):
            i = 0
            while i < len(sps[n]):
                lists= duplicates(sps[n],sps[n][i])
                if len(lists) > 1:
                    for j in lists[1:]:
                        # add ratio
                        rts[n][lists[0]]+=rts[n][j]
                        # remove ratio and products
                        del sps[n][j]
                        del rts[n][j]
                i+=1

    def toSSH(self,Type='all'):
        """transfer to the SSH output"""
        # output in string
        rline=''

        # add reactants
        rcn=self.reactants
        rts=self.ratiosRC

        for j in range(len(rcn)):
            # output format
            if j==0: tmp=''
            elif j%3==2: tmp=' //\n     + '
            else: tmp=' + '

            # record reactant # if ratio is 1, not write ratio
            if rts[j]==1.0: rline+=('{:s}{:s}'.format(tmp,rcn[j]))
            else:
                rline+=('{:s}{:5.4f} {:s}'.format(tmp,rts[j],rcn[j]))
                print(rline)
                raise ValueError('MD: reactant ratio is not 1 ')

        # add arrow
        rline+=' -> '

        # add products and ratio
        rcn=self.products
        rts=self.ratiosPD
        for j in range(len(rcn)):
            # output format
            if j==0: tmp=''
            elif j%3==2: tmp=' //\n     + '
            else: tmp=' + '
            # if ratio is 1, not write ratio
            if rts[j]==1.0: rline+=('{:s}{:s}'.format(tmp,rcn[j]))
            else: rline+=('{:s}{:5.3E} {:s}'.format(tmp,rts[j],rcn[j]))

        if Type == 'simple': return rline
        elif Type == 'all':
            # add Kinetic rate: as comment and spack form
            rcn=self.rate
            for j in ('%'+rcn.str),rcn.SSH: rline+=('\n'+j)
        else:
            raise NameError('MD: type unknown '+Type)

        return rline

    def loadSSH(self,SSHstr):
        """build from the SSH string info"""
        self.status = True
        self.reactants=[]
        self.ratiosRC=[]
        self.products=[] 
        self.ratiosPD=[]       
        self.rate=Kinetic()

        for line in SSHstr.split('\n'):
            if '===' in line: continue # remove comment line
            elif '->' in line: # reactants and products
                # get all items
                items = [i for i in line.split(' ') if i != '']
                n = -1
                for j in items:
                    if j == '+' or j == '': continue
                    elif '->' == j : # n index of products
                        n = items.index(j) 
                        break
                    else: # record reactants
                        self.reactants.append(j)
                        self.ratiosRC.append(1.0)
                if n == -1: # check: n might be referenced before assignment
                    print(line)
                    raise ValueError('MD loadSSH: Not find reactant.')
                # add products
                tmp = 1.0 # branch ratio
                for j in items[n+1:]:
                    if j == '+' or j == '//': continue
                    elif isfloat(j):
                        tmp = float(j)
                    else:
                        self.products.append(j)
                        self.ratiosPD.append(tmp)
                        tmp = 1.0 # reset ratio
            elif line[0] == ' ': # read products
                # rest of products
                items = [i for i in line.split(' ') if i != '']
                tmp = 1.0
                for j in items:
                    if j == '+' or j == '//': continue
                    elif isfloat(j): tmp = float(j)
                    else: # record products
                        self.products.append(j)
                        self.ratiosPD.append(tmp)
                        tmp = 1.0
            # kinetic rate
            elif line[0] == '%':
                # kinetic rate in MCM format
                self.rate.str = line[1:]
                if 'J' in line: self.rate.Photolysis = True
            # kinetic rate
            elif 'KINETIC' in line:
                # kinetic rate in SSH format
                self.rate.SSH = line
        # check
        for j in self.reactants,self.ratiosRC:
            if j == []:
                print(self.reactants,self.ratiosRC)
                raise AttributeError('MD loadSSH: empty reactants or reactant ratio.')
            elif len(self.reactants) != len(self.ratiosRC):
                print(len(self.reactants) != len(self.ratiosRC))
                raise AttributeError('MD loadSSH: the numbers of reactants and ratios are not the same.')
            elif len(self.products) != len(self.ratiosPD):
                print(len(self.products) != len(self.ratiosPD))
                raise AttributeError('MD loadSSH: the numbers of products and ratios are not the same.')
            elif self.rate.SSH is None: 
                print(SSHstr)
                raise AttributeError('MD loadSSH: rate is empty')
            elif self.rate.str is None:
                # generate str format from SSH format
                self.rate.SSHtoStr()
    def isEqual(self,rcn):
        """return if this two reactions are identical"""
        tag = 1
        # check reactants and products
        a = [[self.reactants,self.ratiosRC],[self.products,self.ratiosPD]]
        b = [[rcn.reactants,rcn.ratiosRC],[rcn.products,rcn.ratiosPD]]

        for i in range(2):
            if a[i][0] == b[i][0]: #if self.reactants == rcn.reactants:
                # check ratios
                if a[i][1] != b[i][1]: #if self.ratiosRC != rcn.ratiosRC: tag = 0
                    tag = 0
                    break
            # sequences different
            elif sorted(a[i][0]) == sorted(b[i][0]): #sorted(self.reactants) == sorted(rcn.reactants):
                # check ratios one by one
                for j in range(len(a[i][0])): # for i in range(len(self.reactants)):
                    if a[i][1][j] != b[i][1][b[i][0].index(a[i][0][j])]:
                    #if self.ratiosRC[i] != rcn.ratiosRC[rcn.reactants.index(self.reactants[i])]: 
                        tag = 0
                        break
                if tag == 0: break
            else: #not same reactants 
                tag = 0
                break

        if tag == 0: return False

        # check kinetic rate
        a = [self.rate.str,self.rate.SSH]
        b = [rcn.rate.str,rcn.rate.SSH]

        for i in range(len(a)):
            if a[i] != b[i]:
                tag = 0
                break
        if tag == 0: return False
        else: return True

class Kinetic:
    """store kinetic in MCM string and ssh format"""

    def __init__(self,strin=None):

        self.str=strin
        self.pyformat=None
        self.Photolysis=False
        self.SSH=None

        if strin is not None: self.update()

    def update(self,strin=None):

        if strin is not None: self.str=strin

        self.SSH = MCMtoSSH_kinetic_rate(self.str)

        if self.str[0]=='J': self.Photolysis=True
        else: self.Photolysis=False

        #if self.SSH is None: raise ValueError('Module: kinetic rate for SSH is None, str: '+ self.str)

    def update_value(self,strin=None):
        if strin is not None: self.str=strin
        self.pyformat = kinetic_MCM_to_value(self.str)

    def add_info(self,key,value,check=False):
        if check:
            if not hasattr(self, key):
                print(self.name, key, value)
                raise AttributeErrors('Species add properties: do not have this attr')
        setattr(self, key, value)

