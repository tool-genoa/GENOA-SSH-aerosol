# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  Module.py contains three classes: Species, Reaction, and Kinetic. 
#
#  These classes are utilized to handle chemical species and reactions in GENOA.
#
# ================================================================================

import json
from copy import deepcopy

from Functions import duplicates,trim_items, \
                      isfloat
from Parameters import NokeepSp
from KineticMCMtoSSH import MCMtoSSH_kinetic_rate
from KineticMCMtoPython import kinetic_MCM_to_value

### Need to install openbabel and UManSysProp
# decomment to use some functions in this module
#import MolProperty as mp

sav_species=['name','nameUsed','status',
             'lump','jump',
             'organic','RO2','Radical','condensed',
             'SMILES','InChI','mass','formula','source','groupID',
             'psat_torr','psat_atm','dHvap_KJ','gamma','henry','Kp',
             'temp','DU','ratios','functionalGroups',
             'SOAPStructure','non_volatile','generation']

class Species:
    """record properties of a given species"""
    def __init__(self,name=None):

        self.status=True
         
        # species name
        self.name=name
        # name changed dur to lumping
        self.nameUsed=None

        # properties
        self.organic=False # if it is organic
        self.RO2=False     # if it is a peroxy radical (consider in the RO2 pool)
        self.Radical=False # if it is a radical
        self.condensed=False # if it is condensable

        # smiles, used to compute Psat
        self.SMILES='-'
        # read from MCM, currently not used in reduction
        self.InChI=None
        self.mass=0. # MWs
        self.formula='' # formula in the format: C[i]H[i]N[i]O[i], i is number

        # reduction related
        self.lump=[] # lumped species: output "m"+name in .viz file if it is not []
        self.jump=[] # jumped species
        self.remove=[] # removed species

        # for SSH output - see SSH-aerosol aerosol list for more info
        self.type=4    # Type No. in SSH-aerosol
        self.groupID=2 # group ID: organic
        self.coll_fac=687.
        self.mole_diam=8.39
        self.surf_tens=30.e-03
        self.accomod=1.0
        self.mass_dens=1.30e-06
        self.non_volatile=0 # if it is a non-volatile compound

        self.psat_torr=0.0 # Saturation vapor pressure at Tref (torr)
        self.psat_atm=0.0  # Saturation vapor pressure at Tref (atm)
        self.dHvap_KJ=0.0  # Enthalpy of vaporization (kJ/mol)
        
        # those values may not be updated after reduction: recompute in SSH-aerosol
        self.gamma=None    # Activity coefficient at infinite dilution in water
        self.henry=0.0     # Henry's Law constant (M/atm, mol/atm*l)
        self.Kp=None       # partitioning coefficient

        self.temp= 298.    # Temperature of reference (K)

        self.functionalGroups={} # No. of some functional groups (used to define lumpable species)
        self.ratios=None # atomic ratios: OM/OC, N/C, H/C, O/C
        self.DU=None     # degree of unsaturation
        self.pymol=None  # pymol format: computed by openbabel from smiles

        self.generation = -1 # number of generation: -1 not computed

        self.SOAPStructure = None # No. of non-zero function groups in SSH-aerosol vector format

        self.source = [] # primary VOCs

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


    def fromText(self,text):
        for i in text:
            tmp = i.split('\t')
            if len(tmp) == 2: 
                tmp[1] = json.loads(tmp[1])
                if isinstance(tmp[1],str): tmp[1] = str(tmp[1])
                self.add_info(tmp[0],tmp[1],True)

    def to_species_list_aer(self, Type = 'smiles', firstRun = False):
        if self.organic:
            spr='\t'
            line='P{:s}'.format(self.name)

            if Type == 'smiles': tmp = self.SMILES
            elif Type == 'vectors': tmp = self.VectorinText()
            else:
                raise NameError('MD: unknown type: '+Type)

            content=[self.type,self.groupID, self.mass,self.name,self.coll_fac,
                    self.mole_diam,self.surf_tens,self.accomod,self.mass_dens,
                    self.non_volatile,'BOTH',tmp,self.psat_torr,self.dHvap_KJ]
            # remove precursor if it is the 1st run
            if firstRun: content[3] = '--'

            style=['{:d}','{:d}','{:6.2f}','{:s}','{:5.2E}',
                    '{:5.2E}','{:5.2E}','{:5.2E}','{:5.2E}',
                    '{:d}','{:s}','{:s}','{:5.2E}','{:5.2E}']
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
        if self.psat_atm and self.psat_atm <= AeroCriteria['Psat_NVOC'] : 
            print('species ',self.name,' is non_volatile, with P_sat = ',self.psat_atm,' < ',AeroCriteria['Psat_NVOC'])
            self.non_volatile = 1

        # Hvap
        self.dHvap_KJ = mp.enthalpy_vaporization(self.pymol,AeroCriteria['vpType'])/1000. # unit in KJ

        # functional group
        self.functionalGroups = mp.functional_groups(self.pymol)

        # is RO2? Radicals? Condensed species? 
        self.is_type(AeroCriteria['Psat_aero'])

        # SOAPstructure ?
        self.SOAPStructure = mp.toSOAPStructure(self.pymol)

        # atomic ratios and degree of unsaturation
        self.ratios,self.DU=mp.organic_ratios_and_du(self.pymol)

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
        self.is_type(AeroCriteria['Psat_aero'])

        # SOAPstructure ?
        self.SOAPStructure = mp.toSOAPStructure(self.pymol)

        # activity coefficient
        #self.gamma=mp.activity_coefficient_inf(self.pymol,self.temp)

        # henry's law constant (M/atm)
        # H = 1000*760/(18.0*GAMMAinf*Psat)
        self.henry = 1000/(18*self.gamma*self.psat_atm)

        # partitioning coefficient
        self.Kp=mp.partitioning_coefficient(self,self.temp)

        # atomic ratios and degree of unsaturation
        self.ratios,self.DU=mp.organic_ratios_and_du(self.pymol)

    def update_Psat(self,AeroCriteria):

        # check
        if not self.status: return None
        if self.pymol is None:
            if self.SMILES != '-':
                self.pymol = mp.get_pybelmol(self.SMILES)
            else:
                return None

        # psat
        self.psat_atm = mp.saturation_vapor_pressure(self.pymol,AeroCriteria['vpType'],self.temp)

        self.psat_torr=self.psat_atm*760. # unit

        # Hvap
        self.dHvap_KJ = mp.enthalpy_vaporization(self.pymol,AeroCriteria['vpType'])/1000. # unit in KJ

        # henry's law constant (M/atm)
        # H = 1000*760/(18.0*GAMMAinf*Psat)
        if self.gamma: self.henry = 1000/(18*self.gamma*self.psat_atm)

    def update_mass(self):

        if self.pymol is not None: self.mass = self.pymol.molwt
        elif self.SMILES != '-': mol=self.SMILES
        else: mol=self.name

        # mass
        self.mass=mp.molar_mass(mol)

    def is_type(self,Psat):
        
        #if self.pymol.spin>1 or '[O+]' in self.SMILES or '[O]' in self.SMILES:
        if '[O+]' in self.SMILES or '[O]' in self.SMILES:
            #print('spin==',self.pymol.spin,self.name,self.SMILES)
            self.Radical=True
            if self.functionalGroups == {} : self.functionalGroups = mp.functional_groups(self.pymol)
            #if 'RO2' in self.functionalGroups.keys(): 
            if self.name[-2:] == 'O2': self.RO2=True
            if self.name[-2:] not in ['OO','O2','O3']:
                if self.name[-1:] != 'O':
                    print('MD: radicals not in O, OO, O2, O3: ',self.name,self.pymol.spin)

        if not self.Radical:
            if Psat == 0.0 or self.psat_atm <= Psat : 
                #print(self.name,self.psat_atm,Psat,self.psat_atm <= Psat)
                self.condensed=True

    def toSOAP(self,info):
        if not self.status or self.SOAPStructure is None: return False

        iStr = self.SOAPStructure
        tmp = ("  double group_tmp_"+self.name+" [] = {")
        tmp0 = [0.0] * 60
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
                    if self.psat_atm <= AeroCriteria['Psat_NVOC']: 
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
            for i in range(60):
                if i in self.SOAPStructure[0]: tmp += '&{:5.2E}'.format(self.SOAPStructure[1][self.SOAPStructure[0].index(i)])
                else: tmp += '&{:5.2E}'.format(0.0)
        #else: raise AttributeError('species '+self.name+' has SOAPstructure None')
        if tmp != "": return tmp
        else: return "-"

class Reaction:
    """record a chemical reaction"""
    def __init__(self):
    
        self.status=True
            
        self.reactants=[] # list of reactants
        self.ratiosRC=[]  # list of reactant ratios (default 1)
        self.products=[]  # list of products
        self.ratiosPD=[]  # list of profuct ratios (default 1)

        self.rate=Kinetic() # kinetic contant
        self.type=None      # based on kinetic, used to merge reaucitons

        self.species=None # related species

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
        """ Check the reaction type"""
        
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
            # check length
            if len(sps[n]) != len(rts[n]):
                raise ValueError('len of sps and rts are not the same.',n,sps[n],rts[n])
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

        # change pd if clean version
        if Type == 'clean': # clean version
            # add products and ratio
            rcn=deepcopy(self.products)
            rts=deepcopy(self.ratiosPD)
            
            # remove products
            tmp = list(set(self.products)&set(NokeepSp))
            for i in tmp:
                j = rcn.index(i) # index
                rcn.pop(j)
                rts.pop(j)
                
            # add reactants
            tmp = list(set(self.reactants)&set(NokeepSp))
            for i in tmp:
                rcn.append(i)
                rts.append(1.0)
                
        else: # read but not changed
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

        elif Type in ['all','clean']:
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
                if 'J<' in line: self.rate.Photolysis = True
            # kinetic rate
            elif 'KINETIC' in line:
                # kinetic rate in SSH format
                self.rate.SSH = line
                if 'PHOTOLYSIS' in line: self.rate.Photolysis = True

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

    def check_conservation(self):
        # get reactant list
        rea,rea1 = [],[]
        for i in self.reactants:
            if i in NokeepSp:
                rea.append(i)
                
        # check if inorganics are conserved
        i = 0
        while i <= (len(self.products) - 1):
            s = self.products[i]
            if s in NokeepSp:
                if s in rea:
                    self.ratiosPD[i] = 1.0
                    i += 1
                    rea1.append(s)
                else:
                    print('rcn conservation: del ',self.products[i],' from ',self.toSSH())
                    del self.products[i]
                    del self.ratiosPD[i]
            else: 
                i += 1
        if rea != []:
            for i in rea:
                if i not in rea1:
                    self.products.append(i)
                    self.ratiosPD.append(1.0)
        
class Kinetic:
    """store kinetic in MCM string and ssh format"""

    def __init__(self,strin=None):

        self.str=strin     # input kinetic as a string
        self.pyformat=None # in a format can be executed by python
        self.Photolysis=False # if it is a photolysis reaction
        self.SSH=None         # in SSH-aerosol format

        if strin is not None: self.update()

    def update(self,strin=None):
        """Update kineitc"""
        if strin is not None: self.str=strin
        self.SSH = MCMtoSSH_kinetic_rate(self.str)
        #if 'J<' in self.str or 'PHOTOLYSIS' in self.str: self.Photolysis=True
        #else: self.Photolysis=False

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

