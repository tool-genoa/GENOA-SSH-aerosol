# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  MolProperty.py runs prereduction depending on the given reduction
#   parameters and options. 
#
#  It utilizes two external libraries: OpenBabel and UManSysProp
#
# ================================================================================

import math,sys

# OpenBabel (to be installed) - get smiles structures
from openbabel import pybel

from Parameters import keyels,element_mass
from Functions import isfloat

# UManSysProp (to be downloaded)
sys.path.append('UManSysProp_public-master')
from umansysprop import groups
from umansysprop import boiling_points
from umansysprop import vapour_pressures
from umansysprop import activity_coefficient_models

key_pybel={}
for key,val in keyels.items():key_pybel[key]=pybel.Smarts(val)

# store VP and BP method names used in umansysprop
mthVP=['nannoolal','myrdal_and_yalkowsky']
mthBP=['nannoolal','stein_and_brown','joback_and_reid']

def get_pybelmol(mol):
    """get a molecule to pybel format"""
    if isinstance(mol,pybel.Molecule): return mol
    # if input a SMILES string
    elif isinstance(mol,str):
        try:
            mol=pybel.readstring('smi',mol)
            return mol
        except IOError:
            print('get_pybelmol: not a smi string: '+mol)    
            return False
    else:
        print(mol) 
        sys.exit('get_pybelmol: not recognize data type, not mol and not string.')

def molar_mass(mol):
    """return molar mass: compute by element_mass provided in Parameters.py"""

    line=""
    for i in mol:
        if i in element_mass.keys(): line+="+element_mass['{:s}']".format(i)
        elif isfloat(i): 
            if isfloat(line[-1]): line+=i
            else: line+="*"+i
        else: sys.exit('estimate mass not recognize: '+i+' in '+mol)

    return eval(line[1:])

def functional_groups(mol):
    """return a dictionary of the key functional groups read by pybel
        the key functional group can be modified in Parameter.py"""

    pymol=get_pybelmol(mol)

    functionalGroups={}
    #for key,val in smarts_pybel.items(): 
    for key,val in key_pybel.items(): 
        tmp = val.findall(pymol)
        if tmp !=[]: functionalGroups[key]=len(tmp)*1.0
        else: functionalGroups[key]=0.
    return functionalGroups

def organic_ratios_and_du(mol):
    """get the number of key chemical elements in the given molecule,
        return OM/OC mass ratio, H/C, O/C, N/C atomic ratios and the degree of unsaturation"""
    # For a compound with formula CaHbNcOdXe where X is F, Cl, Br or I, the degree of unsaturation is given by:
    # degree of unsaturation = 1/2 (2 + 2C + N - H - X)
    pymol=get_pybelmol(mol)
    tmp={}
    for key,val in key_pybel.items():
        #print(key,val.findall(pymol),len(val.findall(pymol)))
        tmp[key] = float(len(val.findall(pymol)))
    return {'OM/OC':round(pymol.molwt/(tmp['C']*12.),3),
            'H/C':round(tmp['H']*1./tmp['C'],3),
            'O/C':round(tmp['O']*1./tmp['C'],3),
            'N/C':round(tmp['N']*1./tmp['C'],3)},int(0.5*(2+2*tmp['C']-tmp['H']+tmp['N']))


def toSOAPStructure(mol):
    """get SSH-aerosol strucutre of input compounds.
       input: one SMILES strucutre -> pymol
       output: a dict of SSH-aerosol functional groups {index: number,...}"""
    m = groups.to_SSH(mol)
    # transfer to SOAP format
    mkeys = list(sorted(m.keys()))
    return [mkeys,[float(m[i]) for i in mkeys]]

def saturation_vapor_pressure(mol,Type='evap',temperature=None):
    """return (Psat,T) Psat at the certain temperature T in unit: atm, K"""

    pymol=get_pybelmol(mol)

    if temperature is None: temperature=298.

    # Boiling points [(K)] bp
    # Vapour pressures [log10 (atm) at a specific temperature] vp

    if Type == 'evap': vp=vapour_pressures.evaporation(pymol, temperature)
    elif Type == 'evap2': vp=vapour_pressures.evaporation2(pymol, temperature)
    elif Type =='simpol': vp=vapour_pressures.simpol(pymol, temperature)
    elif 'VP' in Type and 'BP' in Type:
            # example VP0BP0                
            bp=eval('boiling_points.'+mthBP[int(Type[5])]+'(pymol)')
            vp=eval('vapour_pressures.'+mthVP[int(Type[2])]+'(pymol, temperature, bp)')
    else:
        vp=[]
        # obtain vp
        for i in ['evaporation', 'evaporation2', 'simpol']:
            #print mol, i, eval('vapour_pressures.'+i+'(pymol, temperature)')
            vp.append(eval('vapour_pressures.'+i+'(pymol, temperature)'))

        for i in mthBP:
            bp=eval('boiling_points.'+i+'(pymol)')
            for j in mthVP:
                #print(i,j,eval('vapour_pressures.'+j+'(pymol, temperature, bp)'))
                vp.append(eval('vapour_pressures.'+j+'(pymol, temperature, bp)'))

        # compute
        if Type == 'ave':
            num=0
            for i in vp:num+=10**i
            return num/(len(vp)*1.0)

        elif Type == 'ave_log10':
            return 10**(sum(vp)/(len(vp)*1.0))

        elif Type == 'ave_log10_3':
            vp.sort()
            return 10**(sum(vp[0:3])/3.)

        elif Type == 'ave_log10_3_2':
            vp.sort()
            return 10**(sum(vp[0:3])/3.)/2.

        elif Type == 'min' or Type == 'max':
            return 10**eval(Type+'(vp)')

        else: sys.exit('MP: not recognize saturation vapor pressure type') 

    return 10**vp


def enthalpy_vaporization(mol,Type='evap',temp1=298.0,temp2=308.0):
    """return Hvap in unit J.mol-1"""
    if temp1 != temp2:
        psat1=saturation_vapor_pressure(mol,temperature=temp1,Type=Type)
        psat2=saturation_vapor_pressure(mol,temperature=temp2,Type=Type)

        #  Clausius-Clapeyron equation
        # print(psat1,temp1,psat2,temp2,1/temp1-1/temp2)
        return abs(math.log(psat2/psat1)*8.31446/(1.0/temp1-1.0/temp2))

    else: sys.exit('Hvap: temp1 = temp2, '+str(temp1))

# for gamma
water = get_pybelmol('O')

def activity_coefficient_inf(compound,temperature):
    """ This routine compute the activity coefficients at infinite dilution GAMMAinf
         given by UNIFAC for each species and compute the Henry's law constant from
         the saturation vapour pressure
         H = 1000*760/(18.0*GAMMAinf*Psat)"""

    item=get_pybelmol(compound)
    organic = {item:1,water:1e50}#1e10

    dic,m=activity_coefficient_models.calculate_activities_org(organic, temperature)

    dic = {}
    for c in m: dic[c.compound]=c.activity_coefficient

    # check water
    if dic[water] != 1.0 : sys.exit('water is not infinite: {:f}, gamma_inf: {:f}'.format(dic[water],dic[item]))
    else: return dic[item]
        

def activity_coefficient_org(compounds,concs,temperature):
    """input: species name and pymol in a list"""
    organics = {}
    # add conc.
    for i in range(len(compounds)): organics[compounds[i].pymol]=concs[i]
    # compute
    dic,m=activity_coefficient_models.calculate_activities_org(organics, temperature)
    # out infp
    out=[]
    for i in compounds:
        tag=0
        for c in m:
            if i.pymol == c.compound: 
                out.append(c.activity_coefficient)
                #print(i.name,c.activity_coefficient)
                tag=1
                break
        if not tag: sys.exit('activity_coefficient_org: not find gamma for species '+i.name)
    return out

def partitioning_coefficient(compund,temperature,Mavg=200):
    Kp = 8.314*temperature/(compund.gamma*compund.psat_atm*Mavg)
    #print 'Kp',Kp,'psat_atm',compund.psat_atm
    return Kp
    #return 8.314*temperature/(self.gamma*self.psat_atm*Mavg)

def dHvap_simpol(SOAPStructure):
    simpol1 = {
                'bo':-9.0677E+02,
                'C':	-2.3229E+02,
                'C=C':	5.9336E+01,
                'OH':	-8.7901E+02,
                'aldehyde':	-5.2267E+02,
                'ketone':	1.9917E+01,
                'COOH':	-1.1963E+03,
                'nitrate':	-7.8253E+02,
                'peroxide':	4.4567E+02,
                'hydroperoxide':	-7.9762E+02,
                'aromatic ring':	-1.3635E+02,
                'ether':	-2.2814E+02,
                'phenol':	-4.2981E+02,
                'nitrophenol':	2.8685E+02
             }
    simpol0 = {
                'C':[0,1,2,3],
                'C=C':[16,17,18,19,20],
                'OH':[26],
                'aldehyde':	[31],
                'ketone':[29,30],
                'COOH':	[37],
                'nitrate':	[39,40,41],
                'peroxide':	[45,46,47,48,49,50,51,52,53],
                'hydroperoxide':[42,43,44],
                'aromatic ring':[21,22],
                'ether':[34,35,36],
                'phenol':[28],
                'nitrophenol':[38]
             }
    ksim = list(simpol0.keys())
    dH = 0.0
    for i,j in enumerate(SOAPStructure[0]):
        tag = 0
        for k in ksim:
            if j in simpol0[k]:
                dH +=  SOAPStructure[1][i] * simpol1[k]
                tag = 1
                break
        if tag == 0: print('!!!!not find strucutre: ', j)
    if dH != 0.0: return -1 * (dH + simpol1['bo']) * (2.303*8.314)/1000
    else: return 0.0

if __name__ == '__main__':

    mol = [
           'O=CCCC(=C)C(=O)CC(C)(C)O[N+](=O)[O-]', #
           'O=CCC(C)(C)C(O)CCC(=O)C',
           'OCCC(=C)C(=O)CC(C)(C)C(O)CC=O',
           'CC(=O)CCC1C(CC1(C)C)C(=O)OO[N+](=O)[O-]', #C1011PAN
           'OOc1c(C)cccc1O'
           ]
    for i in mol:
        print(i,toSOAPStructure(get_pybelmol(i)))

