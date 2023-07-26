# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  Functions.py is a module that contains commonly used functions intended
#
#    for reuse across multiple scripts within GENOA.
#
# ================================================================================

import os
import re
import csv
import ast
import math
import shutil
import numpy as np

from datetime import datetime, date
from itertools import permutations, combinations

from Parameters import lats, lons

try:
    import angzen
except:
    os.system('python3 -m numpy.f2py -c angzen.edf.f -m angzen')
    import angzen

# operations on numbers
def isint(value):
    """check if a string is an integer."""
    try:
        int(value)
        return True
    except ValueError:
        return False

def isfloat(value):
    """check if a string is a float"""
    try:
        float(value)
        return True
    except ValueError:
        return False

def isFitScale(t1,t2,val):
    t1 = float(t1)
    t2 = float(t2)
    val = float(val)
    if val < 1.0: 
        raise ValueError('FC: isFitScale val should >= 1.0: '+str(val))
    if t1 == 0.0 or t2 == 0.0:
        if t1 == t2 : return True
        else: return False
    if t1/t2 <= val and t1/t2 >= (1./val) : return True
    else: return False

def round_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier

def round_down(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n * multiplier) / multiplier

# operations on list elements
def multiplyList(myList) :       
    """ return the multiply result of all elements in the given list""" 
    result = 1
    for x in myList: 
         result = result * x  
    return result

def list_all_fit(lst1,lst2,val):
    """return True if each element in lists fit the input condition val
       ex. lst1[i] val lst2[i] 
       lst2 can be np.array, list, and string"""

    if isinstance(lst2,str) :
        for i in range(len(lst1)):
            if not eval('lst1[i] '+val+' lst2'): return False
        return True
    else:
        if len(lst1) != len(lst2):
            print(len(lst1), len(lst2))
            raise ValueError('FC list_all_fit: number of list1 and list2 not match')
        for i in range(len(lst1)):
            if not eval('lst1[i] '+val+' lst2[i]'): return False
        return True

def takeSecond(elem):
    """return the second element"""
    return elem[1]

def count_cut(num):
    """return a list of [1,num], [2,num], .... [num, num] if num >= 2"""
    if num >=2 : return [[i+1,num] for i in range(num)]
    else: return [None]

def combinations_number(inval,num0=0):
    """return the number of possibility of combinations."""
    
    if not isinstance(inval,list): 
        if isinstance(inval,int): inval=list(range(inval))
        else: raise ValueError('input should be integer or list.',inval)
    n = 0 # init num
    for i in range(num0,len(inval)):
        n += len(list(combinations(inval,i+1)))
    return n

def combinations_index(inval,num0=0):
    """return the index of all combination possibilities of a list with the input number"""

    if not isinstance(inval,list): 
        if isinstance(inval,int): inval=list(range(inval))
        else: raise ValueError('input must be integer or list.',inval)
    val = []
    for i in range(num0,len(inval)):
        for j in list(combinations(inval,i+1)):
            val.append([k for k in j])
    return val

def compare(s, t):
    """ return True if two lists have same elements (can be unsorted) """
    return sorted(s) == sorted(t)

def duplicates(lst, item):
    """return a list of the index of duplicates"""
    return [i for i, x in enumerate(lst) if x == item]

def isContain(mainlist,sublist):
    """return True if sublist in mainlist"""
    return set(sublist).issubset(set(mainlist))

def addNoDuplicate(vals,myList):
    for i in vals:
        if i not in myList: myList.append(i)
  
def trim_items(sps,items):
    return list(set(sps)-set(items))

def array_to_ratio(a_list):
    """output the ratios of input array"""
    n = len(a_list)
    if n == 1: return [1.0]
    elif n == 0: return []
    else:
        tmp = sum(a_list)
        if tmp == 0.0: tmp1 = [1.0/n]*n
        else:
            tmp1 = [i/tmp for i in a_list]
        if not sum(tmp1) == 1.0: tmp1[-1] = 1.0 - sum(tmp1[:-1])
        return tmp1

# operations on time
def SecondToDate(val, Type = '%d-%Hh'): #,Type="%B %d, %I:%M"):
    #return datetime.fromtimestamp(val).strftime(Type)
    return str(int(val / 3600.))

# operations on files/folders
def is_same_file(file1, file2):
    """return the difference"""
    if not os.path.exists(file1) or not os.path.exists(file2): return False
    with open (file1,'r') as f1:
        with open (file2,'r') as f2:
           if set(f1) == set(f2): return True 
           else: return False

def replace_in_file(item1,item2,filein,fileout):
    """replace item in filein and output as fileout and return the times that replaced"""
    with open (filein,'r') as fin:
        info = fin.read()
        with open (fileout,'w') as fout:
            fout.write(info.replace(item1,item2))
    return info.count(item1)

def create_folder(path, del_exist = False):
    """generate the new folder if it does not exist.
       clean up all files/sub-directories if it exists."""
    if os.path.exists(path):
        if del_exist: shutil.rmtree(path, ignore_errors=True)
        else: return False # no need to create
    os.makedirs(path, exist_ok=True)
    return True # create once

def move_results(path_old, path_new, items = []):
    """remove the results from the old path to the new path"""
    #path_old = '{:s}/{:s}'.format(pathSSH,ResultFolder)
    #path_new = '{:s}/{:s}'.format(pathResult,ResultFolder)
    # check if exists
    if os.path.exists(path_new):
        # file already contains
        exist_locs = os.listdir(path_new)
    else: # generate new file
        os.makedirs(path_new, exist_ok=True)
        exist_locs = []
    # file need to transfer
    tar_locs = os.listdir(path_old)
    # update location folder
    if items == []: 
        results = tar_locs
    else: # remove repeated paths
        results = list(set([i for i in items if i in tar_locs]))
    if results == []:
        print(items, tar_locs, path_old, path_new)
        print('Warning! No find targeted files to be moved.')
    for lc in results:
        # clean old folder if exist
        if lc in exist_locs:
            shutil.rmtree('{:s}/{:s}'.format(path_new,lc), ignore_errors=True)
        # mv new result
        shutil.move('{:s}/{:s}'.format(path_old,lc),
                    '{:s}/{:s}'.format(path_new,lc))
        tar_locs.remove(lc)

    # clean old empty folder
    if os.listdir(path_old) == []: os.rmdir(path_old)

# operations related to atmospheric calculations
def ug_to_ppb(val,MW,TEMP = 298.):
    """transfer the input mass in unit ug/m3 to in unit ppb. The default temperature is 298 K"""
    return val * TEMP / (12.187 * MW)

def ppb_to_ug(val,MW,TEMP = 298.):
    """transfer the input mass in unit ppb to in unit ug/m3. The default temperature is 298 K"""
    return val * 12.187 * MW / TEMP

# PHOTOLYSIS PARAMETERS
def cosx_max(locs):
    """out put the average of the maximum ZSA in the giveb locations"""
    locations = get_locations(locs)
    print('computing ZSA from ',len(locs),' locations...')
    # compute cosx
    cosx = 0.0
    for i in range(len(locs)):
        cosx_i = 0.0
        m  = locs[i][2] # get month
        if m == 11: nm = (date(2016,1,1) - date(2015,m+1,1)).days
        else: nm = (date(2015,m+2,1) - date(2015,m+1,1)).days
        for d in range(nm): # days in the month
            for t in range(24):
                time = (date(2015,m+1,1+d) - date(2015,1,1)).total_seconds() + 3600 * t
                cosx_i = max(cosx_i, angzen.ssh_muzero(time,locations[i][1],locations[i][0]))
        cosx = max(cosx, cosx_i)
    if cosx == 0.0: 
        raise ValueError('check cosx_max: cosx should not be zero')
    else:
        secx=1./cosx
        print('max, COSX,SECX: ',cosx, secx)
        return cosx, secx

def get_locations(locs):
    # get locations
    out = [] # [lat, lon]
    for i in locs: out.append([round(lats[i[0]],2),round(lons[i[1]],1)])
    return out

def get_info(pathorrea, IDchemorsp, Type = 'read'):
    """Output the info of given scheme"""
    if Type == 'read':
        # get number of reactions/species
        tmp = pathorrea+IDchemorsp+'/'+IDchemorsp
        with open (tmp+'.species') as f:
            #print('Number of gas species: ', int(f.read().splitlines()[2].split(' ')[0])-17)
            gas = int(f.read().splitlines()[2].split(' ')[0])-18
        with open (tmp+'.aer.vec') as f:
            #print('Number of aerosol species: ', len(f.read().splitlines())-11)
            aer = len(f.read().splitlines())-11
        with open (tmp+'.reactions') as f:
            rea = f.read().count('KINETIC') - f.read().count('%KINEITC')
        if 0:
            with open (tmp+'.RO2') as f:
                print('Number of RO2: ', len(f.read().splitlines())-1)
            with open (tmp+'.photolysis') as f:
                print('Number of Photolysis: ', len(f.read().splitlines())-1)
    elif Type == 'com':
        rea, gas, aer, ro2, j = 0, 0, 0, 0, 0
        for i in IDchemorsp:
            if i.status and i.organic: 
                gas += 1
                if i.condensed: aer += 1
                if i.RO2: ro2 += 1

        for i in pathorrea:
            if i.status: rea += 1
            if 'J<' in i.rate.str: j += 1
        return rea, gas, aer, ro2, j
    else:
        rea, gas, aer = 0, 0, 0
        for i in IDchemorsp:
            if i.status and i.organic: 
                gas += 1
                if i.condensed: aer += 1
        for i in pathorrea:
            if i.status: rea += 1

    return rea, gas, aer

def get_negative_str(val):
    "get the negative string"
    if val[0] == '-': return val[1:]
    else: return '-' + val
    
def convert_number_format_in_string(string):
    """reformat number in a string in format aEb to a$^{b}$ used in latex"""
    
    pattern = r"(\d+(?:\.\d+)?)E([+-]?\d+)"
    matches = re.findall(pattern, string)
    converted_string = string
    for match in matches:
        a = float(match[0])
        a_str = "{:.{}f}".format(a, len(match[0].split(".")[1].rstrip("0")))
        a_str = a_str.rstrip("0") if "." in a_str else a_str
        a = float(a_str)
        b = int(match[1])
        converted_number = "{:.1f}$\\times 10^{:}$".format(a,"{"+str(b)+"}")
        converted_string = converted_string.replace(match[0] + 'E' + match[1], converted_number)
    return converted_string

def reformat_number_in_string(number_str):
    """reformat number in a string separated by &
       aEb to a$^{b}$ used in latex"""
    
    numbers = number_str.split('&')

    for i, num in enumerate(numbers):
        try:
            num = float(num)
            if 0.01 < num < 1000:
                a_str = "{:.{}f}".format(num, len(str(num).split(".")[1].rstrip("0")))
                numbers[i] = a_str.rstrip("0") if "." in a_str else a_str
            else:
                numbers[i] = '{:.2e}'.format(num).replace('e', '$\\times 10^{').replace('+', '') + '}$'
        except ValueError:
            pass

    return '&'.join(numbers)

# functions to read input parameters
def read_parameters_from_csv(filename, delimiter=','):

    if not os.path.exists(filename):
        raise FileNotFoundError(f"File '{filename}' not found. Please check the file path.")

    input_parameters = {}
    temp_dict = {}
    param_lengths = None

    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=delimiter)
        next(csvreader)  # Assuming the first row is the header row
        for row in csvreader:
            if not row[0].strip().startswith("#"):  # Skip lines starting with #
                param_name = row[0] # first column (parameter names)
                if not param_name: continue
                param_values = [value.strip() for value in row[1:] if value != '']
                
                # check name
                if param_name in list(input_parameters.keys()):
                    raise NameError('Parameter given twice in the csv table',param_name,filename)

                # Check if the lengths of input lists are consistent
                if param_lengths is None:
                    param_lengths = len(param_values)
                    
                elif len(param_values) != param_lengths:
                    raise ValueError(f"The length of input lists for parameter '{param_name}' is inconsistent.", param_lengths, len(param_values))
                
                temp_dict[param_name] = param_values

        for param_name, param_values in temp_dict.items():
            processed_values = []
            for value_str in param_values:
                try:
                    value = ast.literal_eval(value_str)
                    if isinstance(value, int):
                        processed_value = int(value)
                    elif isinstance(value, float):
                        processed_value = float(value)
                    elif isinstance(value, str):
                        if value in temp_dict: 
                            val,n = value,1
                        else:
                            if '/' in value: spr = '/'
                            elif '*' in value: spr = '*'

                            val = value.split(spr)[0].replace(' ','')
                            n = float(value.split(spr)[1])

                        if val in temp_dict:
                            # If the value is a string reference to another parameter, calculate the result
                            ref_param_values = temp_dict[value][len(processed_values)]
                            processed_value = float(ref_param_values) / float(value.split('/')[1])
                    else:
                        processed_value = value_str

                except (SyntaxError, ValueError):
                    processed_value = value_str
                    
                processed_values.append(processed_value)

            input_parameters[param_name] = processed_values

    return input_parameters

