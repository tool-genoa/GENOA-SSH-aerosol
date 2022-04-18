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
import json
import numpy as np

from Parameters import TINYM,NokeepSp,primaryVOCs, \
                       Ttot,DeltaT
from Functions import SecondToDate, compare

# for check inorganic aerosol name
Nsps = ['BiMT','BC','MD','SOAlP','H2O','rganics','NH4','NA','SO4']

def get_time_array(path):
    """return the time array from the input path (only one path)
        MUST HAVE THE FILE: report.txt"""
    return np.arange(0.,432001.,3600.) # for silend mode

    # get time array from report
    with open (path+'/report.txt') as f:
        info=f.read().splitlines()[5].split(' ')
    times = [float(item) for s in info for item in s.split() if isfloat(item)]
    #print(path,times)
    
    # time array
    a_time=np.arange(times[0],times[1]+0.1,times[2])
    nt=len(a_time)
    #print(path,'time array: ',times[0],times[1],times[2],nt)
    if nt-1 != int(times[3]): 
        print('time array length not match: '+str(nt)+' '+str(times[3]),path)
        raise ValueError('check simulation time in report.txt')

    return a_time

def get_conc(paths):
    """return the name and conc. for all input paths (one or more path are accepted)
       MUSt HAVE SUBFOLDERS: /gas/ and /aero/
       return format: a lit of data: [[species name], [gas_conc,aerosol_conc], a_time"""
    # data array
    ns=len(paths)
    dataAll=[]

    # read data
    for i in range(ns):

        # add store 0: name; 1: gas and aero; 2 :time
        data=[[],[]]

        # get time array from report
        data.append(get_time_array(paths[i]))
        nt = len(data[-1])

        # gas
        pa=paths[i]+'/gas/'
        files=os.listdir(pa)

        for j in files: # record all files
            if '.txt' not in j: continue
            sp=j.replace('.txt','')

            # if in gas: create for both gas/aerol conc.
            tmp=np.zeros((2,nt))

            # get gas conc
            with open (pa+'/'+j,'r') as f: info=f.read().splitlines()
            for k in range(nt): 
                tmp[0][k]=float(info[k])
                if tmp[0][k] < TINYM: tmp[0][k] = 0.0

            data[0].append(sp)
            data[1].append(tmp)

        # aero
        pa=paths[i]+'/aero/'
        files=os.listdir(pa)

        for j in files:
            # if not a single species
            if '_1.txt' not in j: continue

            # get precursor name
            sp = j.replace('_1.txt','')[1:]
            if sp in Nsps: continue

            # if precursor recorded
            if sp in data[0]:
                ind = data[0].index(sp)
                # get aero conc
                with open (pa+'/'+j,'r') as f: info=f.read().splitlines()

                for k in range(nt): 
                    data[1][ind][1][k]=float(info[k])
                    if data[1][ind][1][k] < TINYM: data[1][ind][1][k] = 0.0
            # precursor not record
            else:
                tmp=np.zeros((2,nt))

                # get aerosol conc
                with open (pa+'/'+j,'r') as f: info=f.read().splitlines()

                for k in range(nt): 
                    tmp[1][k]=float(info[k])
                    if tmp[1][k] < TINYM: tmp[1][k] = 0.0

                data[0].append(sp)
                data[1].append(tmp)

        dataAll.append(data)

    return dataAll

def ref_conc(paths,species='all',Type='sum',radical='',savname=None,write=False,slist=None):
    """only for gas + aerosol !!! type: sum/ave"""

    if not write and savname is not None: # input paths is a numpy file that saved data
        with open ('{:s}_paths'.format(savname), 'r') as f:
            read_paths = json.loads(f.read())
        if compare(paths,read_paths):
            with open('{:s}_rcConcs.npy'.format(savname), 'rb') as f:
                read_rcConcs = np.load(f)
            with open ('{:s}_rcSps'.format(savname), 'r') as f:
                read_rcSps = json.loads(f.read())
            if species == 'all': # output all species
                return read_rcSps,read_rcConcs
            else: # init
                rcSps, rcLump = [],[]
                for i in species:
                    rcSps.append(i[0])
                    rcLump.append(i[1])

                nrsp=len(rcSps)
                rcConcs=np.zeros(nrsp)

                # get iterated conc of SOA
                n, n0 = [], [] # check species without conc.
                for i in range(nrsp):
                    # check lumped species
                    for j in rcLump[i]:
                        if radical+j in read_rcSps: tmp = radical+j
                        elif j in read_rcSps: tmp = j
                        else: tmp = None
                        if tmp: rcConcs[i] += read_rcConcs[read_rcSps.index(tmp)]
                        else: n0.append(j)

                    # check itself
                    if radical+rcSps[i] in read_rcSps: tmp = radical+rcSps[i]
                    elif rcSps[i] in read_rcSps: tmp = rcSps[i]
                    else: tmp = None
                    if tmp: rcConcs[i] = read_rcConcs[read_rcSps.index(tmp)]
                    else: n.append(rcSps[i])

                # print for check
                if n != []: 
                    print('SSR: not found ref_conc for sps: ',radical,len(n),n,savname)
                    raise NameError
        else: 
            print('SSR ref_conc: paths and read_paths are different.',paths,read_paths,savname)
            raise ValueError
    else: # paths is a list
        if slist is not None:
            MWs = [[],[]]
            # read species name and MWs
            with open (slist,'r') as f: info = f.read().splitlines()[4:]
            for i in info:
                MWs[0].append(i.split(' ')[0]) #name
                MWs[1].append(float(i.split(' ')[-1])) # MWs
                
        datas = get_conc(paths)  #[[species name], [gas_conc,aerosol_conc], a_time]

        # check species
        rcSps, rcLump = [],[]
        if species == 'all': # output all species
            for i in datas:
                for j in i[0]:
                    if j not in rcSps: 
                        rcSps.append(j)
                        rcLump.append([]) 
        else: # output only the requested species
            for i in species:
                rcSps.append(i[0])
                rcLump.append(i[1])

        nrsp=len(rcSps)
        rcConcs=np.zeros(nrsp)

        # get iterated conc of SOA
        for i in range(nrsp):
            for data in datas:
                # itself
                if radical+rcSps[i] in data[0]:
                    rconc = np.sum(data[1][data[0].index(radical+rcSps[i])]) #[0] for gas
                elif rcSps[i] in data[0]:
                    # compute rcConcs TM
                    rconc = np.sum(data[1][data[0].index(rcSps[i])]) #[0] for gas
                else:
                    print('RS: no conc. found in original case for species '+rcSps[i], savname)
                    raise ValueError('check the results of the fake chem.')

                if rconc > 0.0:
                    if slist is not None:
                        nrd = len(radical)
                        if rcSps[i] in MWs[0]: isp = rcSps[i]
                        elif rcSps[i][nrd:] in MWs[0]: isp = rcSps[i][nrd:]
                        else:
                            print('SSR: not found: ',rcSps[i],nrd,radical,'in given list: ',slist)
                            raise NameError('check the results of the fake chem.')
                        rcConcs[i] += rconc / MWs[1][MWs[0].index(isp)]
                    else: rcConcs[i] += rconc
                
                # sum conc of lumped species
                for j in rcLump[i]:
                    if radical+j in data[0]:
                        rconc = np.sum(data[1][data[0].index(radical+j)]) #[0] for gas 
                    elif j in data[0]:
                        rconc = np.sum(data[1][data[0].index(j)]) #[1] 
                    else:
                        rconc = 0.0
                        print('RS: no conc. found in original case for lumped species '+j)
                        
                    if rconc > 0.0:
                        if slist is not None:
                            nrd = len(radical)
                            if j in MWs[0]: isp = j
                            elif j[nrd:] in MWs[0]: isp = j[nrd:]
                            else:
                                print('SSR: not found in given list',j, nrd, radical, slist)
                                raise NameError
                            rcConcs[i] += rconc / MWs[1][MWs[0].index(isp)]
                        else: rcConcs[i] += rconc
                                        
        if Type == 'ave': 
            # time steps
            n = 0
            for data in datas: n += len(data[2])
            for i in range(nrsp): rcConcs[i]/=float(n)

        # save
        if savname is not None:
            with open ('{:s}_rcSps'.format(savname), 'w') as f: json.dump(rcSps,f)
            with open ('{:s}_paths'.format(savname), 'w') as f: json.dump(paths,f)
            np.save('{:s}_rcConcs'.format(savname),rcConcs)
            print('SSR: save ref_conc by the name:', savname)
    return rcSps,rcConcs

def get_organic_conc(paths, Type = 'aero'):
    """return only total SOA conc for auto cases"""
    # data array
    ns=len(paths)
    data=[np.arange(0.,Ttot[0]+0.1,DeltaT[0]),[]]

    # read data
    for i in range(ns):

        # get time array from report
        nt = len(data[0])
        data[1].append(np.zeros(nt))

        # get total aero conc
        if Type == 'aero': 
            with open (paths[i]+'/aero/Organics_1.txt','r') as f: info=f.read().splitlines()
            for k in range(nt): 
                data[1][i][k] = float(info[k])
                # trim data
                if data[1][i][k] < TINYM: data[1][i][k] = 0.0

        elif Type == 'yield': # output yield
            with open (paths[i]+'/aero/Organics_1.txt','r') as f: aero=f.read().splitlines()
            gas = []
            for ivoc in primaryVOCs:
                with open (paths[i]+'/gas/'+ivoc+'.txt','r') as f: gas.append(f.read().splitlines())
            aero0 = float(aero[0])

            for k in range(nt):
                tmp = 0.0
                for i in gas: tmp += float(gas[0]) - float(gas[k])
                data[1][i][k] = (float(aero[k]) - aero0)/(tmp + TINYM) * 100.
                # trim data
                if data[1][i][k] < TINYM: data[1][i][k] = 0.0

        elif Type == 'plain':
            with open (paths[i],'r') as f: info=f.read().splitlines()
            for k in range(nt):
                data[1][i][k] = float(info[k])
                # trim data
                if data[1][i][k] < TINYM: data[1][i][k] = 0.0

        else: 
            print('SSR: get_organic_concs type not find ', Type)
            raise TypeError

    return data
