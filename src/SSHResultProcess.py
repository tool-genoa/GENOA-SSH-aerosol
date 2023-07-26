# ================================================================================
#
#   GENOA v2.0: the GENerator of reduced Organic Aerosol mechanism
#
#    Copyright (C) 2023 CEREA (ENPC) - INERIS.
#    GENOA is distributed under GPL v3.
#
# ================================================================================
#
#  SSHResultProcess.py provides functions that analyze the outputs of the box 
#
#   model SSH-aerosol.
#
# ================================================================================

import os
import json
import numpy as np

from Parameters import TINYM,NokeepSp,primaryVOCs, \
                       Ttot,DeltaT
from Functions import SecondToDate, compare

### For plotting - need to decomment to use functions in this module 
#import matplotlib.pyplot as plt
#import matplotlib.pylab as pylab
#params = {   'legend.fontsize': 'x-large',
#             'axes.labelsize': 'x-large',
#             'axes.titlesize': 'x-large',
#             'xtick.labelsize':'x-large',
#             'ytick.labelsize':'x-large'}
#pylab.rcParams.update(params)

# for check inorganic aerosol name
Nsps = ['BiMT','BC','MD','SOAlP','H2O','rganics','NH4','NA','SO4']
# default color
col10 = ['r','blue','orange','green','cyan','purple','brown','pink','gray','olive']

# add vertical lines for time break points
Add_vlines = [28800, 259500] #8 hour and 72hour # False
Add_xticks = '%d-%Hh' # False

def get_time_array(path):
    """return the time array from the input path (only one path)
        MUST HAVE THE FILE: report.txt"""
    #return np.arange(0.,432001.,3600.) # for silend mode
    return np.arange(0.,Ttot[0]+0.1,DeltaT[0])

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
                data[1][i][k] = float([s for s in info[k].split(' ') if s != ''][0])
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

def display_total_conc(data,labels,inds = [],savpath=None,savname=None,lst=None,Type=None,errType='fe2',col = col10):
    """display the concs"""

    cind=[[0],[1],[0,1]]
    titles=['gas-phase','SOA','TM']
    # paths
    ns = len(data[1])

    # for plotting
    plotData=[]
    plottime=[]

    # read data
    if Type == 'plain':
        plottime = data[0]
        plotData = data[1]

    elif Type == 'yield': # soa yield
        ns = len(data)
        for i in range(ns):

            # timesteps
            nt = len(data[1][2])
            # number of speices
            nsp = len(data[i][0])

            plottime.append(data[1][2])
            plotData.append(np.zeros(nt))

            pvocs= []
            for k in range(nsp):
                if data[i][0][k] in NokeepSp: continue
                elif data[i][0][k] in primaryVOCs: pvocs.append(k)

                if np.sum(data[i][1][k][1]) != 0.0: 
                    for t in range(nt): plotData[-1][t] += data[i][1][k][1][t]

            initvoc, initsoa = sum([data[i][1][k][0][0] for k in pvocs]), plotData[-1][0]

            # generate SOA yield
            for t in range(nt):
                if t > 0: 
                    tmp = sum([data[i][1][k][0][t] for k in pvocs])
                    plotData[-1][t] = (plotData[-1][t]-initsoa)/(initvoc-tmp) * 100. # SOA - VOC

            print('--- max yield: ',round(np.max(plotData[-1]),1),'min: ',round(np.min(plotData[-1][1:]),1),'average: ',round(np.average(plotData[-1][1:]),1),'end: ',round(plotData[-1][-1],1))

    else:
        print('display_total_conc: can only read type plain or yield.',Type)
        raise TypeError

    # compute err
    if Type == 'plain':
        rmse = [0]
        for i in range(ns-1):
            rmse.append(error(data[1][i+1], data[1][0], errType))
    elif Type == 'yield': pass
    else: 
        print('SSR: type not recognized.', Type)
        raise TypeError

    if Type == 'plain':

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            for i in range(ns):
                if lst is None: ax1.plot(plottime,plotData[i],label=labels[i])
                else: ax1.plot(plottime,plotData[i],linestyle = lst[i],label=labels[i],color = col[i])

            if Add_vlines:# vertial line mark time
                for i in Add_vlines:plt.axvline(x=plottime[0]+i,color='gray',linestyle='--')
            ax1.set_xlabel(r"Time (hour)")
            ax1.set_ylabel(r'Concentration $ug/m^3$')

            if Add_xticks:
                a=[]
                a1=[]
                for i in range(6):
                    a.append(plottime[0]+i*86400)
                    a1.append(SecondToDate(i*86400,Type=Add_xticks))
                ax1.set_xlim(a[0],a[-1])
                plt.xticks(a, a1)

            if ns == 2: 
                plt.title('Time variation of total SOA with error {:6.4f}'.format(rmse[-1]))
                print('Total SOA: '+errType+' {:6.4f}\n'.format(rmse[-1]))
            else: plt.title('Time variation of total SOA')
            plt.legend(loc='best',framealpha = 0.5)
            plt.tight_layout(rect=[0, 0.05, 1, 0.95])

            if savpath is None:
                if savname is None: plt.savefig('total_org_{:d}.png'.format(ns),dpi=300)
                else: plt.savefig('total_org_{:s}.png'.format(savname))#,dpi=300)
            else :
                os.makedirs(savpath, exist_ok=True) 
                if savname is None: plt.savefig('{:s}/total_org_{:d}.png'.format(savpath,ns))#,dpi=300)
                else: plt.savefig('{:s}/total_org_{:s}.png'.format(savpath,savname))#,dpi=300)
            #plt.show()
            plt.close()

    elif Type == 'yield':

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            for i in range(ns):
                if lst is None: ax1.plot(plottime[i],plotData[i],label=labels[i])
                else: ax1.plot(plottime[i],plotData[i],linestyle = lst[i],label=labels[i],color = col[i])

            if Add_vlines:# vertial line mark time
                for i in Add_vlines:plt.axvline(x=plottime[0][0]+i,color='gray',linestyle='--')

            ax1.set_xlabel('Time (hour)')
            ax1.set_ylabel('SOA yield (%)')

            if Add_xticks:
                a=[]
                a1=[]
                for i in range(6):
                    a.append(plottime[0][0]+i*86400)
                    a1.append(SecondToDate(i*86400,Type=Add_xticks))
                ax1.set_xlim(a[0],a[-1])
                plt.xticks(a, a1)

            plt.legend(loc='best',framealpha = 0.5)
            plt.tight_layout(rect=[0, 0.05, 1, 0.95])

            if savpath is None:
                if savname is None: plt.savefig('total_yield_{:d}.png'.format(ns),dpi=300)
                else: plt.savefig('total_yield_{:s}.png'.format(savname),dpi=300)
            else :
                os.makedirs(savpath, exist_ok=True)
                if savname is None: plt.savefig('{:s}/total_yield_{:d}.png'.format(savpath,ns),dpi=300)
                else: plt.savefig('{:s}/total_yield_{:s}.png'.format(savpath,savname),dpi=300)

            plt.close()
            return 0.0
    return rmse

def error(val,ref,indicator,threshold=1E-20,points=None):
    """compute the error dependind on the input indicator"""
    
    # check length and type
    if np.array_equal(val, ref): 
        return 0.0

    if len(val) == len(ref):
        n = len(val)
        if n == 0 :
            print('EA: empty lists')
            return None
    else:
        print("EA: length not equal",len(val),len(ref))
        return None
    
    # if only check few points
    if points is not None:
        print('Extracting points with index: ',points)
        n = len(points)
        # extract points
        tmp = np.zeros((2,n))
        for i in range(n):
            tmp[0][i] = val[points[i]]
            tmp[1][i] = ref[points[i]]
        # rebuild val and ref
        val = tmp[0]
        ref = tmp[1]
        print('Built new val and ref.')

    # compute error depending on the indicator
    err = 0.

    # relative error
    if indicator == 're':
        for i in range(n):
           err += abs(val[i]-ref[i])/(ref[i]+threshold)
        err = err / n

    # fractional bias
    elif indicator == 'fb':
        denominator = 0.
        numerator = 0.
        for i in range(n):
            denominator += (val[i] + ref[i])
            numerator += (val[i] - ref[i])
        if denominator != 0. :
            err = 2 * numerator / denominator
        else:
            print('EA: denominator of fb is 0.')
            return None

    # fractional errors
    elif indicator == 'fe':
        denominator = 0.
        numerator = 0.
        for i in range(n):
            denominator += (val[i] + ref[i])
            numerator += abs(val[i] - ref[i])
        if denominator != 0. :
            err = 2 * numerator / denominator
        else:
            print('EA: denominator of fb is 0.')
            return None

    # max fe of day1 and day2 - day5
    elif 'fe' in indicator:

        # obtain number of sections
        nt = int(indicator.replace('fe',''))
        err = 0.0 # err of the section

        # store the number of index of this section
        if abs(nt) <= 5 and abs(nt) >= 2: 
            ns = [25,len(ref)]
            for n in range(abs(nt)-2): ns.insert(n+1,(len(ref)-25)/(abs(nt)-1)+ns[n])
        else:
            print('EA: nt is not 2,3,4,5',nt)
            return None

        ind = 0 # index
        for n in range(abs(nt)):
            denominator = 0.
            numerator = 0.
            tot = 0.0 # total ref value in the section
            while (ind < ns[n]):
                denominator += (val[ind] + ref[ind])
                numerator += abs(val[ind] - ref[ind])
                tot += ref[ind]
                ind += 1
            if denominator != 0. : 
                if nt > 0 : err = max(err, 2. * numerator / denominator * tot)
                else: err = max(err, 2. * numerator / denominator)
            else:
                print('EA: denominator of fb is 0.')
                return None

        if nt > 0 : err /= sum(ref)
    else:
        print('indicator is not found: ',indicator)
        return None

    return err
