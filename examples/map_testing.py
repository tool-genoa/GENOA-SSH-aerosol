#===================================================
#
# This is a file for post-processing the record file
#
#=================================================== 
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

#### user input
# path + name of the record file for testing.
testing_file = '../../../genoa-clean-SQT/toSSHA/B3_chem/'
# name for the output figure.
savname = 'testing'

err_tol = 0.32 # maximum error in the error bar
# number of the color bar for display the error distribution. 
# For each column, the error range will be n * (err_tol/nerr) to (n+1) * (err_tol/nerr)
nerr = 8

#### load map info
lat = np.load('../inputs_bcary/conditions/latitudes.npy')
lon = np.load('../inputs_bcary/conditions/longitudes.npy')
# shape of the data
nm, nt, nz, ny, nx = 12, 24, 9, 153, 143
nps = nm*ny*nx
y_min = lat[0]
x_min = lon[0]
Delta_y = lat[1]-lat[0]
Delta_x = lon[1]-lon[0]


#### read errors and locations
if os.path.exists(testing_file):
    print('find testing file: {:s}'.format(testing_file))
    with open (testing_file,'r') as f: err = f.read().splitlines()
    nout = 0
    ilerr = [[],[]] # [locs, [errs]]

    for i in err:
        if i[0:3] != 'Run' : continue
        tmp= i.split('/')
        vtmp = float(i.split('\t')[-1])

        if vtmp >= err_tol: 
            print(i)
            nout += 1

        for j in tmp: # location
            if 'm' in j and 'y' in j and 'x' in j: il = tmp.index(j)

        m= int(tmp[il].split('m')[1].split('y')[0])
        y= int(tmp[il].split('y')[1].split('x')[0])
        x= int(tmp[il].split('x')[1].split('_')[0])

        if [m,y,x] in ilerr[0]: ilerr[1][ilerr[0].index([m,y,x])].append(vtmp)
        else:
            ilerr[0].append([m, y, x])
            ilerr[1].append([vtmp])
else:
    raise FileNotFoundError('Not find the input testing record file: ',testing_file)

# compute errors
averr, merr = 0.0, 0.0 # init average, max
for i in ilerr[1]: 
    averr += sum(i) / (1.0 * len(i))
    merr = max(merr, max(i))
averr /= len(ilerr[1])
# print info
print('Total {:d} testing conditions. {:d} conditions are with an error >= {:f}.'.format(len(ilerr),nout,err_tol))
print('Average error: {:f}, Max error: {:f}'.format(averr,merr))

# setup Lambert Conformal basemap.
print('Plotting...')
m = Basemap(width=12000000,height=9000000,projection='lcc',
            resolution='c',lat_1=45.,lat_2=55,lat_0=50,lon_0=-107.)
fig = plt.figure()
# get map
ma = Basemap(projection = 'cyl',
                llcrnrlon = x_min - Delta_x / 2.,
                llcrnrlat = y_min - Delta_y / 2.,
                urcrnrlon = x_min + Delta_x / 2. + Delta_x * float(nx - 1),
                urcrnrlat = y_min + Delta_y / 2. + Delta_y * float(ny - 1),
                resolution = 'l', # suppress_ticks = False,
                area_thresh = 1000)
ma.shadedrelief()

parallels = np.arange(35.,70.,10)
ma.drawparallels(parallels,labels=[False,True,True,False], fontsize = 'x-large')
meridians = np.arange(-10,40,10)
ma.drawmeridians(meridians,labels=[True,False,False,True], fontsize = 'x-large')

items = [[],[],[]] # for plotting
nerr_bar = np.arange(0,err_tol * 100. + 0.1, err_tol * 100./ nerr)
for i in range(len(ilerr[0])):
    # month,lat,lon
    m,y,x = ilerr[0][i][0], ilerr[0][i][1], ilerr[0][i][2]
    # y
    tmp =lat[y]
    items[0].append(tmp)
    # x
    tmp =lon[x]
    items[1].append(tmp)
    # err
    tmp = ilerr[1][i]
    #items[2].append(max(tmp) * 100.)
    # get average error
    items[2].append(sum(tmp)/2. * 100.)

# print info
print('Max err: ', max(items[2]),'Ave err: ',sum(items[2])/len(items[2]),'Number of conditions: ',len(items[0]))

# plot errors
plt.scatter(x=items[1], y=items[0],c=items[2], cmap=plt.cm.get_cmap("hot_r", nerr),marker = 'o', lw = 0, s= 10)
cbar = plt.colorbar(orientation="horizontal", fraction=0.09, pad=0.06,label = 'fractional error (%)',ticks = nerr_bar)
plt.clim(vmin=0., vmax = err_tol * 100.)
# layout setting
plt.tight_layout()
plt.subplots_adjust(top = 0.98, right = 0.93, left = 0.04, bottom = 0.08) # without title
# save figure
plt.savefig(savname,dpi=300)
plt.show()
plt.close()

