#!/usr/bin/python3

""" It produces mean seasonal bias plots (WRF-EOBs) for temperature
 (0.11) and checks the statistical significance of the differences. """


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from netCDF4 import num2date
from scipy import stats


seasons = ['DJF', 'MAM', 'JJA', 'SON']
var = ['tg', 'tx', 'tn']
names = ['Tmean', 'Tmax', 'Tmin']

path = '/run/media/irida/linux/ncfiles/temp/'

# tmean, tmax, tmin
for i in range(3):
    
    fullpath1 = path + 'wrf_mean/wrf_mon'+ var[i] + '_d02.nc'
    fullpath2 = path + 'eobs_mean/eobs_mon' + var[i] + '_d02.nc'
    
    # read data
    data1 = Dataset(fullpath1, 'r')
    data2 = Dataset(fullpath2, 'r')
    
    montemp_wrf = np.array(data1.variables[var[i]])
    montemp_eobs = np.array(data2.variables[var[i]])
    
    # nan values
    montemp_eobs[montemp_eobs<=-9999] = np.nan
    montemp_wrf[montemp_wrf<=-9999] = np.nan

    x = data1.variables['XLONG'][:,:]
    y = data1.variables['XLAT'][:,:]
    
    time2 = data2.variables['time']

    #convert time(s) to date andtime
    datetime=num2date(time2[:], time2.units)  
    np.array(datetime)

    ##########################################################
    
    seas_montemp_eobs = np.nan*np.zeros((4,57,len(x[:,0]),len(x[0,:])))
    seas_montemp_wrf = np.nan*np.zeros((4,57,len(x[:,0]),len(x[0,:])))
    
    # seperate season values
    s=0
    j=0
    while j<57:
        
        if j==0:
            seas_montemp_wrf[0,j+1:j+3] = montemp_wrf[s:s+2]
            seas_montemp_eobs[0,j+1:j+3] = montemp_eobs[s:s+2] # djf
        else:
            seas_montemp_wrf[0,j:j+3] = montemp_wrf[s-1:s+2]
            seas_montemp_eobs[0,j:j+3] = montemp_eobs[s-1:s+2] # djf

        # eobs    
        seas_montemp_eobs[1,j:j+3] = montemp_eobs[s+2:s+5] # mam
        seas_montemp_eobs[2,j:j+3] = montemp_eobs[s+5:s+8] # jja
        seas_montemp_eobs[3,j:j+3] = montemp_eobs[s+8:s+11] # son
        # wrf
        seas_montemp_wrf[1,j:j+3] = montemp_wrf[s+2:s+5] # mam
        seas_montemp_wrf[2,j:j+3] = montemp_wrf[s+5:s+8] # jja
        seas_montemp_wrf[3,j:j+3] = montemp_wrf[s+8:s+11] # son
        
        s = s + 12
        j = j + 3
        
    #######################################################################
        
    diff = np.nan*np.zeros((4,57,len(x[:,0]),len(x[0,:])))
    mean_bias = np.nan*np.zeros((4,len(x[:,0]),len(x[0,:])))
    
    for s in range(4):
        for k in range(57):

            # Kelvin
            seas_montemp_eobs[s,k] = seas_montemp_eobs[s,k] + 273.15
            seas_montemp_wrf[s,k] = seas_montemp_wrf[s,k] + 273.15
            
            # model-observations
            diff[s,k,:,:] = seas_montemp_wrf[s,k,:,:]-seas_montemp_eobs[s,k,:,:]

    # calculate mean bias MB
    mean_bias[:,:,:] = np.nanmean(diff[:,:,:,:], axis=1) # seasons,months,y,x
    
    print('MEAN BIAS')
    print(np.nanmax(mean_bias))
    print(np.nanmin(mean_bias))
    
    mean_seas_temp_wrf = np.nanmean(seas_montemp_wrf, axis=1)
    mean_seas_temp_eobs = np.nanmean(seas_montemp_eobs, axis=1)
    seas_std_wrf =  np.nanstd(seas_montemp_wrf, axis=1, ddof=1)
    seas_std_eobs =  np.nanstd(seas_montemp_eobs, axis=1, ddof=1)

    ########################################################
    # weltch test
    
    n = 57
    sp = np.sqrt((seas_std_eobs**2+seas_std_wrf**2)/n)

    t_value = (mean_seas_temp_wrf-mean_seas_temp_eobs)/sp
    df = (((seas_std_eobs**2 + seas_std_wrf**2)/n)**2)/((((seas_std_eobs**2/n)**2))/(n-1)+(((seas_std_wrf**2/n)**2))/(n-1))
    p_value = 2*(1-stats.t.cdf(abs(t_value), df))

    lat_sig = np.zeros((4,len(x[:,0]),len(x[0,:])))*np.nan
    lon_sig = np.zeros((4,len(x[:,0]),len(x[0,:])))*np.nan
    
    for season in range(4):
        for lon in range(len(x[0,:])):
            for lat in range(len(x[:,0])):
                if p_value[season,lat,lon]>=0.05 :
                    # no significant values
                    lon_sig[season,lat,lon] = x[lat,lon]
                    lat_sig[season,lat,lon] = y[lat,lon]
                    
    #################################################################

    # plot seasonal bias
    fig = plt.figure(figsize=(8,5))
 
    map = Basemap(llcrnrlon=-12., urcrnrlon=42., llcrnrlat=32., urcrnrlat=72., resolution='l')
    parallels = np.arange(30., 70.1, 10)
    meridians = np.arange(-10., 40.1, 10)
    map.drawparallels(parallels, labels=[1,1,1,1], fontsize=10, linewidth=0.2, dashes=[10,10])
    map.drawmeridians(meridians, labels=[1,1,1,1], fontsize=10, linewidth=0.2, dashes=[10,10])
    map.drawcoastlines(linewidth=0.5)
    map.drawmapboundary(linewidth=0.5)
    bounds = np.arange(-6., 6.5, 0.5)
    ticks = np.arange(-6, 6.5, 1)
    
    for j in range(4):

        contour = map.contourf(x, y, mean_bias[j,:,:], levels=bounds,  cmap='jet', extend='both')
        
        if j == 0:
            cbar = map.colorbar(contour, ticks=ticks, pad='10%')
            cbar.set_label('['+chr(176)+'C]', rotation=90, labelpad=8, fontsize=10)
            cbar.ax.tick_params(labelsize=10)
    
            
        map.scatter(lon_sig[j,:,:], lat_sig[j,:,:], c='dimgray', s=0.5, marker='+', alpha=0.5)
        plt.title(seasons[j], fontsize=20, y=1.08)
        plt.savefig('./vergina011_eobs_MB_'+ var[i] +'-'+seasons[j]+'.png', bbox_inches='tight', dpi=300)
       
    plt.close()
    
