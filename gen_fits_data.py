#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 14:47:20 2017

@author: sf
"""


from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

import scipy
import scipy.ndimage
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroplan import time_grid_from_range
#import datetime
import pdb
import astropy.units as u

#

#red          blue
#0	7750 	0	4460.57
#1	7850	    1	4518.12
#2	7900	    2	4546.9
#3	7950	    3	4575.68
#4	8000	    4	4604.46
#5	8050	    5	4633.23
#6	8417.6	6	5740.75
#7	8526.21	7	5814.82
#8	8580.52	8	5851.86
#9	8634.83	9	5888.9
#10	8689.13	10	5925.93
#11	8743.44	11	5962.97
#12	SDSS I	12	SDSS g
#13	SDSS z		

RED_FILTER=[
[1,	7750,'R'],
[2,	7850,'R'],
[3,	7900,'R'],
[4,	7950,'R'],
[5,	8000,'R'],
[6,	8050,'R'],
[7,	8417.6,'R'],
[8,	8526.21,'R'],
[9,	8580.52,'R'],
[10 ,8634.83,'R'],
[11,	8689.13,'R'],
[12,	8743.44,'R'],
[13,7625,'R','i'],
[14,	9134,'R','z']
]

BLUE_FILTER=[
[1,	4460.57,'B'],
[2,	4518.12,'B'],
[3,	4546.9,'B'],
[4,	4575.68,'B'],
[5,	4604.46,'B'],
[6,	4633.23,'B'],
[7,	5740.75,'B'],
[8,	5814.82,'B'],
[9,	5851.86,'B'],
[10,	5888.9,'B'],
[11,	5925.93,'B'],
[12,5962.97,'B'],
[13	,4770,'B','g']
]



#################################################################
#def mock_fits(ra,dec,EXP,time,prihdr,data,filter_data,num,path):
def mock_fits(ra,dec,EXP,time,prihdr,data,filter_data,path):
    # 00 01 ...12 13
    if filter_data[0]<10:
        fI = '0'+str(filter_data[0])
    else:
        fI = str(filter_data[0])
    
    # 0~11: N band; >11 broad band
    if filter_data[0]<12:
        fT = 'N'
    else: 
        fT = filter_data[3]
            
    tt = time.iso.replace(' ' ,'_').replace('-','').replace(':','')[2:-4]    
    #    'BN' + '00_'+tt+'_001'+'.fits' 
#    filename = filter_data[2] + fT + fI + '_'+tt + '_'+'{:03d}'.format(num) + '.fits'
    filename = filter_data[2] + fT + fI + '_'+tt  + '.fits'
    # header setting
    prihdr['CHANNEL'] = filter_data[2]   #sting.  B or R
    prihdr['F_T'] = fT            # string.  N or g,i,z
    prihdr['F_I'] = filter_data[0] #int.    filter index
    prihdr['F_B'] = filter_data[1] #float.   wavelength
    prihdr['exp']=EXP
    prihdr['ra']=ra
    prihdr['dec']=dec
    prihdr['ISO']=time.iso
    prihdr['MJD']=time.mjd
        
    prihdu=fits.PrimaryHDU(data,header=prihdr)
    hdulist = fits.HDUList(prihdu)
    hdulist.writeto(path+filename)

#######################################################################

def mock_fits_withN(ra,dec,EXP,time,prihdr,data,filter_data,num,path):
#def mock_fits(ra,dec,EXP,time,prihdr,data,filter_data,path):
    # 00 01 ...12 13
    if filter_data[0]<10:
        fI = '0'+str(filter_data[0])
    else:
        fI = str(filter_data[0])
    
    # 0~11: N band; >11 broad band
    if filter_data[0]<12:
        fT = 'N'
    else: 
        fT = filter_data[3]

  
            
    tt = time.iso.replace(' ' ,'_').replace('-','').replace(':','')[2:-4]    
    #    'BN' + '00_'+tt+'_001'+'.fits' 
    filename = filter_data[2] + fT + fI + '_'+tt + '_'+'{:03d}'.format(num) + '.fits'
#    filename = filter_data[2] + fT + fI + '_'+tt  + '.fits'
    # header setting
    prihdr['CHANNEL'] = filter_data[2]   #sting.  B or R
    prihdr['F_T'] = fT            # string.  N or g,i,z
    prihdr['F_I'] = filter_data[0] #int.    filter index
    prihdr['F_B'] = filter_data[1] #float.   wavelength
    prihdr['exp']=EXP
    prihdr['ra']=ra
    prihdr['dec']=dec
    prihdr['ISO']=time.iso
    prihdr['MJD']=time.mjd
        
    prihdu=fits.PrimaryHDU(data,header=prihdr)
    hdulist = fits.HDUList(prihdu)
    hdulist.writeto(path+filename)
######################################################################



hdu=fits.open('frame-i-000756-2-0427.fits')
#hdu.info()
wcs = WCS(hdu[0].header)
header=hdu[0].header

ra=header['CRVAL1']
dec=header['CRVAL2']

fig = plt.figure()
ax=fig.add_subplot(111, projection=wcs)
plt.imshow(hdu[0].data,origin='lower',vmin=-0.2,vmax=0.2, cmap=plt.cm.viridis)
plt.xlabel('RA')
plt.ylabel('Dec')

#==============================================================================
# fig = plt.figure()
# data=scipy.misc.imrotate(hdu[0].data,20,interp='nearest')
# im=plt.imshow(data, origin='lower', vmin=-0.2,vmax=0.2, cmap=plt.cm.tab20)
# #im.set_cmap('nipy_spectral')
# #ra=ax.coords[0]
# #ra.set_major_formatter('dd:mm:ss.s' )
#==============================================================================


#####################################################
c = SkyCoord(ra,dec,unit=(u.deg, u.deg))
obj='SDSS_J{0}{1}'.format(c.ra.to_string(unit=u.hourangle, sep='', precision=2, pad=True),
           c.dec.to_string(sep='', precision=2, alwayssign=True, pad=True))

obj='obj'
# Initialization
prihdr = fits.Header()
prihdr['TELESCOP']='1.2m'
prihdr['OBSERVER'] = 'szf_test'
prihdr['obj'] =('','obj name')
prihdr['cmd']=('X','command sent to filter-wheel box')
prihdr['channel']=('X','camera type (e.g. B=blue,R=red)')
prihdr['F_T']=('X','filter type (e.g. N=narrow band, g=SDSS g)')
prihdr['F_I']=(0,'index of filter')
prihdr['F_B']=(0,'wavelength (A)')
prihdr['exp']=(0,'Exposure time (seconds)')
prihdr['ra']=(179.755936891,'RA of obj (deg)')
prihdr['dec']=(-0.524622205172,'DEC of obj (deg)')
prihdr['ra_hms'] = (c.ra.to_string(unit=u.hourangle, sep=':'), 'RA of obj (hourangle)')
prihdr['dec_dms'] =( c.dec.to_string(unit=u.deg, sep=':'), 'DEC of obj (deg)')
################################################################



time_range = Time(["2017-09-03 20:00:00", "2017-09-04 08:00"])
#time1=Time("2016-03-04 12:04")
time_grid = time_grid_from_range(time_range, time_resolution=12*u.min)
#tt=datetime.datetime(2016, 11, 07, 10, 49, 53, 239,tzinfo=pytz.timezone('Asia/Urumqi'))


data=scipy.ndimage.rotate(hdu[0].data,-6,reshape=False)
fig = plt.figure()
plt.imshow(data, origin='lower',vmin=-0.2,vmax=0.2, cmap=plt.cm.viridis)


exp=1800
 
#import numpy as np
path='data3/'
for i in xrange(5):
#    II=np.random.randint(11)
    ff=[0,1,2,0,1,1,6,7]
    II=ff[i]
    print II
    prihdr['cmd']=str(II+1)
    exp=1800
    for j in xrange(5):
        jj=i*5+j
        print BLUE_FILTER[II]
        data=scipy.ndimage.rotate(hdu[0].data,i,reshape=False)
        err=np.random.normal(0,0.1,data.size) # generate error from N(0,0.1)
        err=err.reshape(data.shape)
        pdb.set_trace()
        data=data+err
#        mock_fits(ra,dec,exp,time_grid[jj],prihdr,data,RED_FILTER[II],path)
#        mock_fits(ra,dec,exp,time_grid[jj],prihdr,data,BLUE_FILTER[II],path)
        mock_fits_withN(ra,dec,exp,time_grid[jj],prihdr,data,RED_FILTER[II],i,path)
        mock_fits_withN(ra,dec,exp,time_grid[jj],prihdr,data,BLUE_FILTER[II],i,path)
