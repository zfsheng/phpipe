#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 19:43:43 2017

@author: sf
"""

import astropy.io.fits as fits
import subprocess
import os


def mkdir(path):  

    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False


tpath='wcs.tmp' # tmp path which save  *.wcs 
mkdir(tpath)


def add_wcs(ra,dec,radius,input_image,output_wcs,tpath=tpath, arcsecperpix_low = 0.386,arcsecperpix_high = 0.406):
    arguments_def = ["solve-field", "--no-plots", "--dir", tpath, "--overwrite", "--index-xyls","none",
                     "--temp-axy", "--solved", "none", "--match", "none", "--rdls", "none", "--corr","none",
                     "--new-fits", output_wcs, input_image]
    
    arguments = arguments_def + ["--ra", str(ra), "--dec", str(dec), "--radius", str(radius),"--scale-units",
                                 "arcsecperpix", "--scale-low", str(arcsecperpix_low), "--scale-high",str(arcsecperpix_high)]
#    print arguments                                       #,"--cpulimit", "20"]
    subprocess.call(arguments)

### terminal cmd #####
#solve-field --ra 179.75 --dec -0.52 --radius 0.25 --scale-units arcsecperpix --scale-low 0.386 
#--scale-high 0.406 --no-plots --index-xyls none --temp-axy --solved none --match none --rdls none 
#--corr none  --new-fits testP5.new P5.fits

ipath='data3'  # data path
opath='output_wcs' # output path
mkdir(opath)
filelist=os.listdir(ipath)
filelist.sort()

for f in filelist:
    header = fits.getheader(ipath+'/'+f)
    ra_in = header['ra']
    dec_in = header['dec']
    radius = 0.25
   
    input_image = ipath+'/'+ f           # input image
    output_wcs = opath+'/'+f[:-4]+'fits'  # output file with wcs 
    add_wcs(ra_in,dec_in,radius,input_image,output_wcs)

    # delete *.wcs file in wcs.tmp
    try:
        os.remove(tpath+'/'+f[:-4]+'wcs')
    except OSError:
        pass




