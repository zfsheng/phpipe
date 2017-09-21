#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 16:23:36 2017

@author: sf
"""


from astropy.io import fits
import os 
import fnmatch
import numpy as np
#import pdb
import subprocess
#import shutil
def mkdir(path):  

    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False

#RN04_170902_070000_022.fits
OPATH='coadded3'
mkdir(OPATH)

PATH='output_wcs'
Rfile=[]
Bfile=[]
for f in os.listdir(PATH):
    if fnmatch.fnmatch(f, 'B*.fits'):
        Rfile.append(f)
    else:
        Bfile.append(f)





def file_split(filelist):
    filelist.sort()
    f4L=[]

    for f in filelist:
       f4L.append(f[0:4])
       
    split_index_f4=[0]
    k=0
    tmp = f4L[0]
    for f in f4L:
        
       if f == tmp:
           k=k+1
       else:   # find the split index if there're different filters
           split_index_f4.append(k)
           tmp =f 
           k=k+1
    if split_index_f4==[0]: # all the files are from the same filter
        file_List=filelist
        single_filter_flag=1   
    else:        
        file_List=[]
        single_filter_flag=0
        for i in xrange(np.size(split_index_f4)-1):  #[0:i] ....[i:i+1], i+1==np.size() --> i == size-1
           print i     
           file_List.append( filelist[ split_index_f4[i]:split_index_f4[i+1] ]  )  
        
        file_List.append( filelist[ split_index_f4[i+1]: ] )
    return file_List,single_filter_flag





def obj_split(subfileL): #split file list based on different obj
    subfileL.sort()
    
    obj=[]
    for f in subfileL:
        obj.append(f[-8:-5]) # get the obj num
    
    split_obj_index=[0]
    k=0
    tmp = obj[0]
    for f in obj:
       if f == tmp:
           k=k+1
       else:
           split_obj_index.append(k)   # find index of the different obj
           tmp=f
           k=k+1
    if split_obj_index==[0]: # all the obj in subfileL are same
        file_List=subfileL
        single_obj_flag=1
    else:                   # there are different objs in the list
        file_List=[]
        single_obj_flag=0
        for i in xrange(np.size(split_obj_index)-1):
            file_List.append( subfileL[ split_obj_index[i]:split_obj_index[i+1]  ]  )      
    #    pdb.set_trace()
        file_List.append( subfileL[split_obj_index[i+1]: ] )    
    return file_List,single_obj_flag    

#####################################################################################################

fileList,single_filter=file_split(Rfile) 
 
if single_filter==1:  # there is just one filter
    subfileL=fileList   
    ff,flag=obj_split(subfileL)
    
    
    #dstdir =  os.path.join(dstroot, os.path.dirname(srcfile))
    
    if flag==0:
        for f in ff:
            num = np.size(f)
            if num<2:
                print 'less than 2 fits to co-add together'
        #        shutil.copy(PATH+'/'+f[0],O_PATH+'/'+f[0])
                continue
            h = fits.getheader(PATH+'/'+f[0])   
        
            outputName=OPATH+'/'+f[-1][:-4]+'coadd'
            coord_str = str(h['RA_HMS']) + ',' + str(h['DEC_DMS'])
            command = ''
            for i in xrange(num): 
               command=command+PATH+'/'+f[i]+' '
               
            command = command + ' -COMBINE Y -IMAGEOUT_NAME ' + outputName
            command = command + ' -CENTER_TYPE MANUAL -CENTER ' + coord_str
            command = command + ' -COMBINE_TYPE MEDIAN' +' -COPY_KEYWORDS MJD'
        #   command = command + ' -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE ' + str(pixel_scale)
        #   command = command + ' -IMAGE_SIZE ' + string(image_size)
        
            try:
              arg='SWarp '+command
              subprocess.call(arg,shell=True)
            except:
              arg='swarp '+command 
              subprocess.call(arg, shell=True)
    
    elif flag==1:
            f=ff
            num = np.size(f)
            if num<2:
                print 'less than 2 fits to co-add together'
        #        shutil.copy(PATH+'/'+f[0],O_PATH+'/'+f[0])
                exit
            h = fits.getheader(PATH+'/'+f[0])   
        
            outputName=OPATH+'/'+f[-1][:-4]+'coadd'
            coord_str = str(h['RA_HMS']) + ',' + str(h['DEC_DMS'])
            command = ''
            for i in xrange(num): 
               command=command+PATH+'/'+f[i]+' '
               
            command = command + ' -COMBINE Y -IMAGEOUT_NAME ' + outputName
            command = command + ' -CENTER_TYPE MANUAL -CENTER ' + coord_str
            command = command + ' -COMBINE_TYPE MEDIAN' +' -COPY_KEYWORDS MJD'
        #   command = command + ' -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE ' + str(pixel_scale)
        #   command = command + ' -IMAGE_SIZE ' + string(image_size)
        #
            try:
              arg='SWarp '+command
              subprocess.call(arg,shell=True)
            except:
              arg='swarp '+command 
              subprocess.call(arg, shell=True)      
    
    
else:
    for subfile in fileList: # there're more than one filter
    
        ff,flag=obj_split(subfileL)
        
        
        #dstdir =  os.path.join(dstroot, os.path.dirname(srcfile))
        
        if flag==0:
            for f in ff:
                num = np.size(f)
                if num<2:
                    print 'less than 2 fits to co-add together'
            #        shutil.copy(PATH+'/'+f[0],O_PATH+'/'+f[0])
                    continue
                h = fits.getheader(PATH+'/'+f[0])   
            
                outputName=OPATH+'/'+f[-1][:-4]+'coadd'
                coord_str = str(h['RA_HMS']) + ',' + str(h['DEC_DMS'])
                command = ''
                for i in xrange(num): 
                   command=command+PATH+'/'+f[i]+' '
                   
                command = command + ' -COMBINE Y -IMAGEOUT_NAME ' + outputName
                command = command + ' -CENTER_TYPE MANUAL -CENTER ' + coord_str
                command = command + ' -COMBINE_TYPE MEDIAN' +' -COPY_KEYWORDS MJD'
            #   command = command + ' -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE ' + str(pixel_scale)
            #   command = command + ' -IMAGE_SIZE ' + string(image_size)
            
                try:
                  arg='SWarp '+command
                  subprocess.call(arg,shell=True)
                except:
                  arg='swarp '+command 
                  subprocess.call(arg, shell=True)
        
        elif flag==1:
                f=ff
                num = np.size(f)
                if num<2:
                    print 'less than 2 fits to co-add together'
            #        shutil.copy(PATH+'/'+f[0],O_PATH+'/'+f[0])
                    exit
                h = fits.getheader(PATH+'/'+f[0])   
            
                outputName=OPATH+'/'+f[-1][:-4]+'coadd'
                coord_str = str(h['RA_HMS']) + ',' + str(h['DEC_DMS'])
                command = ''
                for i in xrange(num): 
                   command=command+PATH+'/'+f[i]+' '
                   
                command = command + ' -COMBINE Y -IMAGEOUT_NAME ' + outputName
                command = command + ' -CENTER_TYPE MANUAL -CENTER ' + coord_str
                command = command + ' -COMBINE_TYPE MEDIAN' +' -COPY_KEYWORDS MJD'
            #   command = command + ' -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE ' + str(pixel_scale)
            #   command = command + ' -IMAGE_SIZE ' + string(image_size)
            #
                try:
                  arg='SWarp '+command
                  subprocess.call(arg,shell=True)
                except:
                  arg='swarp '+command 
                  subprocess.call(arg, shell=True)             
        
        
        
        
        
