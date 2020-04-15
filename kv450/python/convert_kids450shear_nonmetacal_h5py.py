#!/usr/env python
"""
.. module:: convert_kids450shear_h5py
:synopsis: Script to get the KIDS 450 shear catalog and save as a h5 file.
.. moduleauthor:: Maria Elidaiana <mariaeli@brandeis.edu>
"""

import numpy as np
from astropy.table import Table, vstack
import h5py as h5

kidsdir = '/global/projecta/projectdirs/lsst/groups/WL/projects/wl-txpipe-hack/KIDS/KiDS450_DR3/'

def getdata(filename):
    kids_shear = Table.read(kidsdir + filename, memmap=True)
    return kids_shear

kidsfiles = ['KiDS_DR3.1_G9_ugri_shear.fits', 'KiDS_DR3.1_G12_ugri_shear.fits',
             'KiDS_DR3.1_G15_ugri_shear.fits', 'KiDS_DR3.1_G23_ugri_shear.fits',
             'KiDS_DR3.1_GS_ugri_shear.fits']
tables = []
print('loading data')
for i in range(len(kidsfiles)):
    d = getdata(kidsfiles[i])
    tables.append(d)
print('loaded data')

#Joining the tables
kids_shearall = vstack(tables)

#Sorting by ID
print('sorting data')
kids_shearall.sort('ID')
print('sorted data')

print('getting columns')
#Getting the photometry columns
#tilename  = kids_shearall['KIDS_TILE']

dec        = kids_shearall['DECJ2000']              
T     = (kids_shearall['KRON_RADIUS']* 0.2)**2  #The pixel_scale for KIDS is 0.2 arcsec/pixel                    
flags = kids_shearall['fitclass'] #It's lensfit flag for bad shear measurements       
g1    = kids_shearall['e1']                  
g2    = kids_shearall['e2']         
mag_err_i    = kids_shearall['MAGERR_i'] 
mag_err_r    = kids_shearall['MAGERR_r']
mag_i    = kids_shearall['MAG_i']      
mag_r    = kids_shearall['MAG_r']        
psf_T_mean = (kids_shearall['FWHM_IMAGE']* 0.2)**2 #Is this right? Using pixel_scale = 0.2" 
psf_g1     = kids_shearall['PSF_e1']      
psf_g2     = kids_shearall['PSF_e2']      
s2n        = kids_shearall['pixel_SNratio'] #Is this right?         
mean_z          = kids_shearall['Z_B'] #BPZ best redshift estimate
objectId        = kids_shearall['ID']    
ra              = kids_shearall['RAJ2000'] 
region          = kids_shearall['Patch'] 
snr_i           = 1.086/mag_err_i
snr_r           = 1.086/mag_err_r
tilename        = kids_shearall['KIDS_TILE']

print('loaded columns')
#Dealing with unicode, string needs to be S12
tilename = np.array([a.encode('utf8') for a in tilename])
objectId = np.array([a.encode('utf8') for a in objectId])
region   = np.array([a.encode('utf8') for a in region])

#Saving the data as h5file
data = [dec, T, flags, g1, g2, mag_err_i, mag_err_r, mag_i, mag_r, psf_T_mean, psf_g1, psf_g2, s2n, mean_z, objectId, ra, region, snr_i, snr_r, tilename]

dnames = ['dec', 'T', 'flags', 'g1', 'g2', 'mag_err_i', 'mag_err_r', 'mag_i', 'mag_r', 'psf_T_mean', 'psf_g1', 'psf_g2', 's2n', 'mean_z', 'objectId', 'ra', 'region', 'snr_i', 'snr_r', 'tilename']

print('saving file')
outputdir = '/global/cscratch1/sd/elp25/txpipe-reanalysis/data/kids/'
f = h5.File(outputdir + 'shear_catalog_kids450_nonmetacal.h5', 'w')
g = f.create_group('shear')
for i in range(len(data)):
    g.create_dataset(dnames[i], data=data[i], dtype=data[i].dtype)
metadata = {'catalog_type':'lensfit'}
g.attrs.update(metadata)
# write in a group for information on the catalog type 
f.close()


