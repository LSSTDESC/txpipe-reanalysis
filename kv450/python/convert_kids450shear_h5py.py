#!/usr/env python
"""
.. module:: convert_kids450shear_h5py
:synopsis: Script to get the KIDS 450 shear catalog and save as a h5 file.
.. moduleauthor:: Maria Elidaiana <mariaeli@brandeis.edu>
"""

import numpy as np
from astropy.table import Table, vstack
import h5py as h5

kidsdir = '/global/projecta/projectdirs/lsst/groups/WL/projects/wl-txpipe-hack/KIDS/'

def getdata(filename):
    kids_shear = Table.read(kidsdir + filename, memmap=True)
    return kids_shear

kidsfiles = ['KiDS_DR3.1_G9_ugri_shear.fits', 'KiDS_DR3.1_G12_ugri_shear.fits',
             'KiDS_DR3.1_G15_ugri_shear.fits', 'KiDS_DR3.1_G23_ugri_shear.fits',
             'KiDS_DR3.1_GS_ugri_shear.fits']
tables = []
for i in range(len(kidsfiles)):
    d = getdata(kidsfiles[i])
    tables.append(d)

#Joining the tables
kids_shearall = vstack(tables)

#Sorting by ID
kids_shearall.sort('ID')

#Getting the photometry columns
#tilename  = kids_shearall['KIDS_TILE']

dec        = kids_shearall['DECJ2000']              
mcal_T     = (kids_shearall['KRON_RADIUS']* 0.2)**2  #The pixel_scale for KIDS is 0.2 arcsec/pixel           
mcal_T_1m  = mcal_T     #Placeholder        
mcal_T_1p  = mcal_T     #Placeholder       
mcal_T_2m  = mcal_T     #Placeholder       
mcal_T_2p  = mcal_T     #Placeholder       
mcal_flags = kids_shearall['fitclass'] #It's lensfit flag for bad shear measurements       
mcal_g1    = kids_shearall['e1']         
mcal_g1_1m = mcal_g1    #Placeholder       
mcal_g1_1p = mcal_g1    #Placeholder      
mcal_g1_2m = mcal_g1    #Placeholder     
mcal_g1_2p = mcal_g1    #Placeholder      
mcal_g2    = kids_shearall['e2']      
mcal_g2_1m = mcal_g2    #Placeholder      
mcal_g2_1p = mcal_g2    #Placeholder      
mcal_g2_2m = mcal_g2    #Placeholder      
mcal_g2_2p = mcal_g2    #Placeholder      
mcal_mag_err_i    = kids_shearall['MAGERR_i'] 
mcal_mag_err_i_1m = mcal_mag_err_i #Placeholder   
mcal_mag_err_i_1p = mcal_mag_err_i #Placeholder
mcal_mag_err_i_2m = mcal_mag_err_i #Placeholder
mcal_mag_err_i_2p = mcal_mag_err_i #Placeholder
mcal_mag_err_r    = kids_shearall['MAGERR_r']
mcal_mag_err_r_1m = mcal_mag_err_r #Placeholder
mcal_mag_err_r_1p = mcal_mag_err_r #Placeholder
mcal_mag_err_r_2m = mcal_mag_err_r #Placeholder
mcal_mag_err_r_2p = mcal_mag_err_r #Placeholder
mcal_mag_err_z    = mcal_mag_err_i #Placeholder   
mcal_mag_err_z_1m = mcal_mag_err_i #Placeholder 
mcal_mag_err_z_1p = mcal_mag_err_i #Placeholder 
mcal_mag_err_z_2m = mcal_mag_err_i #Placeholder 
mcal_mag_err_z_2p = mcal_mag_err_i #Placeholder 
mcal_mag_i    = kids_shearall['MAG_i']      
mcal_mag_i_1m = mcal_mag_i #Placeholder       
mcal_mag_i_1p = mcal_mag_i #Placeholder   
mcal_mag_i_2m = mcal_mag_i #Placeholder   
mcal_mag_i_2p = mcal_mag_i #Placeholder   
mcal_mag_r    = kids_shearall['MAG_r']      
mcal_mag_r_1m = mcal_mag_r #Placeholder  
mcal_mag_r_1p = mcal_mag_r #Placeholder  
mcal_mag_r_2m = mcal_mag_r #Placeholder  
mcal_mag_r_2p = mcal_mag_r #Placeholder  
mcal_mag_z    = mcal_mag_i #Placeholder        
mcal_mag_z_1m = mcal_mag_i #Placeholder  
mcal_mag_z_1p = mcal_mag_i #Placeholder  
mcal_mag_z_2m = mcal_mag_i #Placeholder  
mcal_mag_z_2p = mcal_mag_i #Placeholder  
mcal_psf_T_mean = (kids_shearall['FWHM_IMAGE']* 0.2)**2 #Is this right? Using pixel_scale = 0.2" 
mcal_psf_g1     = kids_shearall['PSF_e1']      
mcal_psf_g2     = kids_shearall['PSF_e2']      
mcal_s2n        = kids_shearall['pixel_SNratio'] #Is this right?         
mcal_s2n_1m     = mcal_s2n             #Placeholder       
mcal_s2n_1p     = mcal_s2n             #Placeholder  
mcal_s2n_2m     = mcal_s2n             #Placeholder 
mcal_s2n_2p     = mcal_s2n             #Placeholder 
mean_z          = kids_shearall['Z_B'] #BPZ best redshift estimate
objectId        = kids_shearall['ID']    
ra              = kids_shearall['RAJ2000'] 
region          = kids_shearall['Patch'] 
snr_i           = 1.086/mcal_mag_err_i
snr_r           = 1.086/mcal_mag_err_r
snr_z           = 1.086/mcal_mag_err_z #Placeholder 
tilename        = kids_shearall['KIDS_TILE']


#Dealing with unicode, string needs to be S12
tilename = np.array([a.encode('utf8') for a in tilename])
objectId = np.array([a.encode('utf8') for a in objectId])
region   = np.array([a.encode('utf8') for a in region])

#Saving the data as h5file
data = [dec, mcal_T, mcal_T_1m, mcal_T_1p, mcal_T_2m, mcal_T_2p, mcal_flags, mcal_g1, mcal_g1_1m, mcal_g1_1p, mcal_g1_2m, mcal_g1_2p, mcal_g2, mcal_g2_1m, mcal_g2_1p, mcal_g2_2m, mcal_g2_2p, mcal_mag_err_i, mcal_mag_err_i_1m, mcal_mag_err_i_1p, mcal_mag_err_i_2m, mcal_mag_err_i_2p, mcal_mag_err_r, mcal_mag_err_r_1m, mcal_mag_err_r_1p, mcal_mag_err_r_2m, mcal_mag_err_r_2p, mcal_mag_err_z, mcal_mag_err_z_1m, mcal_mag_err_z_1p, mcal_mag_err_z_2m, mcal_mag_err_z_2p, mcal_mag_i, mcal_mag_i_1m, mcal_mag_i_1p, mcal_mag_i_2m, mcal_mag_i_2p, mcal_mag_r, mcal_mag_r_1m, mcal_mag_r_1p, mcal_mag_r_2m, mcal_mag_r_2p, mcal_mag_z, mcal_mag_z_1m, mcal_mag_z_1p, mcal_mag_z_2m, mcal_mag_z_2p, mcal_psf_T_mean, mcal_psf_g1, mcal_psf_g2, mcal_s2n, mcal_s2n_1m, mcal_s2n_1p, mcal_s2n_2m, mcal_s2n_2p, objectId, ra, snr_i, snr_r, snr_z, tilename, region, mean_z]

dnames = ['dec', 'mcal_T', 'mcal_T_1m', 'mcal_T_1p', 'mcal_T_2m', 'mcal_T_2p', 'mcal_flags', 'mcal_g1', 'mcal_g1_1m', 'mcal_g1_1p', 'mcal_g1_2m', 'mcal_g1_2p', 'mcal_g2', 'mcal_g2_1m', 'mcal_g2_1p', 'mcal_g2_2m', 'mcal_g2_2p', 'mcal_mag_err_i', 'mcal_mag_err_i_1m', 'mcal_mag_err_i_1p', 'mcal_mag_err_i_2m', 'mcal_mag_err_i_2p', 'mcal_mag_err_r', 'mcal_mag_err_r_1m', 'mcal_mag_err_r_1p', 'mcal_mag_err_r_2m', 'mcal_mag_err_r_2p', 'mcal_mag_err_z', 'mcal_mag_err_z_1m', 'mcal_mag_err_z_1p', 'mcal_mag_err_z_2m', 'mcal_mag_err_z_2p', 'mcal_mag_i', 'mcal_mag_i_1m', 'mcal_mag_i_1p', 'mcal_mag_i_2m', 'mcal_mag_i_2p', 'mcal_mag_r', 'mcal_mag_r_1m', 'mcal_mag_r_1p', 'mcal_mag_r_2m', 'mcal_mag_r_2p', 'mcal_mag_z', 'mcal_mag_z_1m', 'mcal_mag_z_1p', 'mcal_mag_z_2m', 'mcal_mag_z_2p', 'mcal_psf_T_mean', 'mcal_psf_g1', 'mcal_psf_g2', 'mcal_s2n', 'mcal_s2n_1m', 'mcal_s2n_1p', 'mcal_s2n_2m', 'mcal_s2n_2p', 'objectId', 'ra', 'snr_i', 'snr_r', 'snr_z', 'tilename', 'region', 'mean_z']

f = h5.File(kidsdir + 'shear_catalog_kids450.h5', 'w')
g = f.create_group('metacal')
for i in range(len(data)):
    g.create_dataset(dnames[i], data=data[i], dtype=data[i].dtype)
f.close()


