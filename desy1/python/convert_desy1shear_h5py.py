#!/usr/env python
"""
.. module:: convert_desy1shear_h5py
:synopsis: Script to convert DES Y1 shear catalog from fits to a h5 file.
.. moduleauthor:: Maria Elidaiana <mariaeli@brandeis.edu>
"""

import numpy as np
from astropy.io import fits
import h5py as h5

#Loading the data
desy1_mcal = fits.open('/global/projecta/projectdirs/lsst/groups/WL/projects/wl-txpipe-hack/DESY1/mcal-y1a1-combined-riz-unblind-v4-matched.fits', memmap=True) 
d = desy1_mcal[1].data

desy1_bpz = fits.open('/global/projecta/projectdirs/lsst/groups/WL/projects/wl-txpipe-hack/DESY1/mcal-y1a1-combined-griz-blind-v3-matched_BPZ.fits', memmap=True)
meanz = desy1_bpz[1].data['mean_z']
desy1_bpz_1p = fits.open('/global/projecta/projectdirs/lsst/groups/WL/projects/wl-txpipe-hack/DESY1/mcal-y1a1-combined-griz-blind-v3-matched_BPZ_1p.fits', memmap=True)
meanz_1p = desy1_bpz_1p[1].data['mean_z']
desy1_bpz_1m = fits.open('/global/projecta/projectdirs/lsst/groups/WL/projects/wl-txpipe-hack/DESY1/mcal-y1a1-combined-griz-blind-v3-matched_BPZ_1m.fits', memmap=True)
meanz_1m = desy1_bpz_1m[1].data['mean_z']
desy1_bpz_2p = fits.open('/global/projecta/projectdirs/lsst/groups/WL/projects/wl-txpipe-hack/DESY1/mcal-y1a1-combined-griz-blind-v3-matched_BPZ_2p.fits', memmap=True)
meanz_2p = desy1_bpz_2p[1].data['mean_z']
desy1_bpz_2m = fits.open('/global/projecta/projectdirs/lsst/groups/WL/projects/wl-txpipe-hack/DESY1/mcal-y1a1-combined-griz-blind-v3-matched_BPZ_2m.fits', memmap=True)
meanz_2m = desy1_bpz_2m[1].data['mean_z']

print(len(meanz), len(meanz_1p), len(meanz_1m), len(meanz_2p), len(meanz_2m))

#Sorting by coadd_id
index_sort = np.argsort(d['coadd_objects_id'])
d = d[index_sort]
meanz = meanz[index_sort]
meanz_1p = meanz_1p[index_sort]
meanz_1m = meanz_1m[index_sort]
meanz_2p = meanz_2p[index_sort]
meanz_2m = meanz_2m[index_sort]

#Getting the fluxes
flux_r     = d['flux_r']
flux_err_r = d['flux_err_r']
flux_r_1m  = d['flux_r_1m']
flux_r_1p  = d['flux_r_1p']
flux_r_2m  = d['flux_r_2m']
flux_r_2p  = d['flux_r_2p']
flux_i     = d['flux_i']
flux_err_i = d['flux_err_i']
flux_i_1m  = d['flux_i_1m']
flux_i_1p  = d['flux_i_1p']
flux_i_2m  = d['flux_i_2m']
flux_i_2p  = d['flux_i_2p']
flux_z     = d['flux_z']
flux_err_z = d['flux_err_z']
flux_z_1m  = d['flux_z_1m']
flux_z_1p  = d['flux_z_1p']
flux_z_2m  = d['flux_z_2m']
flux_z_2p  = d['flux_z_2p']

#Columns that we need for TXpipe
dec        = d['dec']
mcal_T     = d['size'] 
mcal_T_1m  = d['size_1m']  
mcal_T_1p  = d['size_1p']
mcal_T_2m  = d['size_2m']
mcal_T_2p  = d['size_2p']
mcal_flags = d['flags_select'] 
mcal_g1    = d['e1']
mcal_g1_1m = d['e1_1m']
mcal_g1_1p = d['e1_1p']
mcal_g1_2m = d['e1_2m']
mcal_g1_2p = d['e1_2p']
mcal_g2    = d['e2']
mcal_g2_1m = d['e2_1m']
mcal_g2_1p = d['e2_1p']
mcal_g2_2m = d['e2_2m']
mcal_g2_2p = d['e2_2p']

print('converting fluxes')
mcal_mag_err_i    = 1.086*(flux_err_i/flux_i)     #This comes from mag_err_i = (2.5*flux_err_i)/(flux_i*ln(10))
mcal_mag_err_i_1m = 1.086*(flux_err_i/flux_i_1m)  #There is no "flux_err_i_1m" -> I will use flux_err_i instead 
mcal_mag_err_i_1p = 1.086*(flux_err_i/flux_i_1p)  #same
mcal_mag_err_i_2m = 1.086*(flux_err_i/flux_i_2m)  #same  
mcal_mag_err_i_2p = 1.086*(flux_err_i/flux_i_2p)  #same
mcal_mag_err_r    = 1.086*(flux_err_r/flux_r)
mcal_mag_err_r_1m = 1.086*(flux_err_r/flux_r_1m)  #same
mcal_mag_err_r_1p = 1.086*(flux_err_r/flux_r_1p)  #same
mcal_mag_err_r_2m = 1.086*(flux_err_r/flux_r_2m)  #same 
mcal_mag_err_r_2p = 1.086*(flux_err_r/flux_r_2p)  #same 
mcal_mag_err_z    = 1.086*(flux_err_z/flux_z) 
mcal_mag_err_z_1m = 1.086*(flux_err_z/flux_z_1m)  #same  
mcal_mag_err_z_1p = 1.086*(flux_err_z/flux_z_1p)  #same  
mcal_mag_err_z_2m = 1.086*(flux_err_z/flux_z_2m)  #same  
mcal_mag_err_z_2p = 1.086*(flux_err_z/flux_z_2p)  #same  
mcal_mag_i    = 30.0 - 2.5*np.log10(flux_i)  
mcal_mag_i_1m = 30.0 - 2.5*np.log10(flux_i_1m) 
mcal_mag_i_1p = 30.0 - 2.5*np.log10(flux_i_1p)   
mcal_mag_i_2m = 30.0 - 2.5*np.log10(flux_i_2m) 
mcal_mag_i_2p = 30.0 - 2.5*np.log10(flux_i_2p) 
mcal_mag_r    = 30.0 - 2.5*np.log10(flux_r) 
mcal_mag_r_1m = 30.0 - 2.5*np.log10(flux_r_1m)  
mcal_mag_r_1p = 30.0 - 2.5*np.log10(flux_r_1p)
mcal_mag_r_2m = 30.0 - 2.5*np.log10(flux_r_2m) 
mcal_mag_r_2p = 30.0 - 2.5*np.log10(flux_r_2p)
mcal_mag_z    = 30.0 - 2.5*np.log10(flux_z)
mcal_mag_z_1m = 30.0 - 2.5*np.log10(flux_z_1m)  
mcal_mag_z_1p = 30.0 - 2.5*np.log10(flux_z_1p)
mcal_mag_z_2m = 30.0 - 2.5*np.log10(flux_z_2m) 
mcal_mag_z_2p = 30.0 - 2.5*np.log10(flux_z_2p)

print('grabbing rest of terms')
mcal_psf_T_mean = d['psf_size']  
mcal_psf_g1 = d['psf_e1'] 
mcal_psf_g2 = d['psf_e2']
mcal_s2n    = d['snr']
mcal_s2n_1m = d['snr_1m'] 
mcal_s2n_1p = d['snr_1p'] 
mcal_s2n_2m = d['snr_2m']
mcal_s2n_2p = d['snr_2p']
objectId    = d['coadd_objects_id']
ra          = d['ra'] 
snr_i       = 1.086/mcal_mag_err_i   
snr_r       = 1.086/mcal_mag_err_r   
snr_z       = 1.086/mcal_mag_err_z   
tilename    = d['tilename']        #This column does not exist in the example h5
region      = d['region']          #Neither this 

#Dealing with unicode, string needs to be S12  
tilename = np.array([a.encode('utf8') for a in tilename])

print('saving the data')

#Saving the data as h5file
data = [dec, mcal_T, mcal_T_1m, mcal_T_1p, mcal_T_2m, mcal_T_2p, mcal_flags, mcal_g1, mcal_g1_1m, mcal_g1_1p, mcal_g1_2m, mcal_g1_2p, mcal_g2, mcal_g2_1m, mcal_g2_1p, mcal_g2_2m, mcal_g2_2p, mcal_mag_err_i, mcal_mag_err_i_1m, mcal_mag_err_i_1p, mcal_mag_err_i_2m, mcal_mag_err_i_2p, mcal_mag_err_r, mcal_mag_err_r_1m, mcal_mag_err_r_1p, mcal_mag_err_r_2m, mcal_mag_err_r_2p, mcal_mag_err_z, mcal_mag_err_z_1m, mcal_mag_err_z_1p, mcal_mag_err_z_2m, mcal_mag_err_z_2p, mcal_mag_i, mcal_mag_i_1m, mcal_mag_i_1p, mcal_mag_i_2m, mcal_mag_i_2p, mcal_mag_r, mcal_mag_r_1m, mcal_mag_r_1p, mcal_mag_r_2m, mcal_mag_r_2p, mcal_mag_z, mcal_mag_z_1m, mcal_mag_z_1p, mcal_mag_z_2m, mcal_mag_z_2p, mcal_psf_T_mean, mcal_psf_g1, mcal_psf_g2, mcal_s2n, mcal_s2n_1m, mcal_s2n_1p, mcal_s2n_2m, mcal_s2n_2p, objectId, ra, snr_i, snr_r, snr_z, tilename, region, meanz, meanz_1p, meanz_1m, meanz_2p, meanz_2m]

dnames = ['dec', 'mcal_T', 'mcal_T_1m', 'mcal_T_1p', 'mcal_T_2m', 'mcal_T_2p', 'mcal_flags', 'mcal_g1', 'mcal_g1_1m', 'mcal_g1_1p', 'mcal_g1_2m', 'mcal_g1_2p', 'mcal_g2', 'mcal_g2_1m', 'mcal_g2_1p', 'mcal_g2_2m', 'mcal_g2_2p', 'mcal_mag_err_i', 'mcal_mag_err_i_1m', 'mcal_mag_err_i_1p', 'mcal_mag_err_i_2m', 'mcal_mag_err_i_2p', 'mcal_mag_err_r', 'mcal_mag_err_r_1m', 'mcal_mag_err_r_1p', 'mcal_mag_err_r_2m', 'mcal_mag_err_r_2p', 'mcal_mag_err_z', 'mcal_mag_err_z_1m', 'mcal_mag_err_z_1p', 'mcal_mag_err_z_2m', 'mcal_mag_err_z_2p', 'mcal_mag_i', 'mcal_mag_i_1m', 'mcal_mag_i_1p', 'mcal_mag_i_2m', 'mcal_mag_i_2p', 'mcal_mag_r', 'mcal_mag_r_1m', 'mcal_mag_r_1p', 'mcal_mag_r_2m', 'mcal_mag_r_2p', 'mcal_mag_z', 'mcal_mag_z_1m', 'mcal_mag_z_1p', 'mcal_mag_z_2m', 'mcal_mag_z_2p', 'mcal_psf_T_mean', 'mcal_psf_g1', 'mcal_psf_g2', 'mcal_s2n', 'mcal_s2n_1m', 'mcal_s2n_1p', 'mcal_s2n_2m', 'mcal_s2n_2p', 'objectId', 'ra', 'snr_i', 'snr_r', 'snr_z', 'tilename', 'region', 'mean_z', 'mean_z_1p', 'mean_z_1m', 'mean_z_2p', 'mean_z_2m'] 
 
f = h5.File('/global/cscratch1/sd/elp25/DES-Y1-Explore/data/shear_catalog_desy1.h5', 'w')
g = f.create_group('metacal')
for i in range(len(data)):
    g.create_dataset(dnames[i], data=data[i], dtype=data[i].dtype)
f.close()
desy1_mcal.close()

