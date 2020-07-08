#!/usr/env python
"""
.. module:: convert_hscshear_h5py
:synopsis: Script to get the HSC shear info from the shear catalog and save as a h5 file.
"""

import numpy as np
import numpy.random as random
from astropy.table import Table, vstack
import h5py as h5

hscdir = '/global/cfs/cdirs/lsst/groups/LSS/HSC_data/'

def getdata(filename):
    hsc_shear = Table.read(hscdir + filename, memmap=True) 
    return hsc_shear

hscfiles = ['HSC_WIDE_AEGIS.fits', 'HSC_WIDE_GAMA09H.fits',
             'HSC_WIDE_GAMA15H.fits', 'HSC_WIDE_HECTOMAP.fits',
             'HSC_WIDE_VVDS_part1.fits','HSC_WIDE_VVDS_part2.fits',
             'HSC_WIDE_WIDE12H.fits','HSC_WIDE_XMM.fits']
tables = []
print('loading data')
for i in range(len(hscfiles)):
    d = getdata(hscfiles[i])
    tables.append(d)
print('loaded data')
#Joining the tables
hsc_shearall = vstack(tables)


hsc_shearall = hsc_shearall[hsc_shearall['iclassification_extendedness']==1] # Selecting galaxies

#Sorting by ID
print('sorting data')
hsc_shearall.sort('object_id')
print('sorted data')

print('getting columns')

dec        = hsc_shearall['dec']              
T     = hsc_shearall['ishape_hsm_moments_11']+hsc_shearall['ishape_hsm_moments_22']                   
flags = hsc_shearall['ishape_hsm_regauss_derived_shape_weight_isnull'] # hsm regauss failures      
e1_d = hsc_shearall['ishape_hsm_regauss_e1']-hsc_shearall['ishape_hsm_regauss_derived_shear_bias_c1'] # Distortion definition
e2_d = hsc_shearall['ishape_hsm_regauss_e2']-hsc_shearall['ishape_hsm_regauss_derived_shear_bias_c2']
theta = np.arctan2(e2_d, e1_d)
e_mod_d = np.sqrt(e1_d**2+e2_d**2) # Distortion module (Miralda-Escud\'e 91)
e_mod_s = np.arctanh(e_mod_d) # Transforming to shear definition BJ02
e1_s = e_mod_s*np.cos(theta)
e2_s = e_mod_s*np.sin(theta)
g1 = e1_s
g2 = e2_s
mag_err_i    = hsc_shearall['icmodel_mag_err'] 
mag_err_r    = hsc_shearall['rcmodel_mag_err']
mag_i    = hsc_shearall['icmodel_mag']     
mag_r    = hsc_shearall['rcmodel_mag'] 
psf_T_mean = hsc_shearall['ishape_hsm_psfmoments_11']+ hsc_shearall['ishape_hsm_psfmoments_22']
# Conversion of moments to e1, e2
Ixx = hsc_shearall['ishape_hsm_psfmoments_11']
Iyy = hsc_shearall['ishape_hsm_psfmoments_22']
Ixy = hsc_shearall['ishape_hsm_psfmoments_12']
T = Ixx + Iyy
e = (Ixx - Iyy + 2j * Ixy) / (Ixx + Iyy + 2*np.sqrt(Ixx*Iyy-Ixy**2))
e1 = e.real
e2 = e.imag
psf_g1     = e1      
psf_g2     = e2      
s2n        = 1.086/hsc_shearall['icmodel_mag_err']        
mean_z          = hsc_shearall['pz_mean_eab'] 
objectId        = hsc_shearall['object_id']    
ra              = hsc_shearall['ra']
snr_i           = 1.086/mag_err_i
snr_r           = 1.086/mag_err_r
weight = hsc_shearall['ishape_hsm_regauss_derived_shape_weight']
m = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_m']
c1 = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_c1']
c2 = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_c2']
sigma_e = hsc_shearall['ishape_hsm_regauss_derived_sigma_e']

print('loaded columns')
#Dealing with unicode, string needs to be S12
#tilename = np.array([a.encode('utf8') for a in tilename])
#objectId = np.array([a.encode('utf8') for a in objectId])
#region   = np.array([a.encode('utf8') for a in region])

#Saving the data as h5file
data = [dec, T, flags, g1, g2, mag_err_i, mag_err_r, mag_i, mag_r, psf_T_mean, psf_g1, psf_g2, s2n, mean_z, objectId, ra, snr_i, snr_r, weight, m, c1, c2, sigma_e, weight, mean_z]

dnames = ['dec', 'T', 'flags', 'g1', 'g2', 'mag_err_i', 'mag_err_r', 'mag_i', 'mag_r', 'psf_T_mean', 'psf_g1', 'psf_g2', 's2n', 'mean_z', 'objectId', 'ra','snr_i', 'snr_r', 'lensfit_weight', 'm', 'c1', 'c2', 'sigma_e', 'weight', 'redshift_true']

outputdir = '/global/cscratch1/sd/jsanch87/txpipe-reanalysis/hsc/data/'
print('saving file, ',outputdir + 'shear_catalog_hsc_nonmetacal.h5')

f = h5.File(outputdir + 'shear_catalog_hsc_nonmetacal.h5', 'w')
g = f.create_group('shear')
for i in range(len(data)):
    g.create_dataset(dnames[i], data=data[i], dtype=data[i].dtype)
metadata = {'catalog_type':'lensfit'}
g.attrs.update(metadata)
# write in a group for information on the catalog type 
f.close()
