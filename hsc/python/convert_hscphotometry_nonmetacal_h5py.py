#!/usr/env python
"""
.. module:: convert_hscphotometry_h5py
:synopsis: Script to get the HSC photometry from the shear catalog and save as a h5 file.
"""

import numpy as np
#from astropy.io import fits
from astropy.table import Table, vstack
import h5py as h5
import numpy as np

hscdir = '/global/cfs/cdirs/lsst/groups/LSS/HSC_data/'

def getdata(filename):
    hsc_photometry = Table.read(hscdir + filename, memmap=True) 
    return hsc_photometry

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
#Getting the photometry columns
dec       = hsc_shearall['dec']
g_mag     = hsc_shearall['gcmodel_mag']
g_mag_err = hsc_shearall['gcmodel_mag_err']
i_mag     = hsc_shearall['icmodel_mag']
i_mag_err = hsc_shearall['icmodel_mag_err']
objectId  = hsc_shearall['object_id']
r_mag     = hsc_shearall['rcmodel_mag']
r_mag_err = hsc_shearall['rcmodel_mag_err']
ra        = hsc_shearall['ra']
snr_g     = 1.086/g_mag_err
snr_i     = 1.086/i_mag_err
snr_r     = 1.086/r_mag_err
u_mag     = hsc_shearall['gcmodel_mag']
u_mag_err = hsc_shearall['gcmodel_mag_err']
snr_u     = 1.086/u_mag_err
y_mag     = hsc_shearall['ycmodel_mag']   #Placeholder
y_mag_err = hsc_shearall['ycmodel_mag_err'] #Placeholder
snr_y     = 1.086/y_mag_err           #Placeholder
z_mag     = hsc_shearall['zcmodel_mag_err']    #Placeholder
z_mag_err = hsc_shearall['zcmodel_mag_err'] #Placeholder
snr_z     = 1.086/z_mag_err           #Placeholder
mean_z          = hsc_shearall['pz_mean_eab']
sigma_e = hsc_shearall['ishape_hsm_regauss_derived_sigma_e']

print('loaded columns')
#Dealing with unicode, string needs to be S12
#objectId = np.array([a.encode('utf8') for a in objectId])

data   = [dec, g_mag, g_mag_err, i_mag, i_mag_err, objectId, r_mag, r_mag_err, ra, snr_g, snr_i, snr_r, u_mag, u_mag_err, snr_u, y_mag, y_mag_err, snr_y, z_mag, z_mag_err, snr_z, mean_z, sigma_e]

dnames = ['dec', 'g_mag', 'g_mag_err', 'i_mag', 'i_mag_err', 'objectId', 'r_mag', 'r_mag_err', 'ra', 'snr_g', 'snr_i', 'snr_r', 'u_mag', 'u_mag_err', 'snr_u', 'y_mag', 'y_mag_err', 'snr_y', 'z_mag', 'z_mag_err', 'snr_z', 'mean_z', 'sigma_e']

print('saving file')
#Saving the h5 file...
outputdir = '/global/cscratch1/sd/elp25/txpipe-reanalysis/hsc/data/'
f = h5.File(outputdir + 'photometry_catalog_hsc_nonmetacal.h5', 'w')
g = f.create_group('photometry')
for i in range(len(data)):
    g.create_dataset(dnames[i], data=data[i], dtype=data[i].dtype)
f.close()




