#!/usr/env python
"""
.. module:: convert_kids450photometry_h5py
:synopsis: Script to get the KIDS 450 photometry from the shear catalog and save as a h5 file.
.. moduleauthor:: Maria Elidaiana <mariaeli@brandeis.edu>
"""

import numpy as np
#from astropy.io import fits
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
tilename  = kids_shearall['KIDS_TILE']
dec       = kids_shearall['DECJ2000']
g_mag     = kids_shearall['MAG_g']
g_mag_err = kids_shearall['MAGERR_g']
i_mag     = kids_shearall['MAG_i']
i_mag_err = kids_shearall['MAGERR_i']
objectId  = kids_shearall['ID']
r_mag     = kids_shearall['MAG_r']
r_mag_err = kids_shearall['MAGERR_r']
ra        = kids_shearall['RAJ2000']
snr_g     = 1.086/g_mag_err
snr_i     = 1.086/i_mag_err
snr_r     = 1.086/r_mag_err
u_mag     = kids_shearall['MAG_u']
u_mag_err = kids_shearall['MAGERR_u']
snr_u     = 1.086/u_mag_err
y_mag     = kids_shearall['MAG_i']    #Placeholder
y_mag_err = kids_shearall['MAGERR_i'] #Placeholder
snr_y     = 1.086/y_mag_err           #Placeholder
z_mag     = kids_shearall['MAG_i']    #Placeholder
z_mag_err = kids_shearall['MAGERR_i'] #Placeholder
snr_z     = 1.086/z_mag_err           #Placeholder

#Dealing with unicode, string needs to be S12
tilename = np.array([a.encode('utf8') for a in tilename])
objectId = np.array([a.encode('utf8') for a in objectId])

data   = [tilename, dec, g_mag, g_mag_err, i_mag, i_mag_err, objectId, r_mag, r_mag_err, ra, snr_g, snr_i, snr_r, u_mag, u_mag_err, snr_u, y_mag, y_mag_err, snr_y, z_mag, z_mag_err, snr_z]
dnames = ['tilename', 'dec', 'g_mag', 'g_mag_err', 'i_mag', 'i_mag_err', 'objectId', 'r_mag', 'r_mag_err', 'ra', 'snr_g', 'snr_i', 'snr_r', 'u_mag', 'u_mag_err', 'snr_u', 'y_mag', 'y_mag_err', 'snr_y', 'z_mag', 'z_mag_err', 'snr_z']

#Saving the h5 file...
f = h5.File(kidsdir + 'photometry_catalog_kids450.h5', 'w')
g = f.create_group('photometry')
for i in range(len(data)):
    g.create_dataset(dnames[i], data=data[i], dtype=data[i].dtype)
f.close()




