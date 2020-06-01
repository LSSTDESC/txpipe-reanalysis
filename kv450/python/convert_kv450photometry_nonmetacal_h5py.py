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
import numpy as np

kidsdir = '/global/cscratch1/sd/elp25/txpipe-reanalysis/kv450/data/' 

def getdata(filename):
    kids_shear = Table.read(kidsdir + filename, memmap=True) 
    return kids_shear

kidsfiles = ['KV450_G9_reweight_3x4x4_v2_good.cat', 'KV450_G12_reweight_3x4x4_v2_good.cat',
             'KV450_G15_reweight_3x4x4_v2_good.cat', 'KV450_G23_reweight_3x4x4_v2_good.cat',
             'KV450_GS_reweight_3x4x4_v2_good.cat']
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
kids_shearall.sort('SeqNr')
print('sorted data')

print('getting columns')
#Getting the photometry columns
tilename  = kids_shearall['SeqNr_field']
dec       = kids_shearall['DELTA_J2000']
g_mag     = kids_shearall['MAG_GAAP_g']
g_mag_err = kids_shearall['MAGERR_GAAP_g']
i_mag     = kids_shearall['MAG_GAAP_i']
i_mag_err = kids_shearall['MAGERR_GAAP_i']
objectId  = kids_shearall['SeqNr']
r_mag     = kids_shearall['MAG_GAAP_r']
r_mag_err = kids_shearall['MAGERR_GAAP_r']
ra        = kids_shearall['ALPHA_J2000']
snr_g     = 1.086/g_mag_err
snr_i     = 1.086/i_mag_err
snr_r     = 1.086/r_mag_err
u_mag     = kids_shearall['MAG_GAAP_u']
u_mag_err = kids_shearall['MAGERR_GAAP_u']
snr_u     = 1.086/u_mag_err
y_mag     = kids_shearall['MAG_GAAP_Y']   #Placeholder
y_mag_err = kids_shearall['MAGERR_GAAP_Y'] #Placeholder
snr_y     = 1.086/y_mag_err           #Placeholder
z_mag     = kids_shearall['MAG_GAAP_Z']    #Placeholder
z_mag_err = kids_shearall['MAGERR_GAAP_Z'] #Placeholder
snr_z     = 1.086/z_mag_err           #Placeholder

print('loaded columns')
#Dealing with unicode, string needs to be S12
#tilename = np.array([a.encode('utf8') for a in tilename])
#objectId = np.array([a.encode('utf8') for a in objectId])

data   = [tilename, dec, g_mag, g_mag_err, i_mag, i_mag_err, objectId, r_mag, r_mag_err, ra, snr_g, snr_i, snr_r, u_mag, u_mag_err, snr_u, y_mag, y_mag_err, snr_y, z_mag, z_mag_err, snr_z]

dnames = ['tilename', 'dec', 'g_mag', 'g_mag_err', 'i_mag', 'i_mag_err', 'objectId', 'r_mag', 'r_mag_err', 'ra', 'snr_g', 'snr_i', 'snr_r', 'u_mag', 'u_mag_err', 'snr_u', 'y_mag', 'y_mag_err', 'snr_y', 'z_mag', 'z_mag_err', 'snr_z']

print('saving file')
#Saving the h5 file...
outputdir = '/global/cscratch1/sd/elp25/txpipe-reanalysis/kv450/data/DR3/'
f = h5.File(outputdir + 'photometry_catalog_kv450_nonmetacal.h5', 'w')
g = f.create_group('photometry')
for i in range(len(data)):
    g.create_dataset(dnames[i], data=data[i], dtype=data[i].dtype)
f.close()




