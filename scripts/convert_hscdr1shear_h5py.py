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

hscfiles = ['HSC_WIDE_GAMA09H.fits',
             'HSC_WIDE_GAMA15H.fits', 'HSC_WIDE_HECTOMAP.fits',
             'HSC_WIDE_VVDS_part1.fits','HSC_WIDE_VVDS_part2.fits',
             'HSC_WIDE_WIDE12H.fits','HSC_WIDE_XMM.fits']
fields = ['GAMA09H'] #'GAMA15H', 'HECTOMAP', 'VVDS', 'WIDE12H', 'XMM']
tables = []
print('loading data')
for field in fields:
    if field != 'VVDS':
        d = getdata(f'HSC_WIDE_{field}.fits')
    else:
        d1 = getdata(f'HSC_WIDE_{field}_part1.fits')
        d2 = getdata(f'HSC_WIDE_{field}_part2.fits')
        d = vstack([d1, d2])
    #tables.append(d)
    #Joining the tables
    #hsc_shearall = vstack(tables)
    hsc_shearall = d
    ngal = len(hsc_shearall)
    print(f'Read {ngal} sources from field: {field}')
    regauss_flag = (hsc_shearall['ishape_hsm_regauss_flags']==False)
    regaus_sigma_nan = (hsc_shearall['ishape_hsm_regauss_sigma'] !=np.nan)
    extendedness = (hsc_shearall['iclassification_extendedness']!=0)
    flux_cmodel = (hsc_shearall['icmodel_flux']/hsc_shearall['icmodel_flux_err']>=2)
    regauss_resolution = (hsc_shearall['ishape_hsm_regauss_resolution']>=0.3)
    regauss_e = ((hsc_shearall['ishape_hsm_regauss_e1']**2+hsc_shearall['ishape_hsm_regauss_e2']**2)**(1/2)<2)
    efinite = (np.isfinite(hsc_shearall['ishape_hsm_regauss_e1'])) & (np.isfinite(hsc_shearall['ishape_hsm_regauss_e2']))
    regauss_sigma_cut1 = (0<=hsc_shearall['ishape_hsm_regauss_sigma'])
    regauss_sigma_cut2 = (hsc_shearall['ishape_hsm_regauss_sigma']<=0.4)
    imag_cut = (hsc_shearall['icmodel_mag']-hsc_shearall['a_i']<=27)
    blendedness_abs_flux = (hsc_shearall['iblendedness_abs_flux']< 10**(-0.375))
    gflux_cut = (hsc_shearall['gcmodel_flux']/hsc_shearall['gcmodel_flux_err']>=2)
    rflux_cut = (hsc_shearall['rcmodel_flux']/hsc_shearall['rcmodel_flux_err']>=2)
    zflux_cut = (hsc_shearall['zcmodel_flux']/hsc_shearall['zcmodel_flux_err']>=2)
    yflux_cut = (hsc_shearall['ycmodel_flux']/hsc_shearall['ycmodel_flux_err']>=2)

    ishape_hsm_regauss_e1_isnull = hsc_shearall['ishape_hsm_regauss_e1_isnull'] == False
    ishape_hsm_regauss_e2_isnull = hsc_shearall['ishape_hsm_regauss_e2_isnull'] == False
    ishape_hsm_regauss_derived_shape_weight_isnull = hsc_shearall['ishape_hsm_regauss_derived_shape_weight_isnull'] == False
    ishape_hsm_regauss_derived_shear_bias_m_isnull = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_m_isnull'] == False
    ishape_hsm_regauss_derived_shear_bias_c1_isnull = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_c1_isnull'] == False
    ishape_hsm_regauss_derived_shear_bias_c2_isnull = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_c2_isnull'] == False
    ishape_hsm_regauss_derived_sigma_e_isnull = hsc_shearall['ishape_hsm_regauss_derived_sigma_e_isnull'] == False

    all_cuts = (regauss_flag) & (regaus_sigma_nan) & (extendedness) & (flux_cmodel) & (regauss_resolution) & (regauss_e) & (regauss_sigma_cut1) & (regauss_sigma_cut2) & (blendedness_abs_flux) & (gflux_cut) & (rflux_cut) & (zflux_cut) & (yflux_cut) & (imag_cut) & (efinite)
    null_cuts = (ishape_hsm_regauss_e1_isnull) & (ishape_hsm_regauss_e2_isnull) & (ishape_hsm_regauss_derived_shape_weight_isnull) & (ishape_hsm_regauss_derived_shear_bias_m_isnull) & (ishape_hsm_regauss_derived_shear_bias_c1_isnull) & (ishape_hsm_regauss_derived_shear_bias_c2_isnull) &  (ishape_hsm_regauss_derived_sigma_e_isnull)

    hsc_shearall = hsc_shearall[all_cuts&null_cuts]

    #Sorting by ID
    print('sorting data')
    hsc_shearall.sort('object_id')
    print('sorted data')

    print('getting columns')

    dec        = hsc_shearall['dec']              
    T     = hsc_shearall['ishape_hsm_moments_11']+hsc_shearall['ishape_hsm_moments_22']                   
    flags = hsc_shearall['ishape_hsm_regauss_derived_shape_weight_isnull'] # hsm regauss failures
    g1 = hsc_shearall['ishape_hsm_regauss_e1']
    g2 = hsc_shearall['ishape_hsm_regauss_e2']
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
    e1 = (Ixx - Iyy) / T
    e2 = 2*Ixy / T
    psf_g1     = e1      
    psf_g2     = e2      
    s2n        = 1.086/hsc_shearall['icmodel_mag_err']        
    mean_z          = hsc_shearall['pz_best_eab'] 
    objectId        = hsc_shearall['object_id']    
    ra              = hsc_shearall['ra']
    snr_i           = 1.086/mag_err_i
    snr_r           = 1.086/mag_err_r
    weight = hsc_shearall['ishape_hsm_regauss_derived_shape_weight']
    m = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_m']
    c1 = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_c1']
    c2 = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_c2']
    sigma_e = hsc_shearall['ishape_hsm_regauss_derived_rms_e']
    flag = hsc_shearall['wl_fulldepth_fullcolor']
    print('loaded columns')
    #Dealing with unicode, string needs to be S12
    #tilename = np.array([a.encode('utf8') for a in tilename])
    #objectId = np.array([a.encode('utf8') for a in objectId])
    #region   = np.array([a.encode('utf8') for a in region])

    #Saving the data as h5file
    data = [dec, T, flags, g1, g2, mag_err_i, mag_err_r, mag_i, mag_r, psf_T_mean, psf_g1, psf_g2, s2n, mean_z, objectId, ra, snr_i, snr_r, weight, m, c1, c2, sigma_e, weight, mean_z, flag]

    dnames = ['dec', 'T', 'flags', 'g1', 'g2', 'mag_err_i', 'mag_err_r', 'mag_i', 'mag_r', 'psf_T_mean', 'psf_g1', 'psf_g2', 's2n', 'mean_z', 'objectId', 'ra','snr_i', 'snr_r', 'lensfit_weight', 'm', 'c1', 'c2', 'sigma_e', 'weight', 'redshift_true', 'wl_fulldepth_fullcolor']

    outputdir = '/global/cscratch1/sd/jsanch87/txpipe-reanalysis/hsc/data/'
    print('saving file, ',outputdir + f'shear_catalog_hsc_{field}_nonmetacal.h5')

    f = h5.File(outputdir + f'shear_catalog_hsc_{field}_nonmetacal.h5', 'w')
    g = f.create_group('shear')
    for i in range(len(data)):
        g.create_dataset(dnames[i], data=data[i], dtype=data[i].dtype)
    metadata = {'catalog_type':'hsc'}
    g.attrs.update(metadata)
    # write in a group for information on the catalog type 
    f.close()
