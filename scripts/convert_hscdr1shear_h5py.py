#!/usr/env python
"""
.. module:: convert_hscshear_h5py
:synopsis: Script to get the HSC shear info from the shear catalog and save as a h5 file.
"""

import numpy as np
import numpy.random as random
from astropy.table import Table, vstack
import h5py as h5

hscdir = '/global/cfs/projectdirs/lsst/groups/LSS/HSC_reanalysis/catalogs_05_22/'

def getdata(filename):
    hsc_shear = Table.read(hscdir + filename, memmap=True) 
    return hsc_shear

fields = ['GAMA09H', 'GAMA15H', 'HECTOMAP', 'VVDS', 'WIDE12H', 'XMM']
#fields = ['HECTOMAP']
tables = []
print('loading data')
for field in fields:
    d = getdata(f'{field.lower()}_meas.fits')
    hsc_shearall = d
    ngal = len(hsc_shearall)
    print(f'Read {ngal} sources from field: {field}')
    mask_check = (hsc_shearall['photoz_best'] >= 0.3) & (hsc_shearall['photoz_best'] <= 0.6)
    print('z-cut', np.count_nonzero(mask_check))
    regauss_flag = (hsc_shearall['ishape_hsm_regauss_flags']==False)
    mask_check &= regauss_flag
    print('regauss-flag', np.count_nonzero(mask_check))
    regaus_sigma_nan = (np.isnan(hsc_shearall['ishape_hsm_regauss_sigma'])==False)
    mask_check &= regaus_sigma_nan
    print('regauss-sigma-nan', np.count_nonzero(mask_check))
    extendedness = (hsc_shearall['iclassification_extendedness']!=0)
    mask_check &= extendedness
    print('ext-cut', np.count_nonzero(mask_check))
    flux_cmodel = (hsc_shearall['icmodel_flux']/hsc_shearall['icmodel_flux_err']>=10)
    mask_check &= flux_cmodel
    print('snr-i-cut', np.count_nonzero(mask_check))
    regauss_resolution = (hsc_shearall['ishape_hsm_regauss_resolution']>=0.3)
    mask_check &= regauss_resolution
    print('reso-cut', np.count_nonzero(mask_check))
    regauss_e = ((hsc_shearall['ishape_hsm_regauss_e1']**2+hsc_shearall['ishape_hsm_regauss_e2']**2)**(1/2)<2)
    mask_check &= regauss_e
    print('e1-e2-cut', np.count_nonzero(mask_check))
    efinite = (np.isfinite(hsc_shearall['ishape_hsm_regauss_e1'])) & (np.isfinite(hsc_shearall['ishape_hsm_regauss_e2']))
    mask_check &= efinite
    print('finite-e-cut', np.count_nonzero(mask_check))
    regauss_sigma_cut1 = (hsc_shearall['ishape_hsm_regauss_sigma'] >= 0)
    mask_check &= regauss_sigma_cut1
    print('sigma > 0 -cut', np.count_nonzero(mask_check))
    regauss_sigma_cut2 = (hsc_shearall['ishape_hsm_regauss_sigma']<=0.4)
    mask_check &= regauss_sigma_cut2
    print('sigma <= 0.4 -cut', np.count_nonzero(mask_check))
    imag_cut = (hsc_shearall['icmodel_mag']-hsc_shearall['a_i']<=27)
    mask_check &= imag_cut
    print('imag-cut', np.count_nonzero(mask_check))
    blendedness_abs_flux = (hsc_shearall['iblendedness_abs_flux']< 10**(-0.375))
    mask_check &= blendedness_abs_flux
    print('blendedness-cut', np.count_nonzero(mask_check))
    gflux_cut = (hsc_shearall['gcmodel_flux']/hsc_shearall['gcmodel_flux_err'])>=5
    rflux_cut = (hsc_shearall['rcmodel_flux']/hsc_shearall['rcmodel_flux_err'])>=5
    zflux_cut = (hsc_shearall['zcmodel_flux']/hsc_shearall['zcmodel_flux_err'])>=5
    yflux_cut = (hsc_shearall['ycmodel_flux']/hsc_shearall['ycmodel_flux_err'])>=5
    snr_all = 1.0*(gflux_cut) + 1.0*(rflux_cut) + 1.0*(zflux_cut) + 1.0*(yflux_cut)
    snr_all_cut = (snr_all >= 0.0)
    mask_check &= snr_all_cut
    print('snr-all-cut', np.count_nonzero(mask_check))
    ishape_hsm_regauss_e1_isnull = hsc_shearall['ishape_hsm_regauss_e1_isnull'] == False
    mask_check &= ishape_hsm_regauss_e1_isnull
    print('e1-isnull-cut', np.count_nonzero(mask_check))
    ishape_hsm_regauss_e2_isnull = hsc_shearall['ishape_hsm_regauss_e2_isnull'] == False
    mask_check &= ishape_hsm_regauss_e2_isnull
    print('e2-isnull-cut', np.count_nonzero(mask_check))
    ishape_hsm_regauss_derived_shape_weight_isnull = hsc_shearall['ishape_hsm_regauss_derived_shape_weight_isnull'] == False
    mask_check &= ishape_hsm_regauss_derived_shape_weight_isnull
    print('hsm_weight_isnull', np.count_nonzero(mask_check))
    ishape_hsm_regauss_derived_shear_bias_m_isnull = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_m_isnull'] == False
    mask_check &= ishape_hsm_regauss_derived_shear_bias_m_isnull
    print('hsm-m-isnull', np.count_nonzero(mask_check))
    ishape_hsm_regauss_derived_shear_bias_c1_isnull = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_c1_isnull'] == False
    mask_check &= ishape_hsm_regauss_derived_shear_bias_c1_isnull
    print('hsm-c1-isnull', np.count_nonzero(mask_check))
    ishape_hsm_regauss_derived_shear_bias_c2_isnull = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_c2_isnull'] == False
    mask_check &= ishape_hsm_regauss_derived_shear_bias_c2_isnull
    print('hsm-c2-isnull', np.count_nonzero(mask_check))
    ishape_hsm_regauss_derived_sigma_e_isnull = hsc_shearall['ishape_hsm_regauss_derived_sigma_e_isnull'] == False
    mask_check &= ishape_hsm_regauss_derived_sigma_e_isnull
    print('hsm-sigma-e-isnull', np.count_nonzero(mask_check))

    all_cuts = (regauss_flag) & (regaus_sigma_nan) & (extendedness) & (flux_cmodel) & (regauss_resolution) & (regauss_e) & (regauss_sigma_cut1) & (regauss_sigma_cut2) & (blendedness_abs_flux) & (snr_all_cut) & (imag_cut) & (efinite)
    mask_check = (hsc_shearall['photoz_best'] >= 0.3) & (hsc_shearall['photoz_best'] <= 0.6)
    print('rechecking cuts', np.count_nonzero(mask_check))
    print('all-cuts', np.count_nonzero((all_cuts) & (mask_check)))
    null_cuts = (ishape_hsm_regauss_e1_isnull) & (ishape_hsm_regauss_e2_isnull) & (ishape_hsm_regauss_derived_shape_weight_isnull) & (ishape_hsm_regauss_derived_shear_bias_m_isnull) & (ishape_hsm_regauss_derived_shear_bias_c1_isnull) & (ishape_hsm_regauss_derived_shear_bias_c2_isnull) &  (ishape_hsm_regauss_derived_sigma_e_isnull)
    print('null-cuts', np.count_nonzero((null_cuts) & (mask_check)))
    print('all+null-cuts', np.count_nonzero((null_cuts) & (all_cuts) & (mask_check)))
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
    mag_i    = hsc_shearall['icmodel_mag'] - hsc_shearall['a_i']     
    mag_r    = hsc_shearall['rcmodel_mag'] - hsc_shearall['a_r']
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
    s2n        = hsc_shearall['icmodel_flux']/hsc_shearall['icmodel_flux_err']        
    mean_z = hsc_shearall['photoz_best']
    #mean_z          = hsc_shearall['pz_best_eab'] 
    objectId        = hsc_shearall['object_id']    
    ra              = hsc_shearall['ra']
    snr_i           = hsc_shearall['icmodel_flux']/hsc_shearall['icmodel_flux_err']
    snr_r           = hsc_shearall['rcmodel_flux']/hsc_shearall['rcmodel_flux_err']
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
    print(f'Selected {len(dec)} objects')
    outputdir = '/global/cfs/projectdirs/lsst/groups/LSS/HSC_reanalysis/data_javi/2023_reanalysis/'
    print('saving file, ',outputdir + f'shear_sourcecatalog_hsc_{field}_nonmetacal_05_22.h5')

    f = h5.File(outputdir + f'shear_sourcecatalog_hsc_{field}_nonmetacal_05_22.h5', 'w')
    g = f.create_group('shear')
    for i in range(len(data)):
        g.create_dataset(dnames[i], data=data[i], dtype=data[i].dtype)
    metadata = {'catalog_type':'hsc'}
    g.attrs.update(metadata)
    # write in a group for information on the catalog type 
    f.close()
