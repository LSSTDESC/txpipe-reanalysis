import numpy as np
#from astropy.io import fits
from astropy.table import Table, vstack
import h5py as h5
import numpy as np

#hscdir = '/global/cfs/projectdirs/lsst/groups/LSS/HSC_reanalysis/catalogs_05_22/'
hscdir = '/global/cfs/projectdirs/lsst/groups/LSS/HSC_reanalysis/data_javi/lens_sample_2023_pdr1/'

def getdata(filename):
    hsc_photometry = Table.read(hscdir + filename, memmap=True) 
    return hsc_photometry 
pdr1=True
fields = ['GAMA09H', 'GAMA15H', 'VVDS', 'WIDE12H', 'XMM']
tables = []
print('loading data')
for field in fields:
    d = getdata(f'{field.lower()}_pdr1.fits')
    hsc_shearall = d
    ngal = len(hsc_shearall)
    print(f'Read {ngal} sources from field: {field}')

    extendedness = (hsc_shearall['iclassification_extendedness']==1)
    flux_cmodel = (hsc_shearall['icmodel_flux']/hsc_shearall['icmodel_flux_err']>=10)
    imag_cut = (hsc_shearall['icmodel_mag']-hsc_shearall['a_i']<=24.5)
    blendedness_abs_flux = (hsc_shearall['iblendedness_abs_flux']< 10**(-0.375))
    gflux_cut = (hsc_shearall['gcmodel_flux']/hsc_shearall['gcmodel_flux_err']>=5)
    rflux_cut = (hsc_shearall['rcmodel_flux']/hsc_shearall['rcmodel_flux_err']>=5)
    zflux_cut = (hsc_shearall['zcmodel_flux']/hsc_shearall['zcmodel_flux_err']>=5)
    yflux_cut = (hsc_shearall['ycmodel_flux']/hsc_shearall['ycmodel_flux_err']>=5)
    snr_all = 1.0*(gflux_cut) + 1.0*(rflux_cut) + 1.0*(zflux_cut) + 1.0*(yflux_cut)
    #ishape_hsm_regauss_e1_isnull = hsc_shearall['ishape_hsm_regauss_e1_isnull'] == False
    #ishape_hsm_regauss_e2_isnull = hsc_shearall['ishape_hsm_regauss_e2_isnull'] == False
    #ishape_hsm_regauss_derived_shape_weight_isnull = hsc_shearall['ishape_hsm_regauss_derived_shape_weight_isnull'] == False 
    #ishape_hsm_regauss_derived_shear_bias_m_isnull = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_m_isnull'] == False 
    #ishape_hsm_regauss_derived_shear_bias_c1_isnull = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_c1_isnull'] == False 
    #ishape_hsm_regauss_derived_shear_bias_c2_isnull = hsc_shearall['ishape_hsm_regauss_derived_shear_bias_c2_isnull'] == False 
    #ishape_hsm_regauss_derived_sigma_e_isnull = hsc_shearall['ishape_hsm_regauss_derived_sigma_e_isnull'] == False 

    all_cuts = (extendedness) & (flux_cmodel) & (blendedness_abs_flux) & (snr_all>=2) & (imag_cut)
    if not pdr1:
        null_cuts = hsc_shearall['wl_fulldepth_fullcolor'] == True
    else:
        null_cuts = np.ones(len(hsc_shearall), dtype=bool)
    if not pdr1:
        mask_check = (hsc_shearall['pz_best_eab'] >= 0.15) & (hsc_shearall['pz_best_eab'] <= 1.5) & (all_cuts)
    else:
        mask_check = all_cuts
    print('bin 0 + all-cuts', np.count_nonzero((mask_check) & (null_cuts)))
    hsc_shearall = hsc_shearall[all_cuts&null_cuts]

    #Sorting by ID
    print('sorting data')
    hsc_shearall.sort('object_id')
    print('sorted data')

    print('getting columns')
    #Getting the photometry columns
    dec       = hsc_shearall['dec']
    g_mag     = hsc_shearall['gcmodel_mag']
    g_mag_err = hsc_shearall['gcmodel_mag_err']
    i_mag     = hsc_shearall['icmodel_mag'] - hsc_shearall['a_i']
    i_mag_err = hsc_shearall['icmodel_mag_err']
    objectId  = hsc_shearall['object_id']
    r_mag     = hsc_shearall['rcmodel_mag'] - hsc_shearall['a_r']
    r_mag_err = hsc_shearall['rcmodel_mag_err']
    ra        = hsc_shearall['ra']
    snr_g     = hsc_shearall['gcmodel_flux']/hsc_shearall['gcmodel_flux_err']
    snr_i     = hsc_shearall['icmodel_flux']/hsc_shearall['icmodel_flux_err']
    snr_r     = hsc_shearall['rcmodel_flux']/hsc_shearall['rcmodel_flux_err']
    u_mag     = hsc_shearall['gcmodel_mag']
    u_mag_err = hsc_shearall['gcmodel_mag_err']
    snr_u     = hsc_shearall['gcmodel_flux']/hsc_shearall['gcmodel_flux_err'] # Placeholder
    y_mag     = hsc_shearall['ycmodel_mag']   #Placeholder
    y_mag_err = hsc_shearall['ycmodel_mag_err'] #Placeholder
    snr_y     = hsc_shearall['ycmodel_flux']/hsc_shearall['ycmodel_flux_err']           #Placeholder
    z_mag     = hsc_shearall['zcmodel_mag']    #Placeholder
    z_mag_err = hsc_shearall['zcmodel_mag_err'] #Placeholder
    snr_z     = hsc_shearall['zcmodel_flux']/hsc_shearall['zcmodel_flux_err']           #Placeholder
    mean_z          = hsc_shearall['pz_best_eab']
    if not pdr1:
        sigma_e = hsc_shearall['ishape_hsm_regauss_derived_rms_e']
    else:
        sigma_e = np.ones(len(hsc_shearall))*-99.0
    redshift_true = hsc_shearall['pz_best_eab']
    extendedness = hsc_shearall['iclassification_extendedness']
    if not pdr1:
        flag = hsc_shearall['wl_fulldepth_fullcolor']
    else:
        flag = np.ones(len(hsc_shearall), dtype=bool)
    print('loaded columns')
    #Dealing with unicode, string needs to be S12
    #objectId = np.array([a.encode('utf8') for a in objectId])

    data   = [dec, g_mag, g_mag_err, i_mag, i_mag_err, objectId, r_mag, r_mag_err, ra, snr_g, snr_i, snr_r, u_mag, u_mag_err, snr_u, y_mag, y_mag_err, snr_y, z_mag, z_mag_err, snr_z, mean_z, sigma_e, redshift_true,extendedness, flag]

    dnames = ['dec', 'mag_g', 'mag_g_err', 'mag_i', 'mag_i_err', 'objectId', 'mag_r', 'mag_r_err', 'ra', 'snr_g', 'snr_i', 'snr_r', 'mag_u', 'mag_u_err', 'snr_u', 'mag_y', 'mag_y_err', 'snr_y', 'mag_z', 'mag_z_err', 'snr_z', 'mean_z', 'sigma_e','redshift_true','extendedness', 'wl_fulldepth_fullcolor']

    print('saving file')
    #Saving the h5 file...
    outputdir = '/global/cfs/projectdirs/lsst/groups/LSS/HSC_reanalysis/data_javi/2023_reanalysis/'
    f = h5.File(outputdir + f'photometry_lenscatalog_hsc_{field}_nonmetacal_pdr1.h5', 'w')
    g = f.create_group('photometry')
    for i in range(len(data)):
        g.create_dataset(dnames[i], data=data[i], dtype=data[i].dtype)
    f.close()

