import numpy as np
#from astropy.io import fits
from astropy.table import Table, vstack
import h5py as h5
import numpy as np

hscdir = '/global/cfs/cdirs/lsst/groups/LSS/HSC_data/'

def getdata(filename):
    hsc_photometry = Table.read(hscdir + filename, memmap=True) 
    return hsc_photometry 

hscfiles = [ 'HSC_WIDE_GAMA09H.fits', 
             'HSC_WIDE_GAMA15H.fits', 'HSC_WIDE_HECTOMAP.fits',
             'HSC_WIDE_VVDS_part1.fits','HSC_WIDE_VVDS_part2.fits',
             'HSC_WIDE_WIDE12H.fits','HSC_WIDE_XMM.fits']
#'HSC_WIDE_AEGIS.fits'

#for hscfile in hscfiles:
#    print('processing file',hscfile)

#tables = []
#print('loading data')
#for i in range(len(hscfiles)):
#    d = getdata(hscfiles[i])
#    tables.append(d)
#    print('loaded data')
    #Joining the tables
#hsc_shearall = vstack(tables)
fields = ['GAMA09H', 'GAMA15H', 'HECTOMAP', 'VVDS', 'WIDE12H', 'XMM']
tables = []
print('loading data')
for field in fields:
    if field != 'VVDS':
        d = getdata(f'HSC_WIDE_{field}.fits')
    else:
        d1 = getdata(f'HSC_WIDE_{field}_part1.fits')
        d2 = getdata(f'HSC_WIDE_{field}_part2.fits')
        d = vstack([d1, d2])
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
    z_mag     = hsc_shearall['zcmodel_mag']    #Placeholder
    z_mag_err = hsc_shearall['zcmodel_mag_err'] #Placeholder
    snr_z     = 1.086/z_mag_err           #Placeholder
    mean_z          = hsc_shearall['pz_mean_eab']
    sigma_e = hsc_shearall['ishape_hsm_regauss_derived_rms_e']
    redshift_true = hsc_shearall['pz_best_eab']
    extendedness = hsc_shearall['iclassification_extendedness']
    flag = hsc_shearall['wl_fulldepth_fullcolor']
    print('loaded columns')
    #Dealing with unicode, string needs to be S12
    #objectId = np.array([a.encode('utf8') for a in objectId])

    data   = [dec, g_mag, g_mag_err, i_mag, i_mag_err, objectId, r_mag, r_mag_err, ra, snr_g, snr_i, snr_r, u_mag, u_mag_err, snr_u, y_mag, y_mag_err, snr_y, z_mag, z_mag_err, snr_z, mean_z, sigma_e, redshift_true,extendedness, flag]

    dnames = ['dec', 'mag_g', 'mag_g_err', 'mag_i', 'mag_i_err', 'objectId', 'mag_r', 'mag_r_err', 'ra', 'snr_g', 'snr_i', 'snr_r', 'mag_u', 'mag_u_err', 'snr_u', 'mag_y', 'mag_y_err', 'snr_y', 'mag_z', 'mag_z_err', 'snr_z', 'mean_z', 'sigma_e','redshift_true','extendedness', 'wl_fulldepth_fullcolor']

    print('saving file')
    #Saving the h5 file...
    outputdir = '/global/cscratch1/sd/jsanch87/txpipe-reanalysis/hsc/data/'
    #f = h5.File(outputdir + 'hsc_photometry_catalog_'+hscfile[:-5]+'.hdf5', 'w')
    f = h5.File(outputdir + f'photometry_catalog_hsc_{field}_nonmetacal.h5', 'w')
    g = f.create_group('photometry')
    for i in range(len(data)):
        g.create_dataset(dnames[i], data=data[i], dtype=data[i].dtype)
    f.close()

