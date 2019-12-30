
import h5py as h
import numpy as np

# response to an old Stack Overflow question of mine:
# https://stackoverflow.com/questions/33529057/indices-that-intersect-and-sort-two-numpy-arrays
def intersecting_indices(x, y):
    u_x, u_idx_x = np.unique(x, return_index=True)
    u_y, u_idx_y = np.unique(y, return_index=True)
    i_xy = np.intersect1d(u_x, u_y, assume_unique=True)
    i_idx_x = u_idx_x[np.in1d(u_x, i_xy, assume_unique=True)]
    i_idx_y = u_idx_y[np.in1d(u_y, i_xy, assume_unique=True)]
    return i_idx_x, i_idx_y

dirname = '/global/projecta/projectdirs/lsst/groups/WL/projects/wl-txpipe-hack/DESY1/'
#Load the files

s = h.File(dirname+'shear_catalog_desy1.h5', 'r')
p = h.File(dirname+'photometry_catalog.hdf5', 'r')

idxs, idxp = intersecting_indices(s['/metacal/objectId'].value, p['/photometry/objectId'].value)

# only use objects with metacal_flag = 0
mcal_mask = (s['/metacal/mcal_flags'].value[idxs]==0)

# save shear catalog
dnames = ['dec', 'mcal_T', 'mcal_T_1m', 'mcal_T_1p', 'mcal_T_2m', 'mcal_T_2p', 'mcal_flags', 'mcal_g1', 'mcal_g1_1m', 'mcal_g1_1p', 'mcal_g1_2m', 'mcal_g1_2p', 'mcal_g2', 'mcal_g2_1m', 'mcal_g2_1p', 'mcal_g2_2m', 'mcal_g2_2p', 'mcal_mag_err_i', 'mcal_mag_err_i_1m', 'mcal_mag_err_i_1p', 'mcal_mag_err_i_2m', 'mcal_mag_err_i_2p', 'mcal_mag_err_r', 'mcal_mag_err_r_1m', 'mcal_mag_err_r_1p', 'mcal_mag_err_r_2m', 'mcal_mag_err_r_2p', 'mcal_mag_err_z', 'mcal_mag_err_z_1m', 'mcal_mag_err_z_1p', 'mcal_mag_err_z_2m', 'mcal_mag_err_z_2p', 'mcal_mag_i', 'mcal_mag_i_1m', 'mcal_mag_i_1p', 'mcal_mag_i_2m', 'mcal_mag_i_2p', 'mcal_mag_r', 'mcal_mag_r_1m', 'mcal_mag_r_1p', 'mcal_mag_r_2m', 'mcal_mag_r_2p', 'mcal_mag_z', 'mcal_mag_z_1m', 'mcal_mag_z_1p', 'mcal_mag_z_2m', 'mcal_mag_z_2p', 'mcal_psf_T_mean', 'mcal_psf_g1', 'mcal_psf_g2', 'mcal_s2n', 'mcal_s2n_1m', 'mcal_s2n_1p', 'mcal_s2n_2m', 'mcal_s2n_2p', 'objectId', 'ra', 'snr_i', 'snr_r', 'snr_z', 'tilename', 'region', 'mean_z', 'mean_z_1p', 'mean_z_1m', 'mean_z_2p', 'mean_z_2m']

f = h.File('/global/projecta/projectdirs/lsst/groups/WL/projects/wl-txpipe-hack/DESY1/shear_catalog_desy1_merged.h5', 'w')
g = f.create_group('metacal')
for i in range(len(dnames)):
    g.create_dataset(dnames[i], data=s['/metacal/'+dnames[i]].value[idxs][mcal_mask], dtype=s['/metacal/'+dnames[i]].dtype)
f.close()

s.close()


# save photometry catalog
dnames = ['dec', 'g_mag', 'g_mag_err', 'i_mag', 'i_mag_err', 'objectId', 'r_mag', 'r_mag_err', 'ra', 'snr_g', 'snr_i', 'snr_r', 'snr_u', 'snr_y', 'snr_z', 'u_mag', 'u_mag_err', 'y_mag', 'y_mag_err', 'z_mag', 'z_mag_err'] 
 
f = h.File('/global/projecta/projectdirs/lsst/groups/WL/projects/wl-txpipe-hack/DESY1/photometry_catalog_desy1_merged.h5', 'w')
g = f.create_group('photometry')
for i in range(len(dnames)):
    g.create_dataset(dnames[i], data=p['/photometry/'+dnames[i]].value[idxp][mcal_mask], dtype=p['/photometry/'+dnames[i]].dtype)
f.close()
p.close()


