import h5py
import fitsio
import numpy as np
import argparse

def setup_output(path, group, cat, cols, n):
    f = h5py.File(path, 'w')
    g = f.create_group(group)
    for name in cols:
        g.create_dataset(name, shape=(n,), dtype=cat[name].dtype)
    return f

def write_output(output_file, group_name, cols, start, end, data):
    g = output_file[group_name]
    for name in cols:
        g[name][start:end] = data[name]

def shear_catalog_formatter(data_hsc, output_path, delta_gamma=0.02):
    _data = dict()
    data_hsc = data_hsc[data_hsc['iclassification_extendedness']==1] # Selecting galaxies
    _data['id'] = data_hsc['object_id']
    _data['ra'] = data_hsc['ira']
    _data['dec'] = data_hsc['idec']
    _data['mcal_T'] = data_hsc['ishape_hsm_moments_11']+data_hsc['ishape_hsm_moments_22']
    _data['mcal_T_1p'] = data_hsc['ishape_hsm_moments_11']+data_hsc['ishape_hsm_moments_22']
    _data['mcal_T_1m'] = data_hsc['ishape_hsm_moments_11']+data_hsc['ishape_hsm_moments_22']
    _data['mcal_T_2p'] = data_hsc['ishape_hsm_moments_11']+data_hsc['ishape_hsm_moments_22']
    _data['mcal_T_2m'] = data_hsc['ishape_hsm_moments_11']+data_hsc['ishape_hsm_moments_22']
    _data['mcal_g1'] = data_hsc['ishape_hsm_regauss_e1']
    _data['mcal_g2'] = data_hsc['ishape_hsm_regauss_e2']
    _data['mcal_s2n'] = 1.086/data_hsc['icmodel_mag_err']
    _data['mcal_s2n_1p'] = 1.086/data_hsc['icmodel_mag_err']
    _data['mcal_s2n_2p'] = 1.086/data_hsc['icmodel_mag_err']
    _data['mcal_s2n_1m'] = 1.086/data_hsc['icmodel_mag_err']
    _data['mcal_s2n_2m'] = 1.086/data_hsc['icmodel_mag_err']
    _data['i_mag'] = data_hsc['icmodel_mag']
    _data['i_mag_err'] = data_hsc['icmodel_mag_err']
    _data['mcal_g1_1p'] = data_hsc['ishape_hsm_regauss_e1']+delta_gamma
    _data['mcal_g2_1p'] = data_hsc['ishape_hsm_regauss_e2']+delta_gamma
    _data['mcal_g1_1m'] = data_hsc['ishape_hsm_regauss_e1']
    _data['mcal_g2_1m'] = data_hsc['ishape_hsm_regauss_e2']
    _data['mcal_g1_2p'] = data_hsc['ishape_hsm_regauss_e1']+delta_gamma
    _data['mcal_g2_2p'] = data_hsc['ishape_hsm_regauss_e2']+delta_gamma
    _data['mcal_g1_2m'] = data_hsc['ishape_hsm_regauss_e1']
    _data['mcal_g2_2m'] = data_hsc['ishape_hsm_regauss_e2']
    _data['mcal_flags'] = data_hsc['ishape_hsm_regauss_derived_shape_weight_isnull']
    Ixx = data_hsc['ishape_hsm_psfmoments_11']
    Iyy = data_hsc['ishape_hsm_psfmoments_22']
    Ixy = data_hsc['ishape_hsm_psfmoments_12']
    _data['mean_z'] = data_hsc['photoz_mean']
    _data['mean_z_1p'] = _data['mean_z']
    _data['mean_z_1m'] = _data['mean_z']
    _data['mean_z_2p'] = _data['mean_z']
    _data['mean_z_2m'] = _data['mean_z']

    for band in ['g','r','i','z','y']:
        _data[f'mcal_mag_{band}'] = data_hsc[f'{band}cmodel_mag']
        _data[f'mcal_mag_err_{band}'] = data_hsc[f'{band}cmodel_mag_err']
        _data[f'mcal_mag_{band}_1p'] = data_hsc[f'{band}cmodel_mag']
        _data[f'mcal_mag_err_{band}_1p'] = data_hsc[f'{band}cmodel_mag_err']
        _data[f'mcal_mag_{band}_1m'] = data_hsc[f'{band}cmodel_mag']
        _data[f'mcal_mag_err_{band}_1m'] = data_hsc[f'{band}cmodel_mag_err']
        _data[f'mcal_mag_{band}_2m'] = data_hsc[f'{band}cmodel_mag']
        _data[f'mcal_mag_err_{band}_2m'] = data_hsc[f'{band}cmodel_mag_err']
        _data[f'mcal_mag_{band}_2p'] = data_hsc[f'{band}cmodel_mag']
        _data[f'mcal_mag_err_{band}_2p'] = data_hsc[f'{band}cmodel_mag_err']
    # Conversion of moments to e1, e2
    T = Ixx + Iyy
    e = (Ixx - Iyy + 2j * Ixy) / (Ixx + Iyy)
    e1 = e.real
    e2 = e.imag
    _data['mcal_psf_g1'] = e1
    _data['mcal_psf_g2'] = e2
    _data['mcal_psf_T_mean'] = T
    _data['mcal_mag'] = _data['i_mag']
    _data['mcal_mag_err'] = _data['i_mag_err']
    n = len(_data['ra'])
    _f = setup_output(output_path, 'metacal', _data, _data.keys(), n)
    write_output(_f, 'metacal', _data.keys(), 0, n, _data)
    _f.flush()
    _f.close()

parser = argparse.ArgumentParser(description='Convert HSC catalog to TXPipe-friendly catalog')
parser.add_argument('--input-path', '-i', type=str, default=None, dest='input_path',
        help='Input path to convert')
parser.add_argument('--output-path', '-o', type=str, default=None, dest='output_path',
        help='Output path of converted shear catalog')
args = parser.parse_args()
data = fitsio.read(args.input_path)
shear_catalog_formatter(data, args.output_path)


