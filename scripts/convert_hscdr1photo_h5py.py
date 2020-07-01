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

def photometry_catalog_formatter(data_hsc, output_path, photoz):
    _data = dict()
    _data['id'] = data_hsc['object_id']
    _data['ra'] = data_hsc['ra']
    _data['dec'] = data_hsc['dec']
    for band in ['g','r','i','z','y']:
        _data[f'{band}_mag'] = data_hsc[f'{band}cmodel_mag']
        _data[f'{band}_mag_err'] = data_hsc[f'{band}cmodel_mag_err']
        _data[f'snr_{band}'] = 1.086/data_hsc[f'{band}cmodel_mag_err']
    _data['u_mag']=99*np.ones(len(_data['ra']))
    _data['u_mag_err']=99*np.ones(len(_data['ra']))
    #_catalog_out = h5py.File(output_path, 'w')
    #_catalog_out.create_dataset('photometry', data=pd.DataFrame(_data).to_records())
    #_catalog_out.flush()
    _data['redshift_true'] = data_hsc[photoz]
    n = len(_data['ra'])
    _f = setup_output(output_path, 'photometry', _data, _data.keys(), n)
    write_output(_f, 'photometry', _data.keys(), 0, n, _data)
    _f.flush()
    _f.close()

parser = argparse.ArgumentParser(description='Convert HSC catalog to TXPipe-friendly catalog')
parser.add_argument('--input-path', '-i', type=str, default=None, dest='input_path',
        help='Input path to convert')
parser.add_argument('--output-path', '-o', type=str, default=None, dest='output_path',
        help='Output path of converted shear catalog')
parser.add_argument('--pz-col', type=str, default='pz_mean_eab', dest='photoz',
        help='Photoz column name')
args = parser.parse_args()
data = fitsio.read(args.input_path)
photometry_catalog_formatter(data, args.output_path, args.photoz)

