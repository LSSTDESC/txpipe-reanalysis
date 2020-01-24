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

def star_catalog_formatter(data_hsc, output_path):
    star_data = {}
    # We specifically use the stars chosen for PSF measurement
    #star = data_hsc['icalib_psf_used']
    star = data_hsc['iclassification_extendedness'] == 0
    # General columns
    star_data['ra'] = data_hsc['ira'][star]
    star_data['dec'] = data_hsc['idec'][star]
    star_data['id'] = data_hsc['object_id'][star]
    for band in 'grizy':
        star_data[f'{band}_mag'] = data_hsc[f'{band}cmodel_mag'][star]
    star_data['u_mag'] = 99*np.ones(len(star_data['id']))
    star_data['u_mag_err'] = 99*np.ones(len(star_data['id']))
    kinds = [
            ('', 'measured_'),
            ('psf', 'model_')
        ]

    for in_name, out_name in kinds:
        # Pulling out the correct moment columns
        Ixx = data_hsc[f'ishape_hsm_{in_name}moments_11'][star]
        Iyy = data_hsc[f'ishape_hsm_{in_name}moments_22'][star]
        Ixy = data_hsc[f'ishape_hsm_{in_name}moments_12'][star]

        # Conversion of moments to e1, e2
        T = Ixx + Iyy
        e = (Ixx - Iyy + 2j * Ixy) / (Ixx + Iyy)
        e1 = e.real
        e2 = e.imag

        # save to output
        star_data[f'{out_name}e1'] = e1
        star_data[f'{out_name}e2'] = e2
        star_data[f'{out_name}T'] = T
    #_catalog_out = h5py.File(output_path, 'a')
    #_catalog_out.create_dataset('star_catalog', data=pd.DataFrame(star_data).to_records())
    #_catalog_out.flush()
    n = len(star_data['ra'])
    _f = setup_output(output_path, 'stars', star_data, star_data.keys(), n)
    write_output(_f, 'stars', star_data.keys(), 0, n, star_data)
    _f.flush()
    _f.close()

parser = argparse.ArgumentParser(description='Convert HSC catalog to TXPipe-friendly catalog')
parser.add_argument('--input-path', '-i', type=str, default=None, dest='input_path',
        help='Input path to convert')
parser.add_argument('--output-path', '-o', type=str, default=None, dest='output_path',
        help='Output path of converted shear catalog')
args = parser.parse_args()
data = fitsio.read(args.input_path)
star_catalog_formatter(data, args.output_path)

