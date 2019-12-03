from astropy.io import fits
import h5py
import os
import numpy as np

def add_col(name):
    array = gold_sort[name]
    total_length = len(array)
    f.create_dataset('photometry/'+name,maxshape=(total_length,), 
                     shape=(total_length,), chunks=True, data = array) 
    #f['photometry/'+name][:total_length] = array


# List all files for gold catalog (public realease des y1), one file per tile
path = '/global/projecta/projectdirs/lsst/groups/WL/projects/wl-txpipe-hack/DESY1/gold/'
files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path,f)) and 'DES' in f]
print files
print len(files)

# Compile all tiles into a single array
gold = {}
gold['ra'] = []
gold['dec'] = []
gold['objectId'] = []
bands = 'ugrizy'
for band in bands:
    gold['%s_mag'%band] = []
    gold['%s_mag_err'%band] = []
    
for f in files:
    print f
    hdu = fits.open(path+f)
    d = hdu[1].data
    gold['ra'].extend(d['RA'])
    gold['dec'].extend(d['DEC'])
    gold['objectId'].extend(d['COADD_OBJECTS_ID'])
    for band in bands:
        if band == 'u':
            gold['%s_mag'%band].extend(d['MAG_AUTO_G'])
            gold['%s_mag_err'%band].extend(d['MAGERR_AUTO_G'])
        else:
            gold['%s_mag'%band].extend(d['MAG_AUTO_%s'%band.capitalize()])
            gold['%s_mag_err'%band].extend(d['MAGERR_AUTO_%s'%band.capitalize()])

for key in gold:
    gold[key] = np.array(gold[key])

# obtain snr for each band as S/N =  1.086/mag_err
for band in bands:
    gold['snr_%s'%band] = 1.086/gold['%s_mag_err'%band]

# Sort by id
ind = np.argsort(gold['objectId'])
gold_sort ={}
for key in gold:
    gold_sort[key] = gold[key][ind]

# Write output into a h5py file
outfile = path + 'photometry_catalog.hdf5'
f = h5py.File(outfile, 'w')
add_col('ra')
add_col('dec')
add_col('objectId')
for band in bands:
    add_col('%s_mag'%band)
    add_col('%s_mag_err'%band)
    add_col('snr_%s'%band)
f.close()

