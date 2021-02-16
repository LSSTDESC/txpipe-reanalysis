# txpipe-reanalysis

## Brief description:

This will be an extension to the WLPipe paper (1808.07335). We will feed public data from DES Y1, HSC-Y1 and KV450 into TXPipe and attempt to reproduce the measurements and systematics tests. If during the project some newer data becomes available (e.g. DES Y3, KiDS-1000), we will consider switching to the most recent catalogs. We will then feed the measurements into an inference code such as CosmoSIS or firecrown and the output will be compared with the published results. We will evaluate how different analysis choices could change the final constraints.

This project will hopefully help with the development of TXPipe and other DESC infrastructure by piping real data through the analysis pipeline.

## Data and Tools needed:
* TXPipe: https://github.com/LSSTDESC/TXPipe
* DES public data: https://des.ncsa.illinois.edu/home
* HSC public data: https://hsc.mtk.nao.ac.jp/ssp/data-release/
* KiDS public data: http://kids.strw.leidenuniv.nl/availability.php

## Overleaf 

https://www.overleaf.com/9264561458qbqwnkfswtqd 

## Cosmosis Chains 
https://www.dropbox.com/sh/b365wkr7m5w7k4t/AAChu-2nzGdRiNrS8WM7vMfGa?dl=0

## Usage

The scripts here will likely be using TXPipe, so it is likely that you will need to set the path before running 

`export PYTHONPATH=$PYTHONPATH:/path of your TXPipe installation/`

Otherwise running a TXPipe pipeline is the same as the TXPipe example. For convenience we copy the TXPipe example in this repo for people to check. That is, on nersc, run

`shifter --image docker:joezuntz/txpipe bash` 

`cd [path to your txpipe-reanalysis repo]`

`mkdir -p data/example/inputs`

`cd data/example/inputs`

`curl -O https://portal.nersc.gov/project/lsst/WeakLensing/shear_catalog.hdf5`

`curl -O https://portal.nersc.gov/project/lsst/WeakLensing/photometry_catalog.hdf5`

`curl -O https://portal.nersc.gov/project/lsst/WeakLensing/sample_cosmodc2_w10year_errors.dat`

`curl -O https://portal.nersc.gov/project/lsst/WeakLensing/cosmoDC2_trees_i25.3.npy`

`curl -O https://portal.nersc.gov/project/lsst/WeakLensing/exposures.hdf5`

`curl -O https://portal.nersc.gov/project/lsst/WeakLensing/star_catalog.hdf5`

`cd ../../..`

`ceci examples/laptop.yml`
