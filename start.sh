#! /bin/sh

module load python
pip uninstall parsl
pip uninstall ceci

export PYTHONPATH=$PYTHONPATH:/global/homes/c/chihway/TXPipe
shifter --image docker:joezuntz/txpipe bash
#module load python
#pip uninstall parsl
#pip uninstall ceci

