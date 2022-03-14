""""

Mead et al. 2015 shows the two parameters of the model:
eta0 = a0 + a1*A

"""


import numpy as np

from cosmosis.datablock import option_section, names
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    a0 = options.get_double(option_section, "a_0", default=0.98)
    a1 = options.get_double(option_section, "a_1", default=-0.12)
    config = [a0, a1]
    return config

def execute(block, config):
    a0, a1 = config
    A = block[names.halo_model_parameters, "A"]
    eta = a0 + a1*A
    block[names.halo_model_parameters, "eta"] = eta
    return 0
