"""

Errors in cosmic shear measurement can lead to a multiplicative factor
scaling the observed shear spectra.

This module scales the measured C_ell to account for that difference,
assuming model values of the multplicative factor m, either per bin or for all bins.

"""
from __future__ import print_function
from builtins import range
from cosmosis.datablock import names, option_section
import sys

warning_note_displayed = False

def setup(options):
            # This is an option - can set c_per_bin = T to get
        # a different c for each tomographic bin, or F to get
        # one global value
        c_per_bin = options.get_bool(option_section, "c_per_bin", True)
        cal_section = options.get_string(
            option_section, "cal_section", default="c_systematic")
        verbose = options.get_bool(option_section, "verbose", False)
        suffix = options.get_string(option_section, "suffix", "")
        print()
        print("The shear_c_bias module will use calibration values from {} and look for ".format(cal_section))
        
        if suffix:
            suffix = "_" + suffix
        
        sec_names = {
        "shear_shear_plus": "shear_xi_plus" + suffix,
        "shear_shear_minus": "shear_xi_minus" + suffix,
    }
        
        return c_per_bin, sec_names, cal_section, verbose

def get_nbins(block, section):
        print('SECTIONS',block.sections())
        if block.has_value(section, "nbin_a"):
            n_a = block[section, "nbin_a"]
            n_b = block[section, "nbin_b"]
        else:
            n_a = block[section, "nbin"]
            n_b = n_a
        return n_a, n_b


def calibrate_section(block, section, c_a, c_b, verbose):
    n_a = len(c_a)
    n_b = len(c_b)
    for i in range(n_a):
        for j in range(n_b):

            # Get existing Xi_theta
            xi_name = "bin_{}_{}".format(i + 1, j + 1)
            if block.has_value(section, xi_name):
                if verbose:
                    print("Calibrating {} bin {} {} by ({}) * ({}) = {}".format(section, i + 1, j + 1, c_a[i], c_b[j], (c_a[i]) * (c_b[j])))
                block[section, xi_name] += (c_a[i]) * (c_b[j])
            elif verbose:
                print("No {} bin {} {} to calibrate".format(section, i + 1, j + 1))


def calibrate_shear_shear(block, sec_names, cal_section, c_per_bin, verbose):
    shear_shear = sec_names['shear_shear_plus']
    nbin_a, nbin_b = get_nbins(block, shear_shear)
    if c_per_bin:
        c = [block[cal_section, "c{}".format(i + 1)] for i in range(nbin_a)]
    else:
        c0 = block[cal_section, "c0"]
        c = [c0 for i in range(nbin_a)]
    calibrate_section(block, shear_shear, c, c, verbose)
    shear_shear = sec_names['shear_shear_minus']
    calibrate_section(block,shear_shear, c, c, verbose)

def execute(block, config):
    c_per_bin, sec_names, cal_section, verbose=config

    calibrate_shear_shear(
        block, sec_names, cal_section, c_per_bin, verbose)
    
    return 0 
