from __future__ import print_function
from builtins import range
from cosmosis.datablock import option_section, names
import numpy as np

def setup(options):
    perbin = options.get_bool(option_section, "perbin", False)
    psf_template_file = options.get_string(option_section, "psf_template_file")
    alpha_beta = options.get_bool(option_section, "alpha_beta", False)
    do_xi_plus = options.get_bool(option_section, "do_xi_plus", True)
    do_xi_minus = options.get_bool(option_section, "do_xi_minus", True)
    cal_section = options.get_string(
        option_section, "cal_section", default=names.psf_errors)

    suffix = options.get_string(option_section, "suffix", "")

    print()
    print(" - shear-shear.")

    if suffix:
        suffix = "_" + suffix

    sec_names = {
        "shear_shear_plus": "shear_xi_plus" + suffix,
        "shear_shear_minus": "shear_xi_minus" + suffix,
    }   

    return perbin, psf_template_file, do_xi_plus, do_xi_minus, alpha_beta, cal_section


def execute(block, config):
    perbin, psf_template_file, alpha_beta, cal_section, sec_names = config

    if do_xi_plus:
        shear_shear = sec_names['shear_shear_plus']

        if not alpha_beta:
            psf_template = np.loadtxt(psf_template_file, unpack=True, usecols=[1])
        else:
            pp = np.loadtxt(psf_template_file, unpack=True, usecols=1, skiprows=5)
            pq = np.loadtxt(psf_template_file, unpack=True, usecols=2, skiprows=5)
            qq = np.loadtxt(psf_template_file, unpack=True, usecols=3, skiprows=5)

        nbin_shear = block[shear_shear, 'nbin_a']

        if perbin:
            if not alpha_beta:
                A = [block[parameters, "A{}".format(i + 1)]
                     for i in range(nbin_shear)]
            else:
                alpha = [block[parameters, "alpha{}".format(i + 1)]
                     for i in range(nbin_shear)]
                beta = [block[parameters, "beta{}".format(i + 1)]
                     for i in range(nbin_shear)]
        else:
            if not alpha_beta:
                A0 = block[cal_section, "A0"]
                A = [A0 for i in range(nbin_b)]
            else:
                alpha0 = block[cal_section, "alpha0"]
                beta0 = block[cal_section, "beta0"]
                alpha = [alpha0 for i in range(nbin_shear)]
                beta = [beta0 for i in range(nbin_shear)]

        for i in range(nbin_shear):
            for j in range(i + 1):
                bin_ij = 'bin_{0}_{1}'.format(i + 1, j + 1)
                if not alpha_beta:
                    block[shear_shear, bin_ij] +=  A[i] * A[j] * psf_template
                else:
                    block[shear_shear, bin_ij] += (alpha[i] * alpha[j])*pp+(alpha[i] * beta[i])*pq+(beta[i]*beta[j])*qq
    if do_xi_minus:
        shear_shear = sec_names['shear_shear_minus']

        if not alpha_beta:
            psf_template = np.loadtxt(psf_template_file, unpack=True, usecols=[1])
        else:
            pp = np.loadtxt(psf_template_file, unpack=True, usecols=1, skiprows=5)
            pq = np.loadtxt(psf_template_file, unpack=True, usecols=2, skiprows=5)
            qq = np.loadtxt(psf_template_file, unpack=True, usecols=3, skiprows=5)

        nbin_shear = block[shear_shear, 'nbin_a']

        if perbin:
            if not alpha_beta:
                A = [block[parameters, "A{}".format(i + 1)]
                     for i in range(nbin_shear)]
            else:
                alpha = [block[parameters, "alpha{}".format(i + 1)]
                     for i in range(nbin_shear)]
                beta = [block[parameters, "beta{}".format(i + 1)]
                     for i in range(nbin_shear)]
        else:
            if not alpha_beta:
                A0 = block[cal_section, "A0"]
                A = [A0 for i in range(nbin_b)]
            else:
                alpha0 = block[cal_section, "alpha0"]
                beta0 = block[cal_section, "beta0"]
                alpha = [alpha0 for i in range(nbin_shear)]
                beta = [beta0 for i in range(nbin_shear)]

        for i in range(nbin_shear):
            for j in range(i + 1):
                bin_ij = 'bin_{0}_{1}'.format(i + 1, j + 1)
                if not alpha_beta:
                    block[shear_shear, bin_ij] +=  A[i] * A[j] * psf_template
                else:
                    block[shear_shear, bin_ij] += (alpha[i] * alpha[j])*pp+(alpha[i] * beta[i])*pq+(beta[i]*beta[j])*qq


    return 0

