# Copyright (c) 2021 Tilman Troester

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Example setup:
#
# [bin_xipm]
# file = %(REANALYSIS_PATH)s/bin_xi_pm.py
# twopoint_file = %(2PT_FILE)s
# # 9 log-spaced bins between 0.5 and 300.0
# angular_bin_edges = 0.5 1.01777898 2.0717481 4.21716333 8.58428037 17.4738002 35.56893304 72.40262468 147.37973879 300.0
# # How are the datablock bins counted: from 1 or 0?
# bin_offset = 1

import numpy as np
import scipy.integrate
import scipy.interpolate

import astropy.io.fits

from cosmosis.datablock import option_section
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    twopoint_file = options[option_section, "twopoint_file"]
    with astropy.io.fits.open(twopoint_file) as hdu:
        n_angular_bin_2pt = hdu["xim"].header["N_ANG"]
        theta_binned_2pt = hdu["xim"].data["ANG"][:n_angular_bin_2pt]
        ang_index = hdu["xim"].data.names.index("ANG")
        angle_unit = hdu["xim"].header.get('TUNIT{}'.format(ang_index + 1))
        if angle_unit.lower() == "arcmin":
            theta_binned_2pt *= np.pi/180/60

    
    angular_bin_edges = options[option_section, "angular_bin_edges"]
    angular_bin_edges *= np.pi/180/60
    n_angular_bin = len(angular_bin_edges) - 1

    if n_angular_bin != n_angular_bin_2pt:
        raise ValueError("Number of bins in 2pt file and number of bin edges "
                         "do not match.")

    xi_plus_section = "shear_xi_plus"
    xi_minus_section = "shear_xi_minus"
    interpolate = "log"

    # Bin counting from bin_0_0 or bin_1_1?
    bin_offset = options.get_int(option_section, "bin_offset", 1)

    return (angular_bin_edges, n_angular_bin,
            xi_plus_section, xi_minus_section,
            theta_binned_2pt, bin_offset, interpolate)

def execute(block, config):
    (angular_bin_edges, n_angular_bin,
     xi_plus_section, xi_minus_section,
     theta_binned_2pt, bin_offset, interpolate) = config


    for xi_section in [xi_plus_section, xi_minus_section]:
        theta = block[xi_section, "theta"]
        n_z_bin = block[xi_section, "nbin_a"]
        bins = [(i+bin_offset, j+bin_offset) for i in range(n_z_bin)
                for j in range(i+1)]

        for (bin_a, bin_b) in bins:
            bin_name = f"bin_{bin_a}_{bin_b}"
            xi = block[xi_section, bin_name]

            # When using quadrature, interpolate and integrate in log(theta)
            # or theta.
            if interpolate == "log":
                w = scipy.interpolate.InterpolatedUnivariateSpline(
                        x=np.log(theta), y=theta)
                y = scipy.interpolate.InterpolatedUnivariateSpline(
                        x=np.log(theta), y=xi)
            elif interpolate == "lin":
                w = scipy.interpolate.InterpolatedUnivariateSpline(
                        x=theta, y=theta)
                y = scipy.interpolate.InterpolatedUnivariateSpline(
                        x=theta, y=xi)

            xi_binned = np.zeros(n_angular_bin)
            for i, (l, u) in enumerate(zip(angular_bin_edges[:-1],
                                           angular_bin_edges[1:])):
                if interpolate == "log":
                    normalisation = scipy.integrate.quad(
                                        lambda log_x: np.exp(log_x)*w(log_x),
                                        a=np.log(l), b=np.log(u))[0]
                    xi_binned[i] = scipy.integrate.quad(
                                        lambda log_x: np.exp(log_x)*w(log_x)*y(log_x),
                                        a=np.log(l), b=np.log(u))[0]
                    xi_binned[i] /= normalisation
                elif interpolate == "lin":
                    normalisation = scipy.integrate.quad(
                                        lambda x: w(x), a=l, b=u)[0]
                    xi_binned[i] = scipy.integrate.quad(
                                        lambda x: w(x)*y(x),
                                        a=l, b=u)[0]
                    xi_binned[i] /= normalisation
                else:
                    mask = (l <= theta) & (theta < u)
                    x = theta[mask]
                    log_x = np.log(x)
                    w = theta[mask]
                    y = xi[mask]
                    normalisation = np.trapz(x=log_x, y=x*w)
                    xi_binned[i] = np.trapz(x=log_x, y=x*w*y)
                    xi_binned[i] /= normalisation
            
            # Overwrite the xi_pm sections
            block[xi_section, bin_name] = xi_binned
        block[xi_section, "theta"] = theta_binned_2pt

    return 0

def cleanup(config):
    pass
