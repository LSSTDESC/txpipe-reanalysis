"""

Correlating uncorrelated priors, suggestion from Joe Zuntz, multiplying input_parameters
by the inverse cholesky decomposition of the covariance matrix of the correlated priors.
See https://scipy-cookbook.readthedocs.io/items/CorrelatedRandomSamples.html

"""


from __future__ import print_function
from builtins import range
from cosmosis.datablock import names, option_section
import sys
import numpy as np
from scipy.linalg import eigh, cholesky

warning_note_displayed = False

def setup(options):
    uncorrelated_parameters_section = options.get_string(option_section, "uncorrelated_parameters").split(' ')
    correlated_parameters_section = options.get_string(option_section, "output_parameters").split(' ')
    covariance_file = options[option_section, "covariance"]
    verbose = options.get_bool(option_section, "verbose", False)

    print("The correlated priors module will use paramter covariance from {}.".format(covariance_file))

    return uncorrelated_parameters_section, correlated_parameters_section, covariance_file

def correlate_parameters(block, uncorrelated_parameters_section, correlated_parameters_section, covariance_file):
    cov = np.loadtxt(covariance_file).reshape(5,5)
    C = cholesky(cov,lower=True)

    # Extract the uncorrelated_parameters
    print(uncorrelated_parameters_section)
    section_name = uncorrelated_parameters_section[0].split('/')[0]
    parameter_names = [uncorrelated_par.split('/')[1] for uncorrelated_par in uncorrelated_parameters_section]
    x = [block[section_name,parameter_names[i]] for i in range(len(uncorrelated_parameters_section))]

    # Convert the data to correlated random variables.
    y = np.dot(C, x)
    print('y',y)

    return y

def execute(block, config):
    uncorrelated_parameters_section, correlated_parameters_section, covariance_file = config

    correlated_parameters = correlate_parameters(block, uncorrelated_parameters_section, correlated_parameters_section, covariance_file)

    section_name = correlated_parameters_section[0].split('/')[0]
    parameter_names = [correlated_par.split('/')[1] for correlated_par in correlated_parameters_section]
    for i in range(len(parameter_names)):
        block[section_name, parameter_names[i]] = correlated_parameters[i]

    return 0
