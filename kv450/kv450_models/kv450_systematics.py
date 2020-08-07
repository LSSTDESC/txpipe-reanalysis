import numpy as np

from firecrown.ccl.core import Systematic

__all__ = ['CTermSystematic', 'ConstantCTermSystematic']

class CTermSystematic(Systematic):
    """C-term systematic.

    This systematic adds a c-term template to the galaxy_shear_xi_plus and
    galaxy_shear_xi_minus statistic.

    Parameters
    ----------
    template_filename : str
        The filename of c-term template.
    A_c : str
        The name of constant c-term parameter.

    Methods
    -------
    apply : apply the systematic to a source
    """
    def __init__(self, template_filename, A_c):
        self.template_filename = template_filename
        self.c_term = np.loadtxt(template_filename, unpack=True, usecols=[1])
        self.A_c = A_c

    def apply(self, cosmo, params, stat):
        """Apply a linear alignment systematic.

        Parameters
        ----------
        cosmo : pyccl.Cosmology
            A pyccl.Cosmology object.
        params : dict
            A dictionary mapping parameter names to their current values.
        source : a source object
            The source to which apply the shear bias.
        """
        A_c = params[self.A_c]
        stat.predicted_statistic_ += A_c**2 * self.c_term


class ConstantCTermSystematic(Systematic):
    """Constant c-term systematic.

    This systematic adds a constant c-term to the galaxy_shear_xi_plus
    statistic.

    Parameters
    ----------
    delta_c : str
        The name of constant c-term parameter.

    Methods
    -------
    apply : apply the systematic to a source
    """
    def __init__(self, delta_c):
        self.delta_c = delta_c

    def apply(self, cosmo, params, stat):
        """Apply a linear alignment systematic.

        Parameters
        ----------
        cosmo : pyccl.Cosmology
            A pyccl.Cosmology object.
        params : dict
            A dictionary mapping parameter names to their current values.
        stat : a statistics object
            The statistic to which apply the constant c-term.
        """
        delta_c = params[self.delta_c]
        stat.predicted_statistic_ += delta_c**2
