# CFcal  - Crystal Field Module
#
# Copyright (C) 2008-2026 R. Osborn, E. A. Goremychkin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#
"""
NeXpy lmfit model for crystal-field neutron-scattering spectra.

Defines CFModel, an lmfit.model.Model subclass wrapping cfcal.cfcal.CF so
that CF inelastic neutron scattering spectra can be fit interactively from
the NeXpy Fit dialog. Discovered via the nexpy.models entry-point group.
"""

import numpy as np

from lmfit.model import Model

from cfcal.cfcal import CF, REs


class CFModel(Model):
    r"""A model to describe crystal-field inelastic neutron scattering
    spectra of a trivalent rare-earth ion, with Parameters ``B20``,
    ``B22``, ``B40``, ``B42``, ``B43``, ``B44``, ``B60``, ``B62``,
    ``B63``, ``B64``, ``B66`` (Stevens crystal-field parameters, meV),
    ``Hz`` and ``Hx`` (longitudinal and transverse Zeeman terms, meV),
    ``T`` (sample temperature, K), ``sigma`` and ``gamma`` (Gaussian and
    Lorentzian pseudo-Voigt widths, meV), and ``scale`` (overall
    intensity scale factor).

    .. math::

        f(x; \{B_n^m\}, H_z, H_x, T, \sigma, \gamma, s) =
            s \sum_i I_i(\{B_n^m\}, H_z, H_x, T)\,
            V(x - E_i; \sigma, \gamma)

    where the crystal-field Hamiltonian is diagonalized using the Stevens
    Operator formalism to yield transition energies :math:`E_i` and
    thermally-weighted intensities :math:`I_i`, and :math:`V` is a
    pseudo-Voigt profile that reduces to a pure Gaussian of width
    :math:`\sigma` when :math:`\gamma = 0` (the default).

    The rare-earth ion is not a fittable Parameter but is selected at
    construction time via the ``form`` keyword (one of ``cfcal.cfcal.REs``),
    e.g. ``CFModel(form='Nd')``, since it fixes the number of crystal-field
    levels and the ion-specific constants used internally.

    ``T``, ``Hz`` and ``Hx`` describe the experimental conditions rather
    than the crystal field, so they are fixed by default and must be set
    to match the measurement; in particular ``T`` determines the Boltzmann
    weighting of the transition intensities and cannot be inferred from the
    data. ``scale`` is 1.0 by default, which reproduces the absolute cross
    section in mb/sr/meV calculated by ``CF.spectrum``; it should only be
    varied if the data are in arbitrary or detector units.

    If ``sigma`` is fixed at 0 (and ``gamma`` is 0), the model is a line
    spectrum, and NeXpy convolves it analytically with a resolution
    function (see ``lines``), placing the resolution function directly on
    each transition without a numerical convolution.
    """

    valid_forms = tuple(REs)

    def __init__(self, form='Ce', **kwargs):

        cf = CF(RE=form)

        def set_cf(B20, B22, B40, B42, B43, B44, B60, B62, B63, B64, B66,
                   Hz, Hx, T):
            cf.B20, cf.B22, cf.B40, cf.B42, cf.B43 = B20, B22, B40, B42, B43
            cf.B44, cf.B60, cf.B62, cf.B63 = B44, B60, B62, B63
            cf.B64, cf.B66 = B64, B66
            cf.Hz, cf.Hx, cf.T = Hz, Hx, T

        def cf_spectrum(x, B20=0.0, B22=0.0, B40=0.0, B42=0.0, B43=0.0,
                         B44=0.0, B60=0.0, B62=0.0, B63=0.0, B64=0.0,
                         B66=0.0, Hz=0.0, Hx=0.0, T=0.0, sigma=0.1,
                         gamma=0.0, scale=1.0):
            set_cf(B20, B22, B40, B42, B43, B44, B60, B62, B63, B64, B66,
                   Hz, Hx, T)
            if sigma > 0.0:
                return scale * cf.spectrum(eps=x, sigma=sigma,
                                            gamma=gamma if gamma else None)
            elif gamma:
                return scale * cf.spectrum(eps=x, sigma=None, gamma=gamma)
            else:
                return scale * cf.spectrum(eps=x, sigma=1e-6)

        def cf_lines(B20=0.0, B22=0.0, B40=0.0, B42=0.0, B43=0.0,
                     B44=0.0, B60=0.0, B62=0.0, B63=0.0, B64=0.0,
                     B66=0.0, Hz=0.0, Hx=0.0, T=0.0, sigma=0.1,
                     gamma=0.0, scale=1.0, **kwargs):
            if sigma > 0.0 or gamma:
                return None
            set_cf(B20, B22, B40, B42, B43, B44, B60, B62, B63, B64, B66,
                   Hz, Hx, T)
            cf.get_peaks(T)
            peaks = np.array(cf.peaks, dtype='float64').reshape(-1, 2)
            return peaks[:, 0], scale * 72.65 * cf.gJ**2 * peaks[:, 1]

        super().__init__(cf_spectrum, **kwargs)
        self.form = form
        self.lines = cf_lines

        self.set_param_hint('sigma', min=0)
        self.set_param_hint('gamma', min=0, vary=False)
        self.set_param_hint('scale', min=0)
        self.set_param_hint('T', min=0, vary=False)
        self.set_param_hint('Hz', vary=False)
        self.set_param_hint('Hx', vary=False)

    def guess(self, data, x=None, negative=False, **kwargs):
        """Estimate initial model parameter values from data.

        Only the peak width is estimated, using the same convention as
        ``CF.spectrum``. The crystal-field parameters cannot be guessed
        from the data and are left at zero, and ``scale`` is left at 1.0
        so that the model is calculated in absolute units of mb/sr/meV.
        """
        sigma = 0.01 * (np.max(x) - np.min(x)) if x is not None else 0.1
        return self.make_params(sigma=sigma)
