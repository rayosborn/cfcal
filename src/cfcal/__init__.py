__package_name__ = u'cfcal'
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

__documentation_author__ = u'Ray Osborn'
__documentation_copyright__ = u'2017, Ray Osborn'

__license__ = u'BSD'
__author_name__ = u'Ray Osborn'
__author_email__ = u'rayosborn@mac.com'
__author__ = __author_name__ + u' <' + __author_email__ + u'>'

__description__ = u'CFcal: Python calculator for crystal fields'
__long_description__ = \
u"""
This is a Python package for performing calculations of crystal field (CF) 
properties of rare earth ions using the Stevens Operator formalism [K. W. H. 
Stevens, Proc. Phys. Soc. A 65, 209 (1952)]. Once the CF parameters have 
been initialized, the CF Hamiltonian can be diagonalized to determine the 
energies and wavefunctions of all the CF levels. These can be used to determine
the magnetic susceptibility and neutron scattering spectra as a function of 
temperature.
"""
