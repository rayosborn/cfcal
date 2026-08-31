Introduction
============
CFcal is a Python package for performing calculations of crystal field (CF) 
properties of rare earth ions using the Stevens Operator formalism [K. W. H. 
Stevens, Proc. Phys. Soc. A **65**, 209 (1952)]. Once the CF parameters have 
been initialized, the CF Hamiltonian can be diagonalized to determine the 
energies and wavefunctions of all the CF levels. These can be used to determine
the magnetic susceptibility and neutron scattering spectra as a function of 
temperature.

Installing and Running
======================
CFcal requires Python 3.10 or later.

The easiest way to install CFcal is from PyPI:

```
    $ pip install cfcal
```

Alternatively, the latest development version can be installed from the CFcal
[Git repository](https://github.com/rayosborn/cfcal):

```
    $ git clone https://github.com/rayosborn/cfcal.git
    $ cd cfcal
    $ pip install .
```
