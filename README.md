Introduction
============
CFlib is a Python package for performing calculations of crystal field (CF) 
properties of rare earth ions using the Stevens Operator formalism [K. W. H. 
Stevens, Proc. Phys. Soc. A **65**, 209 (1952)]. Once the CF parameters have 
been initialized, the CF Hamiltonian can be diagonalized to determine the 
energies and wavefunctions of all the CF levels. These can be used to determine
the magnetic susceptibility and neutron scattering spectra as a function of 
temperature.

Installing and Running
======================
The latest version of CFlib can be downloaded from the CFlib [Git 
repository](https://github.com/rayosborn/cflib).

From within the main directory, you can install CFlib using:

```
    $ python setup.py install
```

or

```
     $ pip install ./
```

To install in an alternate location, use:

```
    $ python setup.py install --prefix=/path/to/installation/dir
```
