from setuptools import setup, find_packages
import os
import re
import versioneer

# Get the long description from the relevant file
with open('DESCRIPTION.rst') as f:
    long_description = f.read()

setup(
    name="cfcal",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Python calculator for crystal fields",
    long_description=long_description,

    url='http://github.com/rayosborn/cfcal',

    author='Ray Osborn',
    author_email='rayosborn@mac.com',

    license='BSD',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    keywords='rare earth crystal fields',

    package_dir = {'': 'src'},
    packages=find_packages('src'),
    package_data={
        'cfcal': ['*.ini'],
    },

#    entry_points={
#        'console_scripts': [
#            'sample=sample:main',
#        ],
#    },

    requires=['scipy', 'numpy'],
)
