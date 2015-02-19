from setuptools import setup, find_packages
import os
import re

here = os.path.abspath(os.path.dirname(__file__))

# Read the version number from a source file.
# Code taken from pip's setup.py
def find_version(*file_paths):
    with open(os.path.join(here, *file_paths), 'r') as f:
        version_file = f.read()

    # The version line must have the form
    # __version__ = 'ver'
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


# Get the long description from the relevant file
with open('DESCRIPTION.rst') as f:
    long_description = f.read()

setup(
    name="CFlib",
    version=find_version('src', 'cflib', '__init__.py'),
    description="Python calculator for crystal fields",
    long_description=long_description,

    url='http://github.com/rayosborn/CFlib',

    author='Ray Osborn',
    author_email='rayosborn@mac.com',

    license='BSD',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],

    keywords='rare earth crystal fields',

    package_dir = {'': 'src'},
    packages=find_packages('src'),
    package_data={
        'CFlib': ['*.db'],
    },

#    entry_points={
#        'console_scripts': [
#            'sample=sample:main',
#        ],
#    },

    requires=['scipy', 'numpy'],
)
