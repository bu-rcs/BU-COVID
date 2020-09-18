#!/usr/bin/env python3

## ---------------------------
##
## Install the bu_covid library.
##
## Authors: Brian Gregor, Wenrui Li 
##          Boston University
##
## Date Created: 2020-07-31
##
## Email: bgregor@bu.edu
##
## ---------------------------

from setuptools import setup

# This will install the bu_covid code along with covasim 1.6.1
# from PyPI and all necessary additional libraries.

 
setup(name='bu_covid',
      version='0.2',
      description='BU code library for use with Covasim.',
      author='BU',
      author_email='bgregor@bu.edu',
      url='https://github.com/bu-rcs/BU-COVID',
      packages=['bu_covid', ],
      package_dir={'bu_covid': 'bu_covid'},
      install_requires=['covasim==1.6.1','python-igraph','networkx','tqdm'],
)

 
