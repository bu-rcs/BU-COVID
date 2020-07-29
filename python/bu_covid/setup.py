#!/usr/bin/env python3

from setuptools import setup

# This will install the bu_covid code along with covasim 1.4.7 
# from PyPI and all necessary additional libraries.

 
setup(name='bu_covid',
      version='1.0',
      description='BU code library for use with Covasim.',
      author='BU',
      author_email='bgregor@bu.edu',
      url='https://github.com/bu-rcs/BU-COVID',
      packages=['bu_covid', ],
      package_dir={'bu_covid': 'src'},
      install_requires=['covasim==1.4.7',],
)

 
