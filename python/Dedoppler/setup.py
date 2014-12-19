#from distutils.core import setup, Extension
import os, sys, glob
from setuptools import setup

__version__ = '0.0.1'

setup(name = 'Dedoppler',
    version = __version__,
    description = '',
    long_description = '',
    url = 'http://github.com/cs150bf/gbt_seti',
    package_dir = {'dedoppler':'dedoppler'},
    packages = ['dedoppler'],
    package_data={'dedoppler': ['drift_indexes/*.txt']},
    entry_points = {
              'console_scripts': [
                  'dedoppler = dedoppler.main:main',
                  ],
          },
    )


