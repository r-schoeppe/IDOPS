from distutils.core import setup

from setuptools import find_packages
import glob


setup(
    name='idops',
    version='0.2.2',

    url='https://gitlab.gwdg.de/microbiology-bioinformatics-goe/idops',
    license='GPLv3',
    author='Raphael Schöppe',
    author_email='gitlab+microbiology-bioinformatics-goe-idops-12196-issue-@gwdg.de',
    description='IDentification Of Pesticidal Sequences:',

    packages=find_packages(include="idops/*"),
    include_package_data=True,

    entry_points={
      'console_scripts': ['idops=idops.idops:main']
      },
    options={'build_scripts': {'executable': '/usr/bin/env python3.8'}},

    scripts=['scripts/Easyfig_idops.py'],

    requires=[
        "biopython",
        "pandas",
        "matplotlib",
        "numpy",
        "scipy",
    ],
)
