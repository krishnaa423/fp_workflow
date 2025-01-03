#region modules
from setuptools import setup, find_packages
from Cython.Build import cythonize
#endregions

#region variables
#endregions

#region functions
setup(
    name='fp_workflow',
    version='2.2.0',
    description='First priciples workflow and utilities',
    author='Krishnaa Vadivel',
    author_email='krishnaa.vadivel@yale.edu',
    requires=[
        'numpy',
        'scipy',
        'ase',
        'dill',
        'pyyaml',
        'cython',
    ],
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    package_data={'fp': ['data/**/*']},
)
#endregions

#region classes
#endregions
