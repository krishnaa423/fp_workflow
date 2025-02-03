#region modules
from pathlib import Path
from setuptools import setup, find_packages
from Cython.Build import cythonize
#endregions

#region variables
#endregions

#region functions
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()
setup(
    name='fp_workflow',
    version='2.2.5',
    description='First priciples workflow and utilities',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Krishnaa Vadivel',
    author_email='krishnaa.vadivel@yale.edu',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    package_data={'fp': ['data/**/*']},
    requires=[
        'numpy',
        'scipy',
        'ase',
        'dill',
        'pyyaml',
        'cython',
    ],
)
#endregions

#region classes
#endregions
