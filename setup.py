from setuptools import setup, find_packages


setup(
    name="fp_workflow",
    version="0.1.0",
    author="Krishnaa Vadivel",
    author_email="krishnaa.vadivel@yale.edu",
    description="First principles workflow and utilities",
    url="https://gitlab.com/krishnaa42342/fp_workflow.git",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    # install_requires=[
    #     "numpy>=1.26.4",
    #     "ase>=3.23.0",
    #     "package3",
    # ],
    # entry_points={
    #     "console_scripts": [
    #         "your_command=your_module:main_function",
    #     ],
    # },
#     include_package_data=True,
#     package_data={
#         "": ["data/*.dat", "data/*.txt"],
#     },
)
