# setup.py  
# This file is part of gnssr4river.
# Author Lubin Roineau (lubin.roineau@ensg.eu), 2022

import setuptools
from setuptools import find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="gnssr4river",
    author="Lubin Roineau",
    author_email="lubin.roineau@ensg.eu",
    version="1.0.0",
    description="Python librairy to help perform gnss reflectometry",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lroineau/gnssr4river",
    packages=find_packages("."),
    package_dir={"":"."},
    install_requires=['numpy','math','pandas','shapely','matplotlib','astropy','os','geopandas','unlzw3','pathlib','GDAL','Shapely','io','cartopy','urllib'],
    classifiers=["Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering",
        "Intended Audience :: Science/Research",
        "Development Status :: Beta"]
    
)