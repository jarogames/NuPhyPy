import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
##########################################################
#packages=[]   #  directory with __init__.py
#modules=[]    #  single files in root dir
#
############ MANUALS 
#https://pythonhosted.org/an_example_pypi_project/setuptools.html
#
#https://pypi.python.org/pypi?%3Aaction=list_classifiers
#
################# install from git repo #############
#
# pip install git+git://github.com/jarogames/testpyrepo@master
#
############# this included all files ##############
#http://stackoverflow.com/questions/27664504/how-to-add-package-data-recursively-in-python-setup-py

import os

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

extra_files = package_files('NuPhyPy/')

############   if you want to version git and setup be the same #########
# git tag -a v$(python setup.py --version) -m '   '

setup(
    name = "NuPhyPy",
    version = "0.0.3b",
    zip_safe= False,
    author = "jaromir mrazek",
    author_email = "jaromrax@gmail.com",
    description = ("A Nuclear Physics Python package "
                                   "to address many things: kinematics, srim, fresco,....."),
    license = "GPLv2",
    keywords = "nuclear physics",
    url = "http://www.spiral2.cz",
    packages=['NuPhyPy'],
    package_data={'':extra_files},
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    ],
)
