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
setup(
    name = "NuPhyPy",
    version = "0.0.2",
    zip_safe= False,
    author = "jaromir mrazek",
    author_email = "jaromrax@gmail.com",
    description = ("A Nuclear Physics Python package "
                                   "to address many things: kinematics, srim, fresco,....."),
    license = "GPLv2",
    keywords = "nuclear physics",
    url = "http://www.spiral2.cz",
    packages=['NuPhyPy','NuPhyPy/db','NuPhyPy/Reactions',
              'NuPhyPy/Spectra''NuPhyPy/srim','NuPhyPy/db', 'NuPhyPy/fresco_binaries' ],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    ],
)
