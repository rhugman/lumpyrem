# Lumpyrem
Python package to interact with the LUMPREM software suite. Lumpyrem is very much a work in progress. User beware.

# Introduction
Lumpyrem is a Python package that facilitates interaction with the [LUMPREM lumped-parameter recharge model](https://s3.amazonaws.com/docs.pesthomepage.org/software/lumprem.zip) and its ancillary programs that assist in setup, and provide a linkage between LUMPREM and MODFLOW 6.
It is intended to make setting up and running LUMPREM models easier, in particular if integrated into [Flopy](https://github.com/modflowpy/flopy) workflows. 

LUMPREM requires daily inputs of rainfall and potential ET. Nonlinear relationships based on stored soil moisture are used to calculate recharge, runoff and actual ET. The parameters which govern these relationships can be adjusted by PEST. Numerous quantities of interest are recorded on LUMPREM’s output file. These include residual daily ET (which can be reformatted for input to the EVT package of MODFLOW), irrigation demand, and boundary condition head. The latter is calculated using a parameterizable, nonlinear relationship between head and stored soil moisture. Heads calculated in this way can be assigned to MODFLOW boundary conditions after appropriate re-formatting.

# Installation
Lumpyrem requires LUMPREM, LUMPREP and LR2SERIES. Executables can be downloaded from the [PESTHomepage.org](https://pesthomepage.org/software-0). It is recomended these be placed in a folder in the environment path on your local machine. Alternatively, the executables may be placed within the workspace folder.

To install Lupyrem type:
   
    pip install lumpyrem
   
   
# Documentation
Documentation is a work in progress. A Jupyter Notebook is provided with simple examples of code functionality. 

LUMPREM and utilities are extensively documented.

# To do
* deploy documentation to readthedocs
* utility to read lumprem output and convert to Feflow timeseries/.pow files
* utility to write PEST control files



