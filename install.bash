#!/bin/bash
#
# install everything needed to run WFIRST scripts
# 
#
#
# set WFIRST_ROOT to current directory
# All code will use this environment variable to local Python source
# and environment settings

export WFIRST_ROOT=`pwd`
export HOME=$WFIRST_ROOT

# creat bash config file with required paths and
# environment vatiables. All bash scripts will
# source this config file to setup paths and
# python environment options

BASH_RESOURCE_FILE=$WFIRST_ROOT/work/bash_wfirst.bash
echo "#bash resources for wfirst simulation scripts" > $BASH_RESOURCE_FILE

# set the location for the conda environment to be used for Python 
WFIRST_CONDA_ENV=$WFIRST_ROOT/wfirst-env


# append a line to BASH_RESOURCE_FILE defining WFIRST_ROOT.
# Also append line source the bashrc_system.bash file. This is
# a copy of the .bashrc file supplied by NERSC for system users.
# The idea is to make the code run indepedently of any specific
# users. 
 
echo "export WFIRST_ROOT=$WFIRST_ROOT"  >> $BASH_RESOURCE_FILE
echo "source $WFIRST_ROOT/bashrc_system.bash" >> $BASH_RESOURCE_FILE

# Create the conda environment for running Python scripts

module load python/2.7-anaconda
/usr/bin/yes | conda create --no-deps -p $WFIRST_CONDA_ENV numpy
/usr/bin/yes | conda install --no-deps -p $WFIRST_CONDA_ENV -c OpenAstronomy sncosmo
/usr/bin/yes | conda install --no-deps -p $WFIRST_CONDA_ENV -c OpenAstronomy sep
/usr/bin/yes | conda install --no-deps -p $WFIRST_CONDA_ENV -c OpenAstronomy pystan
/usr/bin/yes | conda install --no-deps -p $WFIRST_CONDA_ENV -c OpenAstronomy corner
/usr/bin/yes | conda install --no-deps -p $WFIRST_CONDA_ENV -c CEFCA pyfits
/usr/bin/yes | conda install --no-deps -p $WFIRST_CONDA_ENV -c astropy extinction

# Append lines to BASH_RESOURCE_FILE directing any bash script
# that sources that file to activate the wfirst conda environment

echo "module load python/2.7-anaconda" >> $BASH_RESOURCE_FILE
echo "source deactivate root" >> $BASH_RESOURCE_FILE
echo "source activate $WFIRST_CONDA_ENV" >> $BASH_RESOURCE_FILE

# download, install, and compile the concorde source code
echo "installing concorde"
cd $WFIRST_ROOT
./install_concorde.bash
#
echo ""
echo "Done with installation."
echo ""
echo "*** Be sure to source $BASH_RESOURCE_FILE before using the code ***"

