#!/bin/bash
#config_slow_test.bash
#
# configuration file for running wfirst simulation scripts.
# This is simulates a ~1-square-degree survey using a single
# Wide and a single deep tier. It requires ~3 hours of processing
# time on an 8-core machine.

# D. Rabinowitz 2016 Nov 17
export BASH_RESOURCE_FILE=`pwd`/bash_wfirst.bash
#
# Make sure bash_wfirst.bash is sourced and WFIRST_ROOT is defined
if [ -z $WFIRST_ROOT ]; then
  if [ ! -f $BASH_RESOURCE_FILE ] ; then
     echo "ERROR : can't find BASH_RESOURCE_FILE at $BASH_RESOURCE_FILE"
     exit 1
  fi
  source $BASH_RESOURCE_FILE
  if [ -z $WFIRST_ROOT ]; then
    echo "ERROR : Need to set WFIRST_ROOT in $BASH_RESOURCE_FILE"
    exit 1
  fi
fi

# initialize paths and environment variables
export WFIRST_SCRIPTS=$WFIRST_ROOT/scripts/stan_cosmo
#
#
# use the local paramfiles instead of the ones distributed with the wfirst repository
#export WFIRST_PARAMS=$WFIRST_SCRIPTS/paramfiles
export WFIRST_PARAMS=`pwd`/paramfiles
export WFIRST=$WFIRST_ROOT
export PATH=$PATH:$WFIRST_ROOT/concorde/TSP

#
# initialize name of param file

export PARAM_FILE=paramfile_quick.csv

# initial names of step1 and step2 output pickle files

export PICKLE_FILE_STEP1=pickle_step1.txt
export PICKLE_FILE_STEP2=pickle_step2.txt

export PICKLE_FILE_STEP2=fit_results.pickle

# initialize names of step1 and step2 output directories

export WORKDIR1=`pwd`/step1_output
export WORKDIR2=`pwd`/step2_output

# initialize names of output logs

export STEP1_LOG=step1.out
export STEP1_PLOT_LOG=step1_plots.out
export STEP2_LOG=step2.out
export STEP2_PLOT_LOG=step2_plots.out
export STEP2_GET_FOM_LOG=step2_get_fom.out

# intialize names of python scripts used to run the simulation

export STEP1_PYSCRIPT=STEP1_simulate_survey.py
export STEP1_PLOT_PYSCRIPT=STEP1A_plot_survey.py
export STEP2_PYSCRIPT=STEP2_UNITY.py
export STEP2_PLOT_PYSCRIPT=STEP3_plot_results.py
export STEP2_GET_FOM_PYSCRIPT=STEP3_get_fom.py

#initialize arguments for STEP2_PYSCRIPT

export STEP2_ARGS="-nrestlamb 20 -neigen 1  -gray 0.08 -nredcoeff 2 -IFCIPC 0.02 -crnl 0.003 -fund 0.005 -TTel 260 -mwnorm 0.05 -mwZP 0.003 -mwRV 0.2 -IGext 0.25 -redwave 8600.0"


