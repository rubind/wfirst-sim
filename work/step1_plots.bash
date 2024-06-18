#!/bin/bash
#step1_plots.bash
#
# make plots from step1 output
#
CONFIG_FILE=config.bash

# Make sure the configuration file exists and source it

if [ ! -f $CONFIG_FILE ]; then
  echo "ERROR : can't find config file $CONFIG_FILE"
  exit
fi
source $CONFIG_FILE
#
#
# make sure  python script for making step1 plots exists
if [ ! -f $WFIRST_SCRIPTS/$STEP1_PLOT_PYSCRIPT ]; then
  echo "ERROR : can't find step1 python script at $WFIRST_SCRIPTS/$STEP1_PLOT_PYSCRIPT"
  exit
fi
#
# make sure step1 output exists

if [ ! -f $WORKDIR1/$PICKLE_FILE_STEP1 ]; then
   echo "error : can't find pickle output from step1"
   exit
fi

# change to the step1 output directory
cd $WORKDIR1

# make plots using STEP1_PLOT_PYSCRIPT  
#
date >$STEP1_PLOT_LOG 
echo "python $WFIRST_SCRIPTS/$STEP1_PLOT_PYSCRIPT $PICKLE_FILE_STEP1" >>$STEP1_PLOT_LOG 
python $WFIRST_SCRIPTS/$STEP1_PLOT_PYSCRIPT $PICKLE_FILE_STEP1 >>$STEP1_PLOT_LOG 2>&1
date >>$STEP1_PLOT_LOG
