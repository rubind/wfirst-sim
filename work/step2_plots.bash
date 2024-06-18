#!/bin/bash
#step2_plots.bash
#
# make plots from step2 output
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
# make sure  python script for making step2 plots exists
if [ ! -f $WFIRST_SCRIPTS/$STEP2_PLOT_PYSCRIPT ]; then
  echo "ERROR : can't find step2 python script at $WFIRST_SCRIPTS/$STEP2_PLOT_PYSCRIPT"
  exit
fi
#
# make sure step2 output exists

if [ ! -f $WORKDIR2/$PICKLE_FILE_STEP2 ]; then
   echo "error : can't find pickle output from step2"
   exit
fi

# change to the step1 output directory
cd $WORKDIR2

# erase any previous output log if it exists
if [ -f $STEP2_PLOT_LOG ]; then
  rm $STEP2_PLOT_LOG
fi

# make plots using STEP2_PLOT_PYSCRIPT  
# This also determine the Figure of Merit (FoM) of the survey, which appears as
# a comment in the output log
#
date >$STEP2_PLOT_LOG 
echo "python $WFIRST_SCRIPTS/$STEP2_PLOT_PYSCRIPT $PICKLE_FILE_STEP2" >>$STEP2_PLOT_LOG 
python $WFIRST_SCRIPTS/$STEP2_PLOT_PYSCRIPT $PICKLE_FILE_STEP2 >>$STEP2_PLOT_LOG 2>&1
date >>$STEP2_PLOT_LOG
if [ -f $STEP2_PLOT_LOG ] ; then
  grep Final_FoM $STEP2_PLOT_LOG
fi
