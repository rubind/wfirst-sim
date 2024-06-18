#!/bin/bash
#step2_get_fom.bash
#
# Get FoM from step2 output
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
# make sure  python script for getting fom exists
if [ ! -f $WFIRST_SCRIPTS/$STEP2_GET_FOM_PYSCRIPT ]; then
  echo "ERROR : can't find step2 python script at $WFIRST_SCRIPTS/$STEP2_GET_FOM_PYSCRIPT"
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
if [ -f $STEP2_GET_FOM_LOG ]; then
  rm $STEP2_GET_FOM_LOG
fi

# make plots using STEP2_GET_FOM_PYSCRIPT  
# This also determine the Figure of Merit (FoM) of the survey, which appears as
# a comment in the output log
#
date >$STEP2_GET_FOM_LOG 
echo "python $WFIRST_SCRIPTS/$STEP2_GET_FOM_PYSCRIPT $PICKLE_FILE_STEP2" >>$STEP2_GET_FOM_LOG 
python $WFIRST_SCRIPTS/$STEP2_GET_FOM_PYSCRIPT $PICKLE_FILE_STEP2 >>$STEP2_GET_FOM_LOG 2>&1
date >>$STEP2_GET_FOM_LOG
if [ -f $STEP2_GET_FOM_LOG ] ; then
  grep Final_FoM $STEP2_GET_FOM_LOG
fi
