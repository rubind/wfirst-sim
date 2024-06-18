#!/bin/bash 
#step1.bash
#
# Script to run step1 of WFIRST simulation using
# Python script STEP1_simulate_survey.py

CONFIG_FILE=config.bash

# Make sure the configuration file exists and source it

if [ ! -f $CONFIG_FILE ]; then
  echo "ERROR : can't find config file $CONFIG_FILE"
  exit
fi
source $CONFIG_FILE
#
# make sure step1 python script exists
if [ ! -f $WFIRST_SCRIPTS/$STEP1_PYSCRIPT ]; then
  echo "ERROR : can't find step1 python script at $WFIRST_SCRIPTS/$STEP1_PYSCRIPT"
  exit
fi
#
# make sure paramfile exists
if [ ! -f $WFIRST_PARAMS/$PARAM_FILE ]; then
  echo "ERROR : can't find paramfile  at $WFIRST_PARAMS/$PARAM_FILE"
  exit
fi

# make the output directory and change to it
if [ ! -d $WORKDIR1 ]; then
   mkdir $WORKDIR1
fi
cd $WORKDIR1

# run step1 of the simulation using PARAM_FILE as input params. Output goes to PICKLE_FILE_STEP1
#
cp $WFIRST_PARAMS/$PARAM_FILE .
date >$STEP1_LOG 
echo "python $WFIRST_SCRIPTS/$STEP1_PYSCRIPT $PARAM_FILE $PICKLE_FILE_STEP1" >>$STEP1_LOG 
python $WFIRST_SCRIPTS/$STEP1_PYSCRIPT $PARAM_FILE $PICKLE_FILE_STEP1 >>$STEP1_LOG 2>&1
date >>$STEP1_LOG
