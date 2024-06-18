#!/bin/bash  -l
# run.bash
# See README.run for comments

CONFIG_FILE=config.bash

# get the current directory
current_dir=`pwd`	

# Make sure the configuration file exists and source it

if [ ! -f $CONFIG_FILE ]; then
  echo "ERROR : can't find config file $CONFIG_FILE"
  exit
fi
source $CONFIG_FILE

# delete previous output from earlier runs

if [ -d $WORKDIR1 ]; then
  echo "removing previous output from step1 at $WORKDIR1"
  rm -rf $WORKDIR1
fi

if [ -d $WORKDIR2 ]; then
  echo "removing previous output from step2 at $WORKDIR2"
  rm -rf $WORKDIR2
fi

# Run step1 of the survey simulation

dt=`date`
echo "$dt Starting step1"
./step1.bash 
dt=`date`
echo "$dt Step1 completed"

# Make step1 plots
#dt=`date`
#echo "$dt Making step1 plots"
#./step1_plots.bash
#dt=`date`
#echo "$dt Step1 plots completed"

# Run step2 to get FOM
dt=`date`
echo "$dt Starting step2"
./step2.bash 
dt=`date`
echo "$dt Step2 completed"

# Get fom
dt=`date`
echo "$dt Getting FOM"
./step2_get_fom.bash
dt=`date`
echo "$dt Done getting fom"

