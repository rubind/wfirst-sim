#!/bin/bash
#step2.bash
#
# Run step2 scripts  using output of step1.bash to determine
# the Figure of Merit (FOM) of the survey
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
# make sure step2 python script exists
if [ ! -f $WFIRST_SCRIPTS/$STEP2_PYSCRIPT ]; then
  echo "ERROR : can't find step2 python scriprt at $WFIRST_SCRIPTS/$STEP2_PYSCRIPT"
  exit
fi
#
# make sure output from step1 exists
if [ ! -f $WORKDIR1/$PICKLE_FILE_STEP1 ]; then
  echo "ERROR: Can't find output from step1 at $WORKDIR1/$PICKLE_FILE_STEP1"
  exit
fi

# make a output directory for step2  
if [ ! -d $WORKDIR2 ]; then
   mkdir $WORKDIR2
fi

# copy the pickle output from step1 into the working directory for step2
cd $WORKDIR2
cp $WORKDIR1/$PICKLE_FILE_STEP1 .
 
 
# run step2 pyhon script
#
date >$STEP2_LOG
echo "python $WFIRST_SCRIPTS/$STEP2_PYSCRIPT -pickle $PICKE_FILE $STEP2_ARGS" >>$STEP2_LOG
python $WFIRST_SCRIPTS/$STEP2_PYSCRIPT -pickle $PICKLE_FILE_STEP1 $STEP2_ARGS >>$STEP2_LOG 2>&1
date >>$STEP2_LOG
