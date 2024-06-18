#!/bin/bash
# queue_job.bash
#
# Copy the necessary scripts and config files to a scratch work area and
# start a job to run both step1 and step2 of the WFIRST simulation
#

# first source the config file to determine
# computing environment

source config.bash

# set names of sub-scripts and sub-directories need to run the primary script
SCRIPT_LIST=`ls -d *.bash paramfiles`

# set name of primary script to run
WFIRST_SCRIPT=run.bash

# get name of current directory. The above scripts will be copied from here to a work
# area on the scratch disk. That is where the job will run in batch mode

current_directory=`pwd`

# get the current date in form "yyyymmddhhmmss". THis will be used to tag the work directory
# and queue script name name in the scratch area

dt=`date +"%Y%m%d%H%M%S"`

# start making a script to queue a job on cori to run WFIRST_SCRIPT

# set the name of the script to queue the job
JOB_FILE="queue_job.$dt"

# set the queue partition
PARTITION=shared

# set the scratch area to run the job
SCRATCH_AREA=$CSCRATCH

# set the JOB time limit
RUN_TIME="00:15:00"

# set the repository for charging the run tim
REPOSITORY=m1187

# set the  job name
JOB_NAME="wfirst_sim"

# make a directory at $SCRATCH_AREA with the current time appended to the name (yyyymmddhhmmss)

RUN_DIR=wfirst_sim.$dt

# make sure scratch area exists
if [ ! -d $SCRATCH_AREA ]; then
  echo "ERROR : can't find scratch area at $SCRATCH_AREA"
  exit 1
fi

# copy the bash resource file to SCRATCH_AREA
echo "copying bash and parameter files to $SCRATCH_AREA"
cp $BASH_RESOURCE_FILE $SCRATCH_AREA

# create the RUN_DIR at the SCRATCH_AREA. Make sure it doesn't already exist first.

if [ -d $SCRATCH_AREA/$RUN_DIR ]; then
  echo "ERROR : run directory $RUN_DIR already exists at $SCRATCH_AREA/$RUN_DIR"
  exit 1
fi

mkdir $SCRATCH_AREA/$RUN_DIR

# copy the necessary scripts to the RUNDIR
echo "copying $SCRIPT_LIST to $SCRATCH_AREA/$RUN_DIR"
cp -r $SCRIPT_LIST $SCRATCH_AREA/$RUN_DIR

# change to the RUN_DIR and start creating the job script
#
cd $SCRATCH_AREA/$RUN_DIR

echo '#!/bin/bash' >$JOB_FILE
echo "#SBATCH --partition=$PARTITION" >>$JOB_FILE
echo "#SBATCH -C haswell" >>$JOB_FILE
echo "#SBATCH --time=$RUN_TIME" >>$JOB_FILE
echo "#SBATCH --job-name=$JOB_NAME" >>$JOB_FILE
echo "#SBATCH --mem=15000" >>$JOB_FILE
echo "#SBATCH -A $REPOSITORY" >>$JOB_FILE
echo "cd $SCRATCH_AREA/$RUN_DIR" >> $JOB_FILE
echo "srun -n 1 -c 8 $WFIRST_SCRIPT" >>$JOB_FILE

cd $SCRATCH_AREA/$RUN_DIR
sbatch $JOB_FILE


