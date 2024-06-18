# begin .bashrc/.kshrc

if [[ -z "$SHIFTER_RUNTIME" ]]
then

# ZZ 20150915 to fix the NERSC_HOST undefined error from the batch prologue, e.g., CCM
  if [ -z "$NERSC_HOST" ]; then export NERSC_HOST=`cat /etc/clustername`;fi

# need this to define $PWD, else many scripts will fail when
# attempting "pwd" or ".." & the default login shell is csh

  pwd >/dev/null

# Host-specific settings
#export INTEL_LICENSE_FILE=28518@dmv1.nersc.gov:28518@dmv.nersc.gov

# ZZ 20150915 to fix the undefined module command error for some bash invocations
  if [ "$NERSC_HOST" = "edison" -o "$NERSC_HOST" = "cori" ]; then
    if [ -z "`type module 2>/dev/null`" ]; then 
        . /opt/modules/default/etc/modules.sh
    fi
  fi

##11/16/2015
#if [ "$NERSC_HOST" = "cori" ]; then
#  ulimit -s unlimited
#fi
#for white boxes 4/26/2016
  if [ "$NERSC_HOST" = "carl" ]; then
    . /usr/common/software/etc/modules/modules.sh
    ulimit -s unlimited
  fi

  if [ "$NERSC_HOST" = "cori" ]; then
   if [ -z "$SCRATCH" ]; then
       export SCRATCH=/global/cscratch1/sd/$USER
   fi
  fi

  $WFIRST_ROOT/bashrc_ext_system.bash

fi
# end .bashrc/.kshrc
