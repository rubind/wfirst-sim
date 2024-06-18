#!/bin/bash
# download_sim_data.bash
#
# download data files needed by wfirst-sim scripts

SOURCE_PATH="https://github.com/wfirst-sn/wfirst-sim-data "
WFIRST_SIM_DATA="wfirst-sim-data"
 
d0=`pwd`

# get the concorde source
git clone $SOURCE_PATH

if  [ ! -d $WFIRST_SIM_DATA ]; then
  echo "ERROR : could not download directory $WFIRST_SIM_DATA from $SOURCE_PATH"
  exit 1
fi

echo "$WFIRST_SIM_DATA successfully downloaded "

