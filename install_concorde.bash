#!/bin/bash
# install_concorde.bash
#
# download and install the concorde travelling-salesman-problem solver

SOURCE_PATH="http://www.math.uwaterloo.ca/tsp/concorde/downloads/codes/src"
SOURCE_FILE="co031219.tgz"

d0=`pwd`

# get the concorde source
wget "$SOURCE_PATH/$SOURCE_FILE"


if  [ ! -f $SOURCE_FILE ]; then
  echo "ERROR : could not download concorde source $SOURCE_FILE from $SOURCE_PATH"
  exit 1
fi

tar xvfz  $SOURCE_FILE
cd concorde
./configure
make

cd $d0
if [ ! -f concorde/TSP/concorde ]; then
  echo "ERROR : concorde executable did not compile"
  exit 1
fi

rm -f $SOURCE_FILE

echo "concorde successfully installed"

