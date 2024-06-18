# begin .bashrc.ext

[ -e $HOME/.dbgdot ] && echo "entering .bashrc.ext"
#
# User additions to .bashrc go in this file
#
if [ -f /usr/common/usg/bin/nersc_host ]; then
  export NERSC_HOST=`/usr/common/usg/bin/nersc_host`
else
  export NERSC_HOST="unknown"
fi

if [ $NERSC_HOST == "davinci" ]
then
#  Replace the following line with personal settings for davinci
  touch /dev/null
fi

if [ $NERSC_HOST == "datatran" ]
then
#  Replace the following line with personal settings for datatran
  touch /dev/null
fi

[ -e $HOME/.dbgdot ] && echo "exiting .bashrc.ext"

# SU to project accounts

sup () {
  account=$1
  export GLOBUS_LOCATION=/usr/common/osg/gsissh/globus
  source /usr/common/osg/gsissh/globus/etc/globus-user-env.sh
  myproxy-logon -s nerscca.nersc.gov
  gsissh localhost -l ${account} -p 2222
}

export WFIRST_CORI=/global/common/cori/contrib/wfirst
# end .bashrc.ext
