# Filename: ctools-init.sh
# Description: Bourne-shell flavor initialization for ctools.
#              Runs ctools-setup to generate a sh script tailored
#              specifically to this user and ctools installation,
#              then source that.
# Author/Date: Juergen Knoedlseder, IRAP, July 16, 2011
#
if [ "x$CTOOLS" = x ]; then 
  echo "ctools-init.sh: ERROR -- set CTOOLS before sourcing ctools-init.sh"
elif [ -x "$CTOOLS/bin/ctools-setup" ]; then 
  ctools_init=`$CTOOLS/bin/ctools-setup sh`
  if [ $? -eq 0 -a "x$ctools_init" != x ]; then
    if [ -f "$ctools_init" ]; then
      . $ctools_init
    fi
    rm -f $ctools_init
  fi
  unset ctools_init
else
  echo "ctools-init.sh: ERROR -- cannot execute $CTOOLS/bin/ctools-setup"
fi
