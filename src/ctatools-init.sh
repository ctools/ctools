# Filename: ctatools-init.sh
# Description: Bourne-shell flavor initialization for ctatools.
#              Runs ctatools-setup to generate a sh script tailored
#              specifically to this user and ctatools installation,
#              then source that.
# Author/Date: Juergen Knoedlseder, IRAP, July 16, 2011
#
if [ "x$CTATOOLS" = x ]; then 
  echo "ctatools-init.sh: ERROR -- set CTATOOLS before sourcing ctatools-init.sh"
elif [ -x "$CTATOOLS/bin/ctatools-setup" ]; then 
  ctatools_init=`$CTATOOLS/bin/ctatools-setup sh`
  if [ $? -eq 0 -a "x$ctatools_init" != x ]; then
    if [ -f "$ctatools_init" ]; then
      . $ctatools_init
    fi
    rm -f $ctatools_init
  fi
  unset ctatools_init
else
  echo "ctatools-init.sh: ERROR -- cannot execute $CTATOOLS/bin/ctatools-setup"
fi
