# Filename: ctatools-init.csh
# Description: C-shell flavor initialization for ctatools
#              runs ctatools-setup to generate a temporary csh script
#              tailored specifically to this user and ctatools
#              installation, then source that.
# Author/Date: Juergen Knoedlseder, IRAP, July 16, 2011
#
if(${?CTATOOLS} == 0) then 
  echo "ctatools-init.csh: ERROR -- set CTATOOLS before sourcing ctatools-init.csh"
else if(-x "$CTATOOLS/bin/ctatools-setup") then 
  set ctatools_init=`$CTATOOLS/bin/ctatools-setup csh`
  if($status == 0 && "x$ctatools_init" != x) then
    if(-f "$ctatools_init") then
      source $ctatools_init
    endif
    \rm -f $ctatools_init
  endif
  unset ctatools_init
else
  echo "ctatools-init.csh: ERROR -- cannot execute $CTATOOLS/bin/ctatools-setup"
endif
