# Filename: ctools-init.csh
# Description: C-shell flavor initialization for ctools
#              runs ctools-setup to generate a temporary csh script
#              tailored specifically to this user and ctools
#              installation, then source that.
# Author/Date: Juergen Knoedlseder, IRAP, July 16, 2011
#
if(${?CTOOLS} == 0) then 
  echo "ctools-init.csh: ERROR -- set CTOOLS before sourcing ctools-init.csh"
else if(-x "$CTOOLS/bin/ctools-setup") then 
  set ctools_init=`$CTOOLS/bin/ctools-setup csh`
  if($status == 0 && "x$ctools_init" != x) then
    if(-f "$ctools_init") then
      source $ctools_init
    endif
    \rm -f $ctools_init
  endif
  unset ctools_init
else
  echo "ctools-init.csh: ERROR -- cannot execute $CTOOLS/bin/ctools-setup"
endif
