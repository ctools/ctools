#!/bin/csh

setenv CTOOLS ${CONDA_PREFIX}
source ${CTOOLS}/bin/ctools-init.csh
unsetenv PYTHONPATH
unsetenv LD_LIBRARY_PATH
unsetenv DYLD_LIBRARY_PATH
