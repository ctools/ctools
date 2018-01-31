#!/bin/sh

export CTOOLS=${CONDA_PREFIX}
source ${CTOOLS}/bin/ctools-init.sh
unset PYTHONPATH
unset LD_LIBRARY_PATH
unset DYLD_LIBRARY_PATH
