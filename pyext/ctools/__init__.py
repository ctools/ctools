from .ctools import *
import obsutils

bad_entries = [entry for entry in list(locals())
               if entry.endswith('_swigregister')]
for entry in bad_entries:
    del locals()[entry]
del gammalib
del ctools
