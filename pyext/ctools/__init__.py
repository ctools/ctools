# Import sub modules
from .tools import *
import obsutils

# Remove swig register stuff
bad_entries = [entry for entry in list(locals())
               if entry.endswith('_swigregister')]
for entry in bad_entries:
    del locals()[entry]

# Clean-up
del tools
del ctools
del gammalib
del bad_entries
del entry