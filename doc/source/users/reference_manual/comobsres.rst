.. _comobsres:

comobsres
=========

Generate residuals of COMPTEL observations.


Synopsis
--------

This script generates residual maps for COMPTEL observations.


General parameters
------------------

``inobs [file]``
    Input observation definition XML file.

``inmodel [file]``
    Input model definition XML file.

``outmap [file]``
    Output residual map.

``(outfolder = resmap) [string]``
    Output folder for DRI and PNG files.

``algorithm <SUB|SUBDIV|SUBDIVSQRT|SIGNIFICANCE> [string]``
    Residual map computation algorithm.

``(armmin = -3.0) [real]``
    Minimum ARM (deg).

``(armmax = 3.0) [real]``
    Maximum ARM (deg).

``(margin = 30.0) [real]``
    Sky map margin (deg). This margin is added to the data cube size so that
    also source positions outside the range spanned by the data cube are considered.

``(grouping = 1) [integer]``
    Number of Chi/Psi bins to group for residual computation.

``(dri = no) [bool]``
    Compute DRI residuals?

``(png = no) [bool]``
    Generate residual histogram png files?


Standard parameters
-------------------

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing energy boundaries output file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = comobsres.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

