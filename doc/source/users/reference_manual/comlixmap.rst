.. _comlixmap:

comlixmap
=========

Create TS map using SRCLIX algorithm.


Synopsis
--------

This script generates a Test Statistic map using the SRCLIX algorithm.


General parameters
------------------

``inobs [file]``
    Input observation definition XML file.

``inmodel [file]``
    Input model definition XML file.

``srcname [string]``
    Test source name.

``outmap [file]``
    Output Test Statistic map file.

``(max_iter = 50) [integer]``
    Maximum number of SRCLIX iterations.

``(like_accuracy = 0.05) [real]``
    Absolute accuracy of maximum likelihood value.

``(fix_spat_for_ts = no) [boolean]``
    Fix spatial parameters for TS computation?

``nxpix [integer]``
    Size of the X axis in pixels.

``nypix [integer]``
    Size of the Y axis in pixels.

``binsz [real]``
    Image scale (in degrees/pixel).

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).

``proj <AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN> [string]``
    Projection method.

``xref [real]``
    First coordinate of image center in degrees (RA or galactic l).

``yref [real]``
    Second coordinate of image center in degrees (DEC or galactic b).

``(nrunav = 3) [integer]``
    Number of bins used for running average.

``(navgr = 3) [integer]``
    Number of bins used for averaging.

``(nincl = 13) [integer]``
    Number of Phibar layers to include.

``(nexcl = 0) [integer]``
    Number of Phibar layers to exclude.


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

``(logfile = comlixfit.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`comlixfit`
:doc:`cttsmap`
