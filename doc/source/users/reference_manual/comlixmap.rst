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

``(inmap = NONE) [file]``
    Input Test Statistic map file. If a TS map FITS file is specified, information
    will be extracted from the file for pixels that are contained in that map,
    avoiding recomputation. This allows for example to enlarge an existing TS map
    for which the computations will be done on the enlarged section.

``srcname [string]``
    Test source name.

``outmap [file]``
    Output Test Statistic map file.

``(max_iter = 50) [integer]``
    Maximum number of SRCLIX iterations.

``(like_accuracy = 0.05) [real]``
    Absolute accuracy of maximum likelihood value. The SRCLIX iterations terminate
    if the log-likelihood improvement is below this accuracy.

``(accept_dec = 0.0) [real]``
    Maximum accepted log-likelihood decrease. Setting this parameter to a positive
    value with allow some decrease of the log-likelihood. This may help to get the
    algorithm out of a local minimum. Use this parameter with case since it may
    lead to a solution that is not the maximum likelihood solution.

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
    Number of Chi/Psi bins used for running average (relevant for ``BGDLIXA``
    method).

``(navgr = 9) [integer]``
    Number of Chi/Psi bins used for averaging (relevant for ``BGDLIXA`` and
    ``BGDLIXE`` methods).

``(nincl = 5) [integer]``
    Number of Phibar layers to include (relevant for ``BGDLIXA`` and ``BGDLIXE``
    methods).

``(nexcl = 0) [integer]``
    Number of Phibar layers to exclude (relevant for ``BGDLIXA`` and ``BGDLIXE``
    methods).

``(phi_first = -1)  [integer]``
    First Phibar layer for likelihood fitting, starting from 0. If -1 is specified
    there is no first Phibar layer selection.

``(phi_last = -1)  [integer]``
    Last Phibar layer for likelihood fitting, starting from 0. If -1 is specified
    there is no last Phibar layer selection.


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
