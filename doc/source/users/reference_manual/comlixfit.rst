.. _comlixfit:

comlixfit
=========

Fit model to data using SRCLIX algorithm.


Synopsis
--------

This script fits a model to COMPTEL data using the SRCLIX algorithm.


General parameters
------------------

``inobs [file]``
    Observation definition XML file.

``inmodel [file]``
    Input model definition XML file.

``(suffix = srclix) [string]``
    Suffix for DRB files. This suffix will be appended to each DRB file so that
    distinctive file names can be created in the same outfolder.

``(outfolder = dri) [string]``
    Output folder for DRB files.

``outobs [file]``
    Output observation definition XML file.

``outmodel [file]``
    Output model definition XML file.

``(max_iter = 50) [integer]``
    Maximum number of SRCLIX iterations.

``(like_accuracy = 0.05) [real]``
    Absolute accuracy of maximum likelihood value.

``(fix_spat_for_ts = no) [boolean]``
    Fix spatial parameters for TS computation?

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

:doc:`comlixmap`
:doc:`ctlike`
