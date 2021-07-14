.. _comobsbin:

comobsbin
=========

Bin COMPTEL observations.


Synopsis
--------

This script bins COMPTEL event data into event cubes for analysis.


General parameters
------------------

``inobs [file]``
    Input observation definition XML file.

``outobs [file]``
    Output observation definition XML file.

``outfolder [string]``
    Output folder for files. Put here the location of your COMPTEL data store.

``ebinalg <FILE|LIN|LOG> [string]``
    Algorithm for defining energy bins.

``emin [real]``
    Minimum energy (MeV).

``emax [real]``
    Maximum energy (MeV).

``phase [string]``
    Phase expression in the format phasemin0-phasemax0;phasemin1-phasemax1;...

``phase0 [time]``
    Date of Phase 0 (UTC string, JD, MJD or MET in seconds).

``period [real]``
    Period (days).

``enumbins [integer]``
    Number of energy bins.

``ebinfile [string]``
    Name of the file containing the energy bin definition.

``(dchi = 1.0) [real]``
    Bin size in Chi direction (deg).

``(dpsi = 1.0) [real]``
    Bin size in Psi direction (deg).

``(dphibar = 2.0) [real]``
    Bin size in Phibar direction (deg).

``(nchi = 80) [integer]``
    Number of bins in Chi direction.

``(npsi = 80) [integer]``
    Number of bins in Psi direction.

``(nphibar = 25) [integer]``
    Number of bins in Phibar direction.


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

``(logfile = comobsbin.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

