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

``response <MODEL|SIM2|SIM3> [string]``
    Response type specifying whether ``MODEL`` IAQs or simulated IAQs provided
    in the COMPTEL calibration database should be used.

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

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).

``proj <AIT|AZP|CAR|GLS|MOL|SFL|SIN|STG|TAN> [string]``
    Projection method.

``(usepnt = yes) [boolean]``
    Use COMPTEL pointing direction for DRI centre instead of chi0/psi0 parameters?

``chi0 [real]``
    Right Ascension / Galactic longitude of DRI centre (J2000, in degrees).

``psi0 [real]``
    Declination / Galactic latitude of DRI centre (J2000, in degrees).

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

``(psdmin = 0) [integer]``
    Minimum PSD value.

``(psdmax = 110) [integer]``
    Maximum PSD value.

``(zetamin = 5.0) [real]``
    Minimum Earth horizon - Phibar (zeta) angle (deg).

``(fpmtflag = 0) [integer]``
    Handling of D2 modules with failed PMT flag (0: exclude, 1: include, 2: exclude PMT).

``(d1use = 1111111) [string]``
    D1 module usage (1: use, 0: don't use).

``(d2use = 11111111111111) [string]``
    D2 module usage (1: use, 0: don't use).


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

None
