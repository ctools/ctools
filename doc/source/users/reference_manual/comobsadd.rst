.. _comobsadd:

comobsadd
=========

Combine observations.


Synopsis
--------

This script combines data for multiple viewing period into a single set of
observations.


General parameters
------------------

``inobs [file]``
    Input observation definition XML file.

``inmodel [file]``
    Input model definition XML file.

``(outfolder = dri) [string]``
    Output folder for files.

``outobs [file]``
    Output observation definition XML file.

``(prefix = com) [string]``
    Filename prefix.

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).

``proj <AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN> [string]``
    Projection method.

``ra [real]``
    Right Ascension of DRI centre (deg).

``dec [real]``
    Declination of DRI centre (deg).

``glon [real]``
    Galactic longitude of DRI centre (deg).

``glat [real]``
    Galactic latitude of DRI centre (deg).

``nchi [integer]``
    Number of bins in Chi direction.

``npsi [integer]``
    Number of bins in Psi direction.

``(dchi = 1.0) [real]``
    Bin size in Chi direction (deg).

``(dpsi = 1.0) [real]``
    Bin size in Psi direction (deg).


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

``(logfile = comobsadd.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

