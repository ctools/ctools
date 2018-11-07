.. _ctfindvar:

ctfindvar
==============================

ToDo: Describe in a one liner what the tool is doing.


Synopsis
--------

ToDo: Desribe in detail what the tool is doing, what it takes on input and
what it produces on output. Please do not write more than 80 characters per
line since this file is also used to produce a help text for the terminal.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``(inmodel = NONE) [file]``
    Input model definition file for extracting source positions.

``histtype [string]``
    File type for storing the individual source histograms.

``prefix [string]``
    Output file prefix. The method will save two files:
    - Skymap containing the maximum significance for each spatial pixel
    - File containing the significance as a function of time interval for each 
      source

``minoff [real]``
    Minimum number of off counts to compute the significance

``emin [real]``
    Minimum energy for extracting data

``emax [real]``
    Maximum energy for extracting data

``(smoothkrnl = NONE) <GAUSSIAN|DISK|NONE> [string]``
    Kernel to be used in smoothing the counts maps before computing 
    significances.

``(smoothpar = 0.05) [real]``
    Parameter for the smoothing kernel. This is either the sigma of the GAUSSIAN
    or the radius of the DISK.


Time binning parameters
-----------------------

``tinterval [real]``
    interval in seconds for each time step.

``tmin [time]``
    Start time for interval determination (UTC string, JD, MJD or MET in seconds).
    Start times given in MET seconds are counted with respect to the time
    reference of the input observation(s).
    If ``INDEF``, ``NONE``, ``UNDEF`` or ``UNDEFINED`` is passed as value, the 
    start time will be derived from the earliest observing time.

``tmax [time]``
    Stop time for interval determination (UTC string, JD, MJD or MET in seconds).
    Stop times given in MET seconds are counted with respect to the time
    reference of the input observation(s).
    If ``INDEF``, ``NONE``, ``UNDEF`` or ``UNDEFINED`` is passed as value, the 
    stop time will be derived from the latest observing time.


Spatial binning parameters
--------------------------

``(usepnt = no) [boolean]``
    Use CTA pointing direction for cube centre instead of xref/yref parameters?

``nxpix [integer]``
    Number of cube bins in Right Ascension or Galactic longitude.

``nypix [integer]``
    Number of cube bins in Declination or Galactic latitude.

``binsz [real]``
    Cube bin size (in degrees/pixel).

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).

``proj <AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN> [string]``
    Projection method.

``xref [real]``
    Right Ascension / Galactic longitude of cube centre (J2000, in degrees).

``yref [real]``
    Declination / Galactic latitude of cube centre (J2000, in degrees).

``xsrc [real]``
    Right Ascension / Galactic longitude of source of interest (J2000, in degrees)

``ysrc [real]``
    Declination / Galactic latitude of source of interest (J2000, in degrees)


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
    Specifies whether an existing output file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = ctfindvar.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

None
