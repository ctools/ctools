.. _ctfindvar:

ctfindvar
==============================

ToDo: Describe in a one liner what the tool is doing.


Synopsis
--------
This tool is used to search for source variability. 
The variability is computed by binning the arrival time of the photons 
and applying the on-off method to each bin, integrating the background
over all the other bins. 
A Limit is set to include a bin in the background sample, i.e., its 
significance with respect to all the others bins has to be lower than 4.5
The significance of each bin is computed using Li&Ma eq. 17.
A limit on the background sample can also be included, taking only into 
account bins that have a number of counts which is higher than a threshold
defined by the user via the "minoff" parameter.

Variability is calculated for each pixel of the sky given in the input files. 
The results of the variability is further stored in the srcsig.fits file, 
where the GTIs and the significance of each bin are stored.
Only the sources with a significance > 4.5 sigmas are stored.

:ref:`ctfindvar` takes in input some observation files.
It calculates a cube with several maps, each one with a count number integrated
over a time period defined by the tinterval parameter, for each pixel. 
The tinterval parameter is the time scale over which the variability
will be searched.

Inputs files:
- observation file xml, or event list .fits

Output files: 
- countscube.fits
    provides all the count-skymaps contained in the cube.
- peaksigmap.fits
    a skymap of the highest significance for each pixel
- srcsig.fits
    a table containing the GTIs together with the significance
    for all the sources with a significant variability detected.


The significance of a source can be displayed using the 
``show_variability_evolution.py`` script in the example forlder

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
    Minimum energy (TeV) for extracting events. A value of 0 results in no cut
    being applied.

``emax [real]``
    Maximum energy (TeV) for extracting events. A value of 0 results in no cut
    being applied.

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

``(nthreads = 0) [integer]``
    Number of parallel processes (0=use all available CPUs).

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
