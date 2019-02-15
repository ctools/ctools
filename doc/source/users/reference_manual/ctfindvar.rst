.. _ctfindvar:

ctfindvar
=========

Searches for source variability.


Synopsis
--------

This tool searches for source variability by applying the On/Off method in the
time domain. The variability of sky map pixels is computed by binning the events
as function of the photon arrival time and applying the On/Off method to each
time bin, integrating the background over all the other time bins. Only bins
with a variability significance below a given threshold, defined by the hidden
``threshold`` parameter, are included in the background estimate. Furthermore,
only time bins that exceed a minimum number of events, specified by the hidden
``minoff`` parameter, are considered for the background estimate. The significance
of each time bin is computed using Li & Ma formula.

:ref:`ctfindvar` will search a grid of sky map pixels for source variability. Using
the hidden ``smooth_kernel`` and ``smooth_rad`` parameters, the sky map may be
smoothed before the variability search. The tool will compute sky maps for a
number of time bins. The duration of each time bin is specified by the ``tinterval``
parameter, the start and stop time are defined by the ``tmin`` and ``tmax``
parameters. If the start and stop time are ``NONE``, the start and stop times will
be extracted from the input event list(s). Note that ``tinterval`` specifies the
time scale over which the variability will be searched.

The tool will determine the sky map pixel with the maximum variability and store
the variability histogram and pixel position as ``MAXSIGPIXEL`` in the output file,
specified using the ``outmap`` parameter. In addition, the tool will determine the
variability for a number of source positions. If an input model is specified
using the ``inmodel`` parameter, the source positions will be defined by the centre
positions of all gamma-ray source components in the input model. If no input
model is specified, the tool will use the ``xsrc`` and ``ysrc`` parameters to allow
the specification of a source position. The results will be stored in the output
file using either the source names found in the model definition file, or the
name ``SOURCE`` if the source position was specified by the ``xsrc`` and ``ysrc``
parameters.

All sky map pixels with a variability above the value specified by the ``threshold``
will also added be as a point source model to the ``outmodel`` file. Each source will
have a power law as the spectral component with a total flux that corresponds to
1% of the Crab nebula flux.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``(inmodel = NONE) [file]``
    Input model definition file for extracting source positions.

``(outcube = NONE) [file]``
    Filename for saving counts cube, ``NONE`` will result in no cube being saved.

``outmap [file]``
    Filename for output sky map containing the maximum significance for each
    spatial pixel. The file also contains the variability significance as function
    of time for each source of interest, and the pixel with the maximum
    significance. Furthermore, it contains the Right Ascension and Declination
    of each source of interest and the pixel with the maximum significance.

``outmodel [file]``
    Output model definition file containing sources for significant pixels.

``caldb [string]``
    Calibration database (only required for IRF background subtraction if no
    response information is provided by ``inobs``).

``irf [string]``
    Instrument response function (only required for IRF background subtraction
    if no response information is provided by ``inobs``).


Variability search parameters
-----------------------------

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).

``xsrc [real]``
    Right Ascension / Galactic longitude of source of interest (J2000, in degrees)

``ysrc [real]``
    Declination / Galactic latitude of source of interest (J2000, in degrees)

``emin [real]``
    Minimum energy (TeV) for extracting events.

``emax [real]``
    Maximum energy (TeV) for extracting events.

``(threshold = 4.5) [real]``
    Significance threshold for variability detection.

``(minoff = 0) [real]``
    Minimum number of off counts to compute the significance.

``(smooth_kernel = NONE) <GAUSSIAN|DISK|NONE> [string]``
    Kernel to be used in smoothing the counts maps before computing 
    significances.

``(smooth_rad = 0.05) [real]``
    Radius of smoothing kernel in degrees. This is either the sigma of the
    ``GAUSSIAN`` or the radius of the ``DISK``.


Time binning parameters
-----------------------

``tinterval [real]``
    Interval in seconds for each time step.

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
    Use CTA pointing direction for sky map centre instead of xref/yref parameters?

``nxpix [integer]``
    Number of cube bins in Right Ascension or Galactic longitude.

``nypix [integer]``
    Number of cube bins in Declination or Galactic latitude.

``binsz [real]``
    Sky map pixel size (in degrees/pixel).

``proj <AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN> [string]``
    Projection method.

``xref [real]``
    Right Ascension / Galactic longitude of sky map centre (J2000, in degrees).

``yref [real]``
    Declination / Galactic latitude of sky map centre (J2000, in degrees).


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
