.. _csphasecrv:

csphasecrv
==========

Computes phase dependent spectra for a given source.


Synopsis
--------

This script computes spectra by performing a maximum likelihood fit using
:doc:`ctlike` in a series of phase bins for pulsars. The phase bins can be
either specified in an ASCII file or as an interval divided into equally sized
phase bins. The format of the ASCII file is one row per phase bin, each
specifying the start of stop value of the phase bin, separated by a whitespace.
The phase goes from 0.0 to 1.0.

On output the script writes the fitting results into a FITS file. The script
also produces one XML file per phase bin that contains the fitting results for
that bin.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``inmodel [file]``
    Input model definition XML file.

``srcname [string]``
    Name of the periodic source in the source model XML file.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrumental response function.

``(edisp = no) [boolean]``
    Applies energy dispersion to response computation.

``outfile [file]``
    Name of the XML output file. The phase interval will be automatically
    appended to the name.

``phbinalg <FILE|LIN> [string]``
    Algorithm for defining phase bins.

``phbins [integer]``
    Number of phase bins.

``phbinfile [file]``
    File defining the phase binning.

``emin [real]``
    Lower energy limit of events (in TeV).

``emax [real]``
    Upper energy limit of events (in TeV).

``enumbins [integer]``
    Number of energy bins per phase bin (0=unbinned).

``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).

``proj <AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN> [string]``
    Projection method.

``xref [real]``
    Right Ascension / Galactic longitude of image centre (J2000, in degrees).

``yref [real]``
    Declination / Galactic latitude of image centre (J2000, in degrees).

``nxpix [integer]``
    Size of the Right Ascension / Galactic longitude axis (in pixels).

``nypix [integer]``
    Size of the Declination / Galactic latitude axis (in pixels).

``binsz [real]``
    Pixel size (in degrees/pixel).

``(statistic = DEFAULT) <DEFAULT|CSTAT|WSTAT|CHI2> [string]``
    Optimization statistic. ``DEFAULT`` uses the default statistic for all
    observations, which is ``CSTAT`` or the statistic specified in the
    observation definition XML file. ``CSTAT`` uses the C statistic for
    all observations, ``WSTAT`` uses the W statistic for all On/Off
    observations, and ``CHI2`` uses the Chi squared statistic for all
    binned or stacked observations.


Standard parameters
-------------------

``(publish = no) [boolean]``
    Specifies whether the phase curve should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing light curve output file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file
    output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csphasecrv.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctlike`
