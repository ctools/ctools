.. _ctmodel:

ctmodel
=======

Generate model cube.


Synopsis
--------

This tool generates a model cube based on an input model. A model cube is
a 3-dimensional cube providing the number of predicted counts for a model as 
function of Right Ascension or Galactic longitude, Declination or Galactic
latitude, and energy. The energy binning may be either linear, logarithmic,
or custom defined using an input file.

ctmodel requires on input either a counts cube, an event list or an observation
definition file.

If a counts cube is provided, ctmodel will use the definition of this cube
(such as sky coordinates and projection, number of pixels, pixel scale,
energy binning) to compute a model cube. In case that the counts cube combines
multiple observations (i.e. for a so-called "stacked cube"), an exposure cube,
a point spread function cube and a background cube have to be provided
(otherwise you may just enter ``NONE`` when the names of these files are
queried).

If an event list is provided, ctmodel will query for a counts cube to 
extract the model cube definition (parameter ``incube``). If no counts cube
is  provided (``incube=NONE``), ctmodel tools will query for cube definition
parameters.

If an observation definition file is provided, ctmodel will query for a counts
cube to extract the model cube definition, unless the observation definition
file contains a single binned observation (in that case, the counts cube of
that observation will be used to extract the model cube definition).

ctmodel generates a model cube FITS file comprising three extensions. The
primary extension contains a 3-dimensional image that contains the model 
cube values. The next extension named ``EBOUNDS`` contains a binary table
that defines the energy boundaries of the background cube. The last extension
named ``GTI`` contains a binary table that defines the Good Time Intervals
that are covered by the model cube.


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [string]``
    Input model XML file.

``incube [file]``
    Input counts cube file to extract model cube definition.

``expcube [file]``
    Input exposure cube file (only needed for stacked analysis).

``psfcube [file]``
    Input PSF cube file (only needed for stacked analysis).

``bkgcube [file]``
    Input background cube file (only needed for stacked analysis).

``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrument response function.
 	 	 
``(edisp = no) [boolean]``
    Apply energy dispersion to response computation.

``outcube [file]``
    Output model cube file.
 	 	 
``ra [real]``
    Right Ascension of pointing (J2000, in degrees).
 	 	 
``dec [real]``
    Declination of pointing (J2000, in degrees).

``rad [real]``
    Radius of field of view (in degrees).
 	 	 
``tmin [real]``
    CTA mission elapsed start time (in seconds).
 	 	 
``tmax [real]``
    CTA mission elapsed stop time (in seconds).
 	 	 
``(deadc = 0.95) [real]``
    Average deadtime correction factor.

``(ebinalg = LOG) <FILE|LIN|LOG> [string]``
    Algorithm for defining energy bins.
 	 	 
``emin [real]``
    Lower energy value for first energy bin (in TeV).
 	 	 
``emax [real]``
    Upper energy value for last energy bin (in TeV).
 	 	 
``enumbins [integer]``
    Number of energy bins.
 	 	 
``(ebinfile = NONE) [file]``
    Name of the file containing the energy bin definition.

``(usepnt = no) [boolean]``
    Use pointing instead of xref/yref parameters?
 	 	 
``nxpix [integer]``
    Size of the Right Ascension / Galactic longitude axis (in pixels).
 	 	 
``nypix [integer]``
    Size of the Declination / Galactic latitude axis (in pixels).
 	 	 
``binsz [real]``
    Pixel size (in degrees/pixel).
 	 	 
``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).
 	 	 
``xref [real]``
    Right Ascension / Galactic longitude of image centre (J2000, in degrees).
 	 	 
``yref [real]``
    Declination / Galactic latitude of image centre (J2000, in degrees).
 	 	 
``proj <AIT|AZP|CAR|MER|MOL|STG|TAN> [string]``
    Projection method.


Standard parameters
-------------------

``(publish = no) [boolean]``
    Specifies whether the model cube should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing output model cube file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctmodel.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctbin`
:doc:`ctexpcube`
:doc:`ctpsfcube`
:doc:`ctbkgcube`
