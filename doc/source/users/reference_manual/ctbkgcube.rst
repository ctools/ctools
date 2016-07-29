.. _ctbkgcube:

ctbkgcube
=========

Generate background cube for a stacked analysis.


Synopsis
--------

This tool generates a background cube for a stacked analysis based on an
input model. A background cube is a 3-dimensional cube spanned by Right
Ascension or Galactic longitude, Declination or Galactic latitude, and energy.
The energy binning may be either linear, logarithmic, or custom defined using
an input file. The input model is used to predict the expected number of
background counts in each background cube bin.

ctbkgcube requires on input the event list or observation definition file 
that has been used in the generation of the counts cube using :doc:`ctbin`.

To assure consistency between an existing counts cube and the 
corresponding background cube, it is recommended to specify the counts 
cube as ``incube`` parameter. This will instruct ctbkgcube to extract the 
background cube definition (such as sky coordinates and projection, number 
of pixels, pixel scale, energy binning) from the counts cube.

Note that there should at least be 25 bins per energy decade in the background
cube to assure that the energy dependence of the background rate is sufficiently
well sampled (in particular at low energies).

ctbkgcube generates a background cube FITS file comprising two extensions.
The primary extension contains a 3-dimensional image that contains the 
background cube values. The next extension named ``ENERGIES`` contains a
binary table that defines the energies of the background cube.

ctbkgcube generates also an output model XML file that can serve as input 
for a maximum likelihood analysis. The output model XML file is a copy of
the input model XML file where the input background model has been replaced
by a background model of type ``CTACubeBackground``. The ``CTACubeBackground``
background model instructs any tool analysing binned data to extract 
background information from a background cube. The ``CTACubeBackground``
model has a spectral component that can be adjusted in a maximum 
likelihood fit to accomodate for uncertainties in the prediction of the 
energy dependence of the background rate. ctbkgcube will use a power law
as spectral component, but you can replace this by any component of your
choice.



General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``inmodel [file]``
    Input model XML file.

``incube [file]``
    Input counts cube file to extract background cube definition.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrument response function.

``outcube [file]``
    Output background cube file.

``outmodel [file]``
    Output model XML file.

``ebinalg <FILE|LIN|LOG> [string]``
    Algorithm for defining energy bins.
 	 	 
``emin [real]``
    Lower energy value for first energy bin (in TeV).
 	 	 
``emax [real]``
    Upper energy value for last energy bin (in TeV).
 	 	 
``enumbins [integer]``
    Number of energy bins.
 	 	 
``ebinfile [file]``
    Name of the file containing the energy bin definition.

``(addbounds = no) [boolean]``
    Add energies to the background cube at the observation energy boundaries?

``(usepnt = no) [boolean]``
    Use CTA pointing direction for cube centre instead of xref/yref parameters?
 	 	 
``nxpix [integer]``
    Number of cube bins in Right Ascension or Galactic longitude.
 	 	 
``nypix [integer]``
    Number of cube bins in Declination or Galactic latitude.
 	 	 
``binsz [real]``
    Cube bin size (in degrees/pixel).
 	 	 
``coordsysL <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).
 	 	 
``xref [real]``
    Right Ascension / Galactic longitude of cube centre (J2000, in degrees).
 	 	 
``yref [real]``
    Declination / Galactic latitude of cube centre (J2000, in degrees).
 	 	 
``proj <AIT|AZP|CAR|MER|MOL|STG|TAN> [string]``
    Projection method.
 	 	 

Standard parameters
-------------------

``(publish = no) [boolean]``
    Specifies whether the background cube should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing output background cube file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctbkgcube.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctbin`
