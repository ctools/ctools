.. _ctmapcube:

ctmapcube
=========

Generate a map cube from an input model.


Synopsis
--------

This tool generates a map cube from all sky model components contained in
an input model. A map cube is a series of sky maps, each spanned by Right
Ascension or Galactic longitude and Declination or Galactic latitude. Each
map corresponds to a specific photon energy. The energy binning of the map
cube may be either linear, logarithmic, or custom defined using an input
file.

ctmapcube generates a map cube FITS file comprising two extensions. The
primary extension contains a 3-dimensional image containing the map cube
values. A second extension named ``ENERGIES`` contains a binary table that
defines the energy values of each map in the cube.


General parameters
------------------

``inmodel [file]``
    Input model XML file.

``outcube [file]``
    Output map cube file.

``(ptsrcsig = 1.0) [real]``
    Sigma of Gaussian to be used for point sources (in arcmin).

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
 	 	 
``nxpix [integer]``
    Number of cube bins in Right Ascension or Galactic longitude.
 	 	 
``nypix [integer]``
    Number of cube bins in Declination or Galactic latitude.
 	 	 
``binsz [real]``
    Image scale (in degrees/pixel).
 	 	  	 	 
``coordsys <CEL|GAL> [string]``
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
    Specifies whether the map cube should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing output counts cube file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctmapcube.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

None
