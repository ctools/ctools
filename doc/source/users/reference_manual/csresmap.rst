.. _csresmap:

csresmap
========

Generates residual map.


Synopsis
--------

This scripts generates a residual map for a given model. It works for
event lists, counts cubes or observation definition files. For event
lists, parameters that define the spatial and spectral binning need to
be provided so that the script can bin the data internally. The model
is then convolved with the instrumental response function for that
binning and used for residual computation. Before residual computation,
the counts and model cubes are collapsed into maps by summing over all
energies. Three options exist then for residual computation:

* the subtraction of the model from the counts (SUB)
* the subtraction and division by the model (SUBDIV)
* the subtraction and division by the square root of the model (SUBDIVSQRT)

The residual map is written into a FITS file.  


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.

``modcube [file]``
    Input model cube file (generated with ctmodel).

``expcube [file]``
    Input exposure cube file (only needed for stacked analysis).

``psfcube [file]``
    Input PSF cube file (only needed for stacked analysis).

``bkgcube [file]``
    Input background cube file (only needed for stacked analysis).

``inmodel [file]``
    Input model XML file.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrument response function.

``(edisp = no) [boolean]``
    Apply energy dispersion?

``outmap [file]``
    Output residual counts map file.

``(ebinalg = LOG) <FILE|LIN|LOG> [string]``
    Algorithm for defining energy bins.
 	 	 
``emin [real]``
    Lower energy value for first energy bin (in TeV).
 	 	 
``emax [real]``
    Upper energy value for last energy bin (in TeV).
 	 	 
``(enumbins = 20) [integer]``
    Number of model cube energy bins.
 	 	 
``ebinfile [file]``
    Name of the file containing the energy bin definition.
 	 	 
``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).
 	 	 
``proj <AIT|AZP|CAR|MER|MOL|STG|TAN> [string]``
    Projection method.

``xref [real]``
    Right Ascension / Galactic longitude of cube centre (J2000, in degrees).
 	 	 
``yref [real]``
    Declination / Galactic latitude of cube centre (J2000, in degrees).
 	 	 
``nxpix [integer]``
    Number of cube bins in Right Ascension or Galactic longitude.
 	 	 
``nypix [integer]``
    Number of cube bins in Declination or Galactic latitude.
 	 	 
``binsz [real]``
    Cube bin size (in degrees/pixel).
 	 	 
``(algorithm = SUBDIV) <SUB|SUBDIV|SUBDIVSQRT> [string]``
    Algorithm used to generate the residual map.
 	 	 

Standard parameters
-------------------

``(publish = no) [boolean]``
    Specifies whether the residual map should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing residual map file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csresmap.log) [string]``
    Log filename.


Related tools or scripts
------------------------

None
