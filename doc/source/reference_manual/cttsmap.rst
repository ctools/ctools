.. _cttsmap:

cttsmap
=======

Generate a Test Statistic map.


Synopsis
--------

This tool generates a Test Statistic map.


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [file]``
    Input source model XML file.

``srcname [string]``
    Name of the source in the source model XML file for which the Test
    Statistic map should be computed.

``expcube [file]``
    Input exposure cube file (only needed for stacked analysis).

``psfcube [file]``
    Input PSF cube file (only needed for stacked analysis).

``bkgcube [file]``
    Input background cube file (only needed for stacked analysis).

``caldb [string]``
    Calibration database.

``irf [string]``
    Response function.

``(edisp = no) [boolean]``
    Applies energy dispersion to response computation.

``outmap [file]``
    Output Test Statistic map file.
 	 	 
``(usepnt = no) [boolean]``
    Use CTA pointing direction for map centre instead of xref/yref parameters?
 	 	 
``nxpix [integer]``
    Number of map pixels in Right Ascension or Galactic longitude.
 	 	 
``nypix [integer]``
    Number of map pixels in Declination or Galactic latitude.
 	 	 
``binsz [real]``
    Map pixel size (in degrees/pixel).
 	 	 
``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).
 	 	 
``xref [real]``
    Right Ascension / Galactic longitude of map centre (J2000, in degrees).
 	 	 
``yref [real]``
    Declination / Galactic latitude of map centre (J2000, in degrees).
 	 	 
``proj <AIT|AZP|CAR|MER|MOL|STG|TAN> [string]``
    Projection method.

``(binmin = -1) [integer]``
    First bin to compute.

``(binmax = -1) [integer]``
    Last bin to compute.

``(logLO = -1) [real]``
    LogLikelihood value of null hypothesis.
 	 	 

Standard parameters
-------------------

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing output counts cube should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = cttsmap.log) [string]``
    Name of log file.


Related tools
-------------

None
