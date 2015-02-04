.. _cttsmap:

cttsmap
=======

Generate a Test Statistics map.


Synopsis
--------

This tool generates a test statistics map.


General parameters
------------------

``inobs = events [file]``
    Input event list, cube or observation definition XML file.

``inmodel = $CTOOLS/share/models/crab.xml [file]``
    Input source model XML file.

``srcname = Crab [string]``
    Name of the source in the source model XML file for which the Test
    Statistic map should be computed.

``incube = NONE [file]``
    Counts cube for background cube definition.

``caldb = dummy [string]``
    Calibration database.

``irf = cta_dummy_irf [string]``
    Response function.

``outmap = tsmap.fits [file]``
    Output Test Statistics map file.
 	 	 
``(usepnt = no) [boolean]``
    Use CTA pointing direction for map centre instead of xref/yref parameters?
 	 	 
``nxpix = 200 [integer]``
    Number of map pixels in Right Ascension or Galactic longitude.
 	 	 
``nypix = 200 [integer]``
    Number of map pixels in Declination or Galactic latitude.
 	 	 
``binsz = 0.02 [double]``
    Map pixel size (in degrees/pixel).
 	 	 
``coordsys = CEL <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).
 	 	 
``xref = 83.63 [double]``
    Right Ascension / Galactic longitude of map centre (J2000, in degrees).
 	 	 
``yref = 22.01 [double]``
    Declination / Galactic latitude of map centre (J2000, in degrees).
 	 	 
``proj = CAR <AIT|AZP|CAR|MER|MOL|STG|TAN> [string]``
    Projection method.
 	 	 

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
