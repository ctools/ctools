.. _ctskymap:

ctskymap
========

Computes sky map from CTA data.


Synopsis
--------

This tool generates a sky map from CTA data.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.
 	 	 
``outmap [file]``
    Output sky map file.
 	 	 
``emin [real]``
    Minimum energy in map (in TeV).
 	 	 
``emax [real]``
    Maximum energy in map (in TeV).

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
 	 	 
``proj <AIT|AZP|CAR|MER|STG|TAN> [string]``
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

``(logfile = ctskymap.log) [string]``
    Name of log file.


Related tools
-------------

None
