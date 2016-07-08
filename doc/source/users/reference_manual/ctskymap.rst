.. _ctskymap:

ctskymap
========

Generate sky map from event list(s).


Synopsis
--------

This tool creates a sky map from either a single event list or event lists
provided in an observation definition file. The tool will loop over all event
lists that are provided and fill all events into a single sky map. Only events
within an energy interval spanned by ``emin`` and ``emax`` are considered.

ctskymap generates a FITS file comprising a sky map as primary extension.


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

``(publish = no) [boolean]``
    Specifies whether the sky map should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing output sky map file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctskymap.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

None
