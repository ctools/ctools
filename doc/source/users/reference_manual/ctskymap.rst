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

Optionally, the tool will subtract a background model from the sky map. The
background subtraction method can be selected using the ``bkgsubtract``
parameter. By default, no background model is subtracted (method ``NONE``).
If ``IRF`` is selected, the background template that are shipped with the
Instrument Response Functions will be used for background subtraction.

ctskymap generates a FITS file comprising a sky map as primary extension.
If a background subtraction method was selected, the FITS file will contain
the additional extensions ``BACKGRAOUND`` and ``SIGNIFICANCE`` that contain
the background map and a significance map, respectively. For the latter, the
significance of the signal is computed for each sky map pixel.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``caldb [string]``
    Calibration database (only required for IRF background subtraction if no
    response information is provided by ``inobs``).

``irf [string]``
    Instrument response function (only required for IRF background subtraction
    if no response information is provided by ``inobs``).

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
 	 	 
``proj <AIT|AZP|CAR|MER|STG|TAN> [string]``
    Projection method.

``xref [real]``
    Right Ascension / Galactic longitude of image centre (J2000, in degrees).
 	 	 
``yref [real]``
    Declination / Galactic latitude of image centre (J2000, in degrees).
 	 	 
``bkgsubtract <NONE|IRF> [string]``
    Background subtraction method.


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
    Enables debug mode. In debug mode the executable will dump any log file
    output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctskymap.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

None
