.. _ctbin:

ctbin
=====

Generate counts cube from event list(s).


Synopsis
--------

This tool creates a counts cube and fills it with events. A counts cube is 
a 3-dimensional data cube spanned by Right Ascension or Galactic longitude,
Declination or Galactic latitude, and energy. The energy binning may be either
linear, logarithmic, or custom defined using an input file. The events are 
either taken from a single event list file or from the event lists that are 
specified in an observation definition file. In case that multiple event 
lists are given in an observation definition file, the tool will loop over
all event lists and fill all events into a single counts cube.

ctbin generates a counts cube FITS file comprising three extensions. The
primary extension contains a 3-dimensional image that contains the counts
cube values. The next extension named ``EBOUNDS`` contains a binary table
that defines the energy boundaries of the counts cube. The last extension
named ``GTI`` contains a binary table that defines the Good Time Intervals
of all event lists that have been filled into the counts cube.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``outcube [file]``
    Output counts cube file.
 	 	 
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
 	 	 
``(usepnt = no) [boolean]``
    Use CTA pointing direction for cube centre instead of xref/yref parameters?
 	 	 
``nxpix [integer]``
    Number of cube bins in Right Ascension or Galactic longitude.
 	 	 
``nypix [integer]``
    Number of cube bins in Declination or Galactic latitude.
 	 	 
``binsz [real]``
    Cube bin size (in degrees/pixel).
 	 	  	 	 
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
    Specifies whether the counts cube should be published on VO Hub.

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

``(logfile = ctbin.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctexpcube`
:doc:`ctpsfcube`
:doc:`ctbkgcube`
