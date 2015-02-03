ctselect
========

Performs sub-selections of event data.


Synopsis
--------

Creates a new FITS file of selected events from an input event list based
on spatial, temporal, and energy cuts.


Links
-----


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file
 	 	 
``outobs [file]``
    Output event list or observation definition XML file.
 	 	 
``(prefix = "selected_") [string]``
    Prefix for output event lists in observation definition file.
 	 	 
``(usepnt = no) [boolean]``
    Use pointing instead of RA/DEC parameters?
 	 	 
``ra [double]``
    Right Ascension of acceptance cone (or ROI) centre (J2000, in degrees).
 	 	 
``dec [double]``
    Declination of acceptance cone (or ROI) centre (J2000, in degrees).
 	 	 
``rad [double]``
    Radius of acceptance cone (or ROI) centre (in degrees).
 	 	 
``tmin [double]``
    CTA mission elapsed start time defining the earliest event (in seconds).
 	 	 
``tmax [double]``
    CTA mission elapsed stop time defining the latest event (in seconds).
 	 	 
``emin [double]``
    Lower energy limit of events (in TeV).
 	 	 
``emax [double]``
    Upper energy limit of events (in TeV).
 	 	 
``(expr = "") [string]``
    Additional event selection expression (cfitsio syntax).

``(usethres = NONE) [string]``
    Energy threshold type (one of NONE, DEFAULT or USER).


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

``(logfile = ctlike.log) [string]``
    Name of log file.


Related tools
-------------

None
