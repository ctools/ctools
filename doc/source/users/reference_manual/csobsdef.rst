.. _csobsdef:

csobsdef
========

Generates an observation definition file.


Synopsis
--------

The csobsdef script generates an observation definition file from
a pointing list. The pointing list is a comma-separated value (CSV)
ASCII file with header keywords in the first row followed by a list
of pointings (one pointing per row). The following header keywords
are supported (case sensitive, column order irrelevant):
    
* name     - Observation name string
* id       - Unique observation identifier string
* ra       - Right Ascension of pointing (deg)
* dec      - Declination of pointing (deg)
* lon      - Galactic longitude of pointing (deg)
* lat      - Galactic latitude of pointing (deg)
* duration - Duration of pointing (seconds)
* emin     - Lower energy limit (TeV)
* emax     - Upper energy limit (TeV)
* rad      - Radius of region of interest (deg)
* deadc    - Deadtime correction factor [0-1]
* caldb    - Calibration database
* irf      - Response function name
            
Only the pairs (ra,dec) or (lon,lat) are mandatory header keywords.
All other keywords are optional and can be specified when calling
csobsdef as user parameters. The only exception is the "duration"
keyword that will automatically be queried.
    
Here some usage examples:
    
``csobsdef``
      Creates minimal observation definition file.

``csobsdef emin=0.1 emax=100.0``
      Creates observation definition file with an energy range 100 GeV - 100 TeV.

``csobsdef rad=5``
      Creates observation definition file with a ROI radius of 5 deg.

``csobsdef caldb=prod2 irf=South_50h``
      Creates observation definition file using the ``South_50h`` IRF in the
      ``prod2`` calibration database.


General parameters
------------------

``inpnt [file]``
    Pointing definition ASCII file in comma-separated value (CSV) format.

``outobs [file]``
    Output Observation Definition XML file.

``duration [real]``
    Duration of observation (in seconds).

``(name = "") [string]``
    Observation name.

``(caldb = "") [string]``
    Calibration database.
 	 	 
``(irf = "") [string]``
    Instrumental response function.

``(emin = UNDEF) [real]``
    Lower energy limit of events (in TeV).
 	 	 
``(emax = UNDEF) [real]``
    Upper energy limit of events (in TeV).
 	 	 
``(rad = UNDEF) [real]``
    ROI radius (in degrees).

``(deadc = 0.95) [real]``
    Average deadtime correction factor.
 	 	 

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
    Specifies whether an existing Observation Definition XML file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csobsdef.log) [string]``
    Log filename.


Related tools or scripts
------------------------

None
