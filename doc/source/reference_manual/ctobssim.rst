ctobssim
========

Simulates CTA events.


Synopsis
--------

Generate photon events from astrophysical sources and process those photons
according to the specified instrument response functions and generate
instrumental background events.

Links
-----


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file. If provided (i.e. the
    parameter is not blank or NONE), the pointing definition and eventually the
    response information will be extracted from the input file for event
    simulation.

``inmodel [file]``
    XML file that describes the astrophysical sources and the instrumental
    background.
 	 	 
``outevents [file]``
    Output event list or observation definition XML file.
 	 	 
``(prefix = sim_events_) [string]``
    Prefix for event list in observation definition XML file.
 	 	 
``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrumental response function.
 	 	 
``(seed = 1) [integer]``
    Integer seed value to be used for Monte Carlo simulations.
 	 	 
``ra [double]``
    Right Ascension of CTA pointing (J2000, in degrees).
 	 	 
``dec [double]``
    Declination of CTA pointing (J2000, in degrees).
 	 	 
``rad [double]``
    Radius of CTA field of view (simulation cone radius) (in degrees).
 	 	 
``tmin [double]``
    Mission elapsed start time of observation (in seconds).
 	 	 
``tmax [double]``
    Mission elapsed stop time of observation (in seconds).
 	 	 
``emin [double]``
    Lower energy limit of events (in TeV).
 	 	 
``emax [double]``
    Upper energy limit of events (in TeV).
 	 	 
``(edisp = no) [boolean]``
    Apply energy dispersion?
 	 	 
``(deadc = 0.95) [double]``
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
    Specifies whether an existing output counts cube should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctobssim.log) [string]``
    Name of log file.


Related tools
-------------

None
