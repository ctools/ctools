.. _cstssplit:

cstssplit
=======

Generate commands to split TS map computation


Synopsis
--------

This tool allows to generate a command sequence in order to compute a 
Test Statistic (TS) map using e.g. several machines. The command sequence
is stored in an ASCII file. Each line contains an individual command.
Each command will call the tool :ref:`cttsmap` with the appropriate parameters.
To split the computation, each TS map that is computed has a suffix to the filename
according to the number of the command. The resulting individual TS maps can then be
merged via :ref:`cstsmapmerge`.


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [file]``
    Input model XML file.

``srcname [string]``
    Name of source model for which the TS map should be computed.

``expcube [file]``
    Input exposure cube file (only needed for stacked analysis).

``psfcube [file]``
    Input PSF cube file (only needed for stacked analysis).

``bkgcube [file]``
    Input background cube file (only needed for stacked analysis).

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrumental response function.

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

``bins_per_job [integer]``
    Number of bins of the TS map that should be computed within one command

``compute_null [boolean]``
    Pre-computes the null hypothesis globally and passes the value to the individual commands
    If 'false', the null hypothesis will be computed in each run of :ref:`cttsmap` individually.
    
``outfile [file]``
	Output ASCII file name where the commands are written to.
 	 	 

Standard parameters
-------------------

``(publish = no) [boolean]``
    Specifies whether the TS map should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing output TS map file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = cttsmap.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:ref:`cttsmap`
:ref:`cstsmapmerge`


