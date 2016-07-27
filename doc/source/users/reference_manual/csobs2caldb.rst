.. _csobs2caldb:

csobs2caldb
===========

Creates a caldb entry from an input observation.


Synopsis
--------

This script creates an entry in the local calibration database from an
input observation container. The caldb entry can be used for simulations
with :doc:`ctobssim` later on. On default csobs2caldb uses the first observation
in the given container. This can be modified using the hidden parameter "index"
(e.g. index=3 instructs the tool to use the third observation in the container).
Note that for the moment, the script only works for unbinned observations which 
use the latest CTA IRF format (i.e. 2D/3D FITS binary tables).  

In particular this script might be useful for simulations of current IACT data.
Provided a specific observation with IRFs e.g. at certain zenith/azimuth angle, 
the user can create a caldb entry using this tool and subsequently simulate e.g. 
50 hours using this configuration. 

General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrument response function (e.g. Standard).

``(outirfs = irf_file.fits) [file]``
    Output IRF file name.

``(rootdir = ) [string]``
    Optional CALDB root directory.

``(index = 0) [integer]``
    Index of observation to be used from XML file.


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
    Specifies whether an existing output file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csobs2caldb.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`cscaldb`
