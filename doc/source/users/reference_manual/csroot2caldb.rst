.. _csroot2caldb:

csroot2caldb
============

Creates a caldb entry from a ROOT performance file.


Synopsis
--------

This script creates an entry in the local calibration database from a ROOT
performance file.


General parameters
------------------

``infile [file]``
    Input ROOT performance file file.

``(outdir = ) [string]``
    Output CALDB directory.

``inst [string]``
    Instrument name.

``id [string]``
    Instrument Response Function identifier.

``version [string]``
    Instrument Response Function version.

``analysis [string]``
    Instrument Response Function analysis.

``zenith [real]``
    Zenith angle (in degrees).

``azimuth [real]``
    Azimuth angle (in degrees).

``(emin = UNDEF) [real]``
    Minimum energy threshold (TeV).

``(emax = UNDEF) [real]``
    Maximum energy threshold (TeV).

``(psftype = Gauss) [string]``
    Point spread function type.

``(split = no) [boolean]``
    Split IRF components into different files.

``(norm1d = no) [boolean]``
    Normalize on 1D histograms.

``(rebin = no) [boolean]``
    Rebin.

``(eascale = 1.0) [real]``
    Scaling factor for effective areas.

``(bgdscale = 1.0) [real]``
    Scaling factor for background rates.

``(bgdoversample = 1) [integer]``
    Spatial oversampling factor for background template.

``(bgdethres = 1000.0) [real]``
    Energy above which to replace background by power law (TeV).

``(bgdinfill = no) [boolean]``
    Infill of background template?.


Standard parameters
-------------------

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged
     
     ``chatter = 1``: only errors will be logged
     
     ``chatter = 2``: errors and actions will be logged
     
     ``chatter = 3``: report about the task execution
     
     ``chatter = 4``: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing output file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csroot2caldb.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`cscaldb`
