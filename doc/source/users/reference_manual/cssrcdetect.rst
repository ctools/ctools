.. _cssrcdetect:

cssrcdetect
===========

Detects sources in sky map and writes results into a model definition XML file.


Synopsis
--------

This script detects sources in a sky map using a peak detection method and
writes the detected sources into a model definition XML file. In addition,
the script produces also a DS9 region file of the detected sources that can
be overlayed over a sky map displayed with DS9.

For the moment, the script writes out all sources as point sources. Optionally,
a background model can be added to the model definition XML file so that the
output file can be used directly as input of a :doc:`ctlike` run.


General parameters
------------------

``inmap [file]``
    Input sky map file.

``outmodel [file]``
    Output model definition XML file.

``outds9file [file]``
    Output DS9 region file.

``srcmodel <POINT> [string]``
    Source model type. For the moment, all detected sources will be added as
    point sources with power law spectral shapes to the model definition
    XML file.

``bkgmodel <NONE|IRF|AEFF|CUBE> [string]``
    Background model type. Using ``NONE`` no background model will be added
    to the model definition XML file. Using ``IRF`` a background model based
    on the template information in the Instrument Response Function will be
    added. Using ``AEFF`` a background model based on the shape of the effective
    area will be added. Using ``CUBE`` a background model for stacked analysis
    will be added.

``threshold [real]``
    Detection threshold (Gaussian sigma). Only sources above this detection
    threshold will be added to the model definiton XML file.

``(maxsrcs = 20) [integer]``
    Maximum number of sources that should be detected.

``(exclrad = 0.2) [real]``
    Radius around a detected source that is excluded from further source
    detection (degrees).

``(fit_pos = yes) [boolean]``
    Enable fitting of source positions in model definiton XML file?

``(fit_shape = yes) [boolean]``
    Enable fitting of source shapes in model definiton XML file?

    
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
    Specifies whether an existing output model file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file
    output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = cssrcdetect.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

None
