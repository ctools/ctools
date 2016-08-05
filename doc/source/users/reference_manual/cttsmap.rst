.. _cttsmap:

cttsmap
=======

Generate a Test Statistic (TS) map.


Synopsis
--------

This tool generates a Test Statistic (TS) map for a specific source model.
cttsmap works for point, radial and elliptical source models. The tool
displaces the specified source on a grid of sky directions and computes for
each direction the TS value for the model. The TS value is defined as twice
the difference between the log-likelihood obtained for the full model 
(including the specified source) and the log-likelihood obtained for a model
that excludes the specified source. The square-root of the TS values 
corresponds roughly to the pre-trial detection significance of the source (in
Gaussian sigma).

cttsmap generates a FITS file comprising a sky map of TS values followed by 
one extension per free parameter that contains sky maps of the fitted 
parameter values. To save computation time, the errors are not computed by
default. Specifying the hidden parameter ``errors=yes`` will switch on the
error computation. Accordingly, for each free parameter, the errors are also
stored in a separate sky map.


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
    
``(errors = no) [boolean]``
    Compute and store parameter errors?
 	 	 
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

``(binmin = -1) [integer]``
    First bin to compute.

``(binmax = -1) [integer]``
    Last bin to compute.

``(logLO = -1) [real]``
    LogLikelihood value of null hypothesis.
 	 	 

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

:ref:`cstsmapmerge`

