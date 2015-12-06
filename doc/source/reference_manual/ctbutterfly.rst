.. _ctbutterfly:

ctbutterfly
===========

Computes butterfly diagram for a power law model.


Synopsis
--------

Calculates the confidence band of a power law model, taking into account the
covariance matrix resulting from a maximum likelihood fit. The tool derives
the envelope of all power law models whose prefactor and spectral index fall
within the error ellipse that is associated with the covariance of the 
parameters. It outputs an ASCII file in column separated value (CSV) format
containing the results. Each row in the result file corresponds to an energy.
The columns give the energy in MeV, the fitted intensity in ph/cm2/s/MeV,
and the minimum and maximum intensity that envelope the power law model 
with a given confidence (by default a confidence level of 68% is used, but 
the level can be adjusted using the ``confidence`` parameter).


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.
 	 	 
``inmodel [file]``
    Input model XML file.
 	 	 
``srcname [string]``
    Name of model component for which butterfly should be computed.
 	 	 
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
 	 	 
``outfile [file]``
    Output butterfly ASCII file name.

``(refit = no) [boolean]``
    Performs refitting of solution ignoring any provided covariance matrix.
 	 	 
``(confidence = 0.68) [real]``
    Confidence level for error computation.
    
``(matrix = "NONE") [file]``
    Input covariance matrix file (not used)

``(ebinalg = "LOG") <FILE|LIN|LOG> [string]``
    Butterfly energy binning algorithm.
 	 	 
``emin [real]``
    Minimum energy of butterfly (in TeV).
 	 	 
``emax [real]``
    Maximum energy of butterfly (in TeV).
 	 	 
``(enumbins = 100) [string]``
    Number of energy bins of butterfly.
 	 	 
``(ebinfile = "NONE") [string]``
    File for energy binning algorithm.


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

``(logfile = ctbutterfly.log) [string]``
    Name of log file.


Related tools
-------------

:ref:`ctulimit`
:ref:`cterror`
