.. _ctbutterfly:

ctbutterfly
===========

Computes butterfly diagram for a power law model.


Synopsis
--------

This tool calculates a butterfly diagram for a specific source with power law
spectral model. The butterfly diagram is the envelope of all power law models
that are within a given confidence limit compatible with the data. By default
a confidence level of 68% is used, but the level can be adjusted using the 
hidden ``confidence`` parameter. ctbutterfly computes this envelope by
evaluating for each energy the minimum and maximum intensity of all power law
models that fall within the error ellipse of the prefactor and index parameters.
The error ellipse is derived from the covariance matrix of a maximum likelihood
fit.

ctbutterfly assumes that the input model (parameter ``inmodel``) has been 
adjusted using :doc:`ctlike` to the data, but if this is not the case you 
can request a maximum likelihood fit by setting the hidden parameter ``fit=yes``.

ctbutterfly writes the butterfly diagram into an ASCII file with 4 columns 
separated by a whitespace. Each row in the result file corresponds to a specific
energy. The meaning of the columns are:

* Energy in MeV
* Fitted intensity for that energy in ph/cm2/s/MeV
* Minimum intensity for that energy in ph/cm2/s/MeV
* Maximum intensity for that energy in ph/cm2/s/MeV

The butterfly diagram can be displayed using the ``show_butterfly.py`` script
in the example folder.


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.
 	 	 
``inmodel [file]``
    Input model XML file.
 	 	 
``srcname [string]``
    Name of source model for which the butterfly diagram should be computed.
 	 	 
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

``(fit = no) [boolean]``
    Performs maximum likelihood fitting of input model ignoring any provided
    covariance matrix.
 	 	 
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
    Specifies whether an existing output file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctbutterfly.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctlike`
:ref:`ctulimit`
:ref:`cterror`
