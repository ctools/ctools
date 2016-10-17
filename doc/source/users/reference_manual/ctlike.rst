.. _ctlike:

ctlike
======

Perform maximum likelihood model fitting.


Synopsis
--------

This tool performs a maximum likelihood fitting of a model to unbinned or 
binned data. All model parameters marked with the attribute ``free="1"`` 
in the input model XML file will be adjusted to the value that maximises 
the probability that the observed data have been drawn from the model.
The tool can mix the analysis of unbinned and binned data, and is also able
to perform a joint maximum likelihood analysis of data collected in 
separate observations or with different instruments.

ctlike will automatically switch between unbinned, binned and joint analysis
on basis of the file provided on input (parameter ``inobs``). Providing an 
event list will lead to an unbinned analysis while providing a counts cube 
will lead to a binned analysis. If several observations have been stacked
into a single counts cube, an exposure cube, a point spread function cube and
a background cube need also to be provided. ctlike will automatically prompt
for these parameters in case that a counts cube is given. Answering ``NONE``
will ignore the exposure cube, point spread function cube and background cube
and handle the data as a single observation using the specified instrument
response function (parameters ``caldb`` and ``irf``).

Multiple observations, including data collected with different instruments,
can be handled by specifying an observation definition file on input. Each of
the observations will be kept separately and  associated with its appropriate
instrument response functions, as opposed to a stacked analysis where average
response functions, computed using :doc:`ctexpcube`, :doc:`ctpsfcube` and :doc:`ctbkgcube`,
are used. If an observation definition file is provided, ctlike will use the
joint likelihood of all the observations for parameter optimisation.

By default, ctlike will use the Poisson statistics for likelihood computation,
but for binned analysis also Gaussian statistics can be specified (parameter
``stat``). Optionally, a refit can be requested after an initial fit (parameter
``refit``), but normally this is not needed. For all model components that
have the ``tscalc="1"`` attribute set, ctlike will also compute the Test
Statistics value that is a measure of the source significance. Optionally,
the spatial model parameters can be kept fixed during that computation
(parameter ``fix_spat_for_ts``).

ctlike generates an output model XML file that contains the values of the 
best fitting model parameters. For all free parameters, an ``error`` attribute
is added that provides the statistical uncertainty in the parameter estimate.
In addition, the entire output of the covariance matrix to a separate FITS file
can be queried (parameter ``outcovmat``). The parameter order given in the file
corresponds to the arrangement of the covariance matrix entries. If computation
of the Test Statistics has been requested for a model component, a ``ts``
attribute providing the Test Statistics value is added. The output model can be
used as an input model for other ctools.



General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [string]``
    Input model XML file.
 	 	 
``expcube [file]``
    Input exposure cube file (only needed for stacked analysis).

``psfcube [file]``
    Input PSF cube file (only needed for stacked analysis).

``bkgcube [file]``
    Input background cube file (only needed for stacked analysis).

``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrument response function.
 	 	 
``(edisp = no) [boolean]``
    Applies energy dispersion to response computation.

``outmodel [string]``
    Output model XML file with values and uncertainties updated by
    the maximum likelihood fit.

``(outcovmat = NONE) [string]``
    Output FITS file to store covariance matrix.

``(stat = POISSON) [string]``
    Fitting statistics (POISSON or GAUSSIAN; only affects binned analysis).
 	 	 
``(refit = no) [boolean]``
    Perform refitting of solution after initial fit.

``(fix_spat_for_ts = no) [boolean]``
    Fix spatial parameters for TS computation.


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
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctlike.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

None
