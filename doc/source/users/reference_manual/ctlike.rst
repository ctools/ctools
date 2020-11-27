.. _ctlike:

ctlike
======

Perform maximum likelihood model fitting.


Synopsis
--------

This tool performs a maximum likelihood fitting of a model to unbinned, binned
or On/Off data. All model parameters marked with the attribute ``free="1"``
in the input model XML file will be adjusted to the value that maximises 
the probability that the observed data have been drawn from the model. The tool
can mix the analysis of unbinned, binned and On/Off data, and is also able to
perform a joint maximum likelihood analysis of data collected in separate
observations or with different instruments.

:ref:`ctlike` will automatically switch between unbinned, binned, On/Off and joint
analysis on basis of the file provided on input (parameter ``inobs``). Providing
an event list will lead to an unbinned analysis while providing a counts cube
will lead to a binned analysis. If several observations have been stacked into a
single counts cube, an exposure cube, a point spread function cube and a
background cube need also to be provided. :ref:`ctlike` will automatically prompt
for these parameters in case that a counts cube is given. Answering ``NONE``
will ignore the exposure cube, point spread function cube and background cube
and handle the data as a single observation using the specified instrument
response function (parameters ``caldb`` and ``irf``).

Providing an observation definition XML file containing an On/Off observation
will trigger an On/Off analysis. Multiple observations, including data collected
with different instruments, can also be handled by specifying an observation
definition XML file on input. Each of the observations will be kept separately
and associated with its appropriate instrument response functions, as opposed to
a stacked analysis where average response functions, computed using :doc:`ctexpcube`,
:doc:`ctpsfcube` and :doc:`ctbkgcube`, are used. If an observation definition XML
file is provided, :ref:`ctlike` will use the joint likelihood of all the
observations for parameter optimisation.

By default, :ref:`ctlike` will use the Poisson statistic for likelihood computation,
but for binned analysis also Gaussian statistic can be specified (parameter
``statistic``). In case of an On/Off analysis, data can be fitted using either the
``CSTAT`` or the ``WSTAT`` statistic. The former will adjust a background model to
the On and Off data, while the latter will assume that the background rates per
solid angle are identical in the On and Off regions, and will marginalise over
the corresponding background rates. The fit statistic can also be specified in
the observation definition XML file for each observation using the ``statistic``
attribute.

For all model components that have the ``tscalc="1"`` attribute set in the input
model XML file, :ref:`ctlike` will also compute the Test Statistic value that is a
measure of the source significance. Optionally, the spatial model parameters can
be kept fixed during that computation (parameter ``fix_spat_for_ts``).

Optionally, a refit can be requested after an initial fit (parameter ``refit``),
but normally this is not needed.

:ref:`ctlike` generates an output model XML file that contains the values of the
best fitting model parameters. For all free parameters, an ``error`` attribute
is added that provides the statistical uncertainty in the parameter estimate.
If computation of the Test Statistic has been requested for a model component,
a ``ts`` attribute providing the Test Statistic value is added. The output model
can be used as an input model for other ctools.

In addition, the covariance matrix can be saved into a file by specifying a
file name for the ``outcovmat`` parameter. If the file name ends with ``.fits``
or ``.fit`` a FITS file is written, otherwise a CSV file is written.


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [string]``
    Input model XML file.

``expcube [file]``
    Input exposure cube file.

``psfcube [file]``
    Input PSF cube file.

``bkgcube [file]``
    Input background cube file.

``edispcube [file]``
    Input energy dispersion cube file.

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
    Output FITS or CSV file to store covariance matrix.

``(statistic = DEFAULT) <DEFAULT|CSTAT|WSTAT|CHI2> [string]``
    Optimization statistic. ``DEFAULT`` uses the default statistic for all
    observations, which is ``CSTAT`` or the statistic specified in the
    observation definition XML file. ``CSTAT`` uses the C statistic for
    all observations, ``WSTAT`` uses the W statistic for all On/Off
    observations, and ``CHI2`` uses the Chi squared statistic for all
    binned or stacked observations.

``(refit = no) [boolean]``
    Perform refitting of solution after initial fit.

``(refit_if_failed = yes) [boolean]``
    Perform refitting of solution in case that the initial fit failed. Failures
    considered are a stalled fit, an exhaustion of the maximum number of fit
    iterations, or a significant difference between the number of observed and
    predicted events.

``(like_accuracy = 0.005) [real]``
    Absolute accuracy of maximum likelihood value. Reducing this value will
    increase the number of iterations and provide a more accurate maximum
    log likelihood value. Converserly, decreasing the value will result in less
    iterations at the expense of a less accurate maximum likelihood value.

``(max_iter = 50) [integer]``
    Maximum number of fit iterations.

``(fix_spat_for_ts = no) [boolean]``
    Fix spatial parameters for TS computation.


Standard parameters
-------------------

``(nthreads = 0) [integer]``
    Number of parallel processes (0=use all available CPUs).

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing output model file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = ctlike.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

None
