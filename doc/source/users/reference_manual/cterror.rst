.. _cterror:

cterror
=======

Computes parameter errors for a source model from the likelihood profiles.


Synopsis
--------

This tool computes the parameter errors for a specific source model using
the likelihood profiles. Starting from the maximum likelihood model parameters,
the tool finds the minimum and maximum model parameters that lead to a decrease
of the likelihood that corresponds to a given confidence level. By default a
confidence level of 68% is used, but this level can be adjusted using the hidden
``confidence`` parameter.

:ref:`cterror` generates an output model XML file that contains the values of the
best fitting model parameters. For all free parameters, an ``error`` attribute
is added that provides the statistical uncertainty in the parameter estimate
as obtained from the likelihood profile. While cterror computes asymmetrical
errors, which are written into the log file, the XML file will contain the 
mean error that is obtained by computing the mean of the negative and the
positive parameter errors.


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [file]``
    Input model XML file.

``srcname [string]``
    Name of source model for which the parameter errors should be computed.

``expcube [file]``
    Input exposure cube file.

``psfcube [file]``
    Input PSF cube file.

``edispcube [file]``
    Input energy dispersion cube file.

``bkgcube [file]``
    Input background cube file.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrumental response function.

``(edisp = no) [boolean]``
    Applies energy dispersion to response computation.

``outmodel [file]``
    Output model XML file with updated error information.

``(confidence = 0.68) [real]``
    Confidence level for error computation.

``(tol = 1e-3) [real]``
    Computation tolerance.

``(statistic = DEFAULT) <DEFAULT|CSTAT|WSTAT|CHI2> [string]``
    Optimization statistic. ``DEFAULT`` uses the default statistic for all
    observations, which is ``CSTAT`` or the statistic specified in the
    observation definition XML file. ``CSTAT`` uses the C statistic for
    all observations, ``WSTAT`` uses the W statistic for all On/Off
    observations, and ``CHI2`` uses the Chi squared statistic for all
    binned or stacked observations.

``(like_accuracy = 0.005) [real]``
    Absolute accuracy of maximum likelihood value. Reducing this value will
    increase the number of iterations and provide a more accurate maximum
    log likelihood value. Converserly, decreasing the value will result in less
    iterations at the expense of a less accurate maximum likelihood value.

``(max_iter = 50) [integer]``
    Maximum number of iterations before stopping the likelihood profile
    computations.


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
    Specifies whether an existing output file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = cterror.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:ref:`ctlike`
:ref:`ctulimit`
:ref:`ctbutterfly`
