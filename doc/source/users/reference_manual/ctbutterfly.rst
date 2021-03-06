.. _ctbutterfly:

ctbutterfly
===========

Computes butterfly diagram for a given spectral model.


Synopsis
--------

This tool calculates a butterfly diagram for a specific source according to 
its spectral model. The butterfly diagram is the envelope of all spectral
models that are within a given confidence limit compatible with the data.
The default method used for the calculation is Gaussian error propagation
using the covariance matrix from a maximum likelihood fit. By default a
confidence level of 68% is used, but the level can be adjusted using the
hidden ``confidence`` parameter. For power law models, an alternative
calculation method can be specified by setting the hidden parameter
``method=ENVELOPE``. By using this method, the envelope is computed by
evaluating for each energy the minimum and maximum intensity of all power
law models that fall within the error ellipse of the prefactor and index
parameters. The error ellipse is derived from the covariance matrix of a
maximum likelihood fit.

:ref:`ctbutterfly` assumes that the input model (parameter ``inmodel``) has been
adjusted using :doc:`ctlike` to the data, but if this is not the case you 
can request a maximum likelihood fit by setting the hidden parameter ``fit=yes``.

:ref:`ctbutterfly` writes the butterfly diagram into a FITS file with a binary table
extension.

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

``outfile [file]``
    Output butterfly FITS file.

``(fit = no) [boolean]``
    Performs maximum likelihood fitting of input model ignoring any provided
    covariance matrix.

``(method = "GAUSSIAN") <GAUSSIAN|ENVELOPE> [string]``
    Computation method.

``(confidence = 0.68) [real]``
    Confidence level for error computation.

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
    Maximum number of fit iterations.

``(matrix = "NONE") [file]``
    Input covariance matrix file (not used)

``(ebinalg = LOG) <FILE|LIN|LOG> [string]``
    Algorithm for defining energy bins. For ``FILE``, the energy bins are defined
    in a FITS file that is specified by the ``ebinfile`` parameter, for ``LIN``
    and ``LOG`` there will be ``enumbins`` energy bins spaced linearly or
    logarithmically between ``emin`` and ``emax``, respectively.

``emin [real]``
    Lower energy value for first energy bin (in TeV) if ``LIN`` or ``LOG``
    energy binning algorithms are used.

``emax [real]``
    Upper energy value for last energy bin (in TeV) if ``LIN`` or ``LOG``
    energy binning algorithms are used.

``(enumbins = 100) [integer]``
    Number of energy bins if ``LIN`` or ``LOG`` energy binning algorithms are
    used.

``ebinfile [file]``
    Name of the file containing the energy binning definition if ``ebinalg=FILE``.
    You may use :ref:`csebins` to generate a file with appropriate energy binning.


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
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = ctbutterfly.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctlike`
:ref:`ctulimit`
:ref:`cterror`
