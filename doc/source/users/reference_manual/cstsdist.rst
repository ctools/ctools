.. _cstsdist:

cstsdist
========

Generates the TS distribution for a given source.


Synopsis
--------

This script generates the Test Statistics (TS) distribution for a given 
source by repeatedly computing the TS value for ``ntrials`` simulated data 
sets. The Test Statistics is defined as twice the log-likelihood difference
that is obtained when fitting simulated data with and without a given
source model component.

:ref:`cstsdist` supports multiprocessing. By default each simulation/analysis will
be performed in parallel over as many processes as the number of CPUs available on your
machine. The maximum number of parallel processes can be set by the user through the
``nthreads`` hidden parameter.

:ref:`cstsdist` will create on output a FITS file containing a table with one row
per TS computation. The columns in the table comprise the TS value, the number
of observed and fitted events, and the values and errors for all fitted
parameters.

From the output file, TS distribution plots can be generated using for
example the ``show_ts_distribution.py`` script in the examples folder. The
script requires matplotlib for plotting.

.. warning::
   This script does not work for On/Off observations. If an observation
   definition XML file is specified the script assumes that all observations
   are event lists.


General parameters
------------------

``inobs [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [file]``
    Input model XML file.

``srcname [string]``
    Name of the source in the source model XML file which should be used
    for Test Statistics computation.

``expcube [file]``
    Input exposure cube file.

``psfcube [file]``
    Input PSF cube file.

``edispcube [file]``
    Input energy dispersion cube file.

``bkgcube [file]``
    Input background cube file.

``(instrument = CTA) [string]``
    Name of Cherenkov telescope.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrumental response function.

``(edisp = no) [boolean]``
    Apply energy dispersion to response computation.

``(deadc = 0.98) [real]``
    Average deadtime correction factor.

``outfile [file]``
    Output TS distribution FITS file.

``ntrials [integer]``
    Number of Monte Carlo samples.

``(statistic = DEFAULT) <DEFAULT|CSTAT|WSTAT|CHI2> [string]``
    Optimization statistic. ``DEFAULT`` uses the default statistic for all
    observations, which is ``CSTAT`` or the statistic specified in the
    observation definition XML file. ``CSTAT`` uses the C statistic for
    all observations, ``WSTAT`` uses the W statistic for all On/Off
    observations, and ``CHI2`` uses the Chi squared statistic for all
    binned or stacked observations.


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
    Enables debug mode. In debug mode the executable will dump any log file
    output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = cstsdist.log) [string]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctlike`
:doc:`cspull`
