.. _csspec:

csspec
======

Computes spectrum for a given source.


Synopsis
--------

This script computes the source spectrum by fitting a model in a given set
of energy bins. The model fit per energy bin is performed using :doc:`ctlike`
and the script provides the possibility to fix sources other than the
source of interest (hidden parameter ``fix_srcs``) or to fix the background
model component(s) (hidden parameter ``fix_bkg``). The script computes the
source flux and its uncertainty in each energy bin, as well as the significance
of the source detection. Optionally, it also computes an upper flux limit
that is particularily useful in case that the source is not significantly
detected within an energy bin (hidden parameter ``calc_ulim``). The script 
works on both, binned and unbinned observation containers. In case of binned input,
the script will determine from the ``enumbins`` parameter, which energy layers 
of the input cube will be merged for one spectral bin. Of course, the script
cannot create more spectral bins as available energy bins in the cube. 
Also note that for the moment, the Npred column in the output file is not 
filled for binned analyses.

On output, the script will provide a FITS file with the fitted source 
spectrum.


General parameters
------------------

``(inobs = NONE) [file]``
    Input event list, counts cube or observation definition XML file.

``inmodel [file]``
    Input model XML file.

``srcname [string]``
    Name of the source in the source model XML file which should be used
    for sensitivity computation.

``outfile [file]``
    Name of the source spectrum output file.

``expcube [file]``
    Exposure cube file (only needed for stacked analysis).

``psfcube [file]``
    PSF cube file (only needed for stacked analysis).

``bkgcube [file]``
    Background cube file (only needed for stacked analysis).

``caldb [string]``
    Calibration database.
 	 	 
``irf [string]``
    Instrumental response function.

``(edisp = no) [boolean]``
    Apply energy dispersion to response computation?

``emin [real]``
    Lower energy limit of events (in TeV).
 	 	 
``emax [real]``
    Upper energy limit of events (in TeV).
 	 	 
``enumbins [integer]``
    Number of energy bins (0=unbinned).
 	 	 
``(ebinalg = LOG) <FILE|LIN|LOG> [string]``
    Algorithm for defining energy bins.
 	 	 
``(calc_ts = yes) [boolean]``
    Compute TS for each spectral point?

``(calc_ulim = yes) [boolean]``
    Compute upper limit for each spectral point?

``(fix_srcs = yes) [boolean]``
    Fix other sky model parameters?

``(fix_bkg = no) [boolean]``
    Fix background model parameters?


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

``(logfile = csspec.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctlike`
