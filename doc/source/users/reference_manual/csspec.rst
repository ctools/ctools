.. _csspec:

csspec
======

Computes spectrum for a given source.


Synopsis
--------

This script computes the source spectrum by fitting a model in a given set
of spectral bins. The model fit per spectral bin is performed using :doc:`ctlike`
and the script provides the possibility to fix sources other than the
source of interest (hidden parameter ``fix_srcs``) or to fix the background
model component(s) (hidden parameter ``fix_bkg``). The script computes the
source flux and its uncertainty in each spectral bin, as well as the
significance of the source detection. Optionally, it also computes an upper
flux limit that is particularly useful in case that the source is not
significantly detected within a spectral bin (hidden parameter ``calc_ulim``).

The script works on binned and unbinned data.

For unbinned data, the spectral binning is either defined by a FITS file
containing the energy boundaries of each bin (option ``ebinalg=FILE``) or
as ``enumbins`` bins spread linearly  (option ``ebinalg=LIN``) or
logarithmically (option ``ebinalg=LOG``) from a minimum energy given by
``emin`` to a maximum energy given by ``emax``.

For binned data, all energy bins within the counts cube that overlap with
the energy range spanned by ``emin`` and ``emax`` are considered. The number
of spectral bins is only approximately determined by the ``enumbins`` parameter.
Naturally, the script cannot create more spectral bins than the number of
energy bins that are available in the cube. Also, in case that ``enumbins``
is smaller than the number of energy bins in the cube, the script will fit
several layers of the counts cube for each spectral bin. The number of 
layers fit is determined by the total number of energy bins divided by the
``enumbins`` parameter.

On output, the script will provide a FITS file with the fitted source 
spectrum in form of a binary table. Each row corresponds to a spectral bin.
The columns are the mean as well as the boundaries of the spectral bin, 
the fitted flux and flux error, the Test Statistics value (option
``calc_ts=yes``), the upper flux limit (option ``calc_ulim=yes``) and the
predicted number of events (only for unbinned data).


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
    Number of energy bins.
 	 	 
``ebinalg <FILE|LIN|LOG> [string]``
    Algorithm for defining energy bins.
 	 	 
``ebinfile [file]``
    Name of the file containing the energy bin definition.

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

``(publish = no) [boolean]``
    Specifies whether the spectrum should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing source spectrum output file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csspec.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctlike`
