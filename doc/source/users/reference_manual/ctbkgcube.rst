.. _ctbkgcube:

ctbkgcube
=========

Generate background cube for a stacked analysis.


Synopsis
--------

This tool generates a background cube for a stacked analysis based on an
input model. A background cube is a 3-dimensional cube spanned by Right
Ascension or Galactic longitude, Declination or Galactic latitude, and energy.
The energy binning may be either linear, logarithmic, or custom defined using
an input file. The input model is used to predict the expected number of
background counts in each background cube bin.

:ref:`ctbkgcube` requires on input the event list or observation definition file
that has been used in the generation of the counts cube using :doc:`ctbin`.
To assure consistency between an existing counts cube and the
corresponding background cube, the counts cube needs also to be specified via
the ``incube`` parameter.

:ref:`ctbkgcube` generates a background cube FITS file comprising two extensions.
The primary extension contains a 3-dimensional image that contains the 
background cube values. The next extension named ``EBOUNDS`` contains a
binary table that defines the energy boundaries of the background cube.

:ref:`ctbkgcube` generates also an output model XML file that can serve as input 
for a maximum likelihood analysis. The output model XML file is a copy of
the input model XML file where the input background model has been replaced
by a background model of type ``CTACubeBackground``. The ``CTACubeBackground``
background model instructs any tool analysing binned data to extract 
background information from a background cube. The ``CTACubeBackground``
model has a spectral component that can be adjusted in a maximum 
likelihood fit to accomodate for uncertainties in the prediction of the 
energy dependence of the background rate. ctbkgcube will use a power law
as spectral component, but you can replace this by any component of your
choice.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``incube [file]``
    Input counts cube file to extract background cube definition.

``inmodel [file]``
    Input model XML file.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrument response function.

``outcube [file]``
    Output background cube file.

``outmodel [file]``
    Output model XML file.


Standard parameters
-------------------

``(publish = no) [boolean]``
    Specifies whether the background cube should be published on VO Hub.

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing output background cube file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = ctbkgcube.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctbin`
