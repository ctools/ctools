.. _sec_high_level:

High level analysis tools for IACT data
=======================================

This section describes and explains how to create high level output. The
following scripts and tools are explained:

* :ref:`csobsinfo` (inspect and observation XML files)
* :ref:`csmodelinfo` (inspect model XML files)
* :ref:`ctulimit` (compute upper limits)
* :ref:`cterror` (compute asymmetric errors)
* :ref:`csspec` (compute spectral points)
* :ref:`cslightcrv` (compute a light curve)
* :ref:`csresmap` (compute a residual map)
* :ref:`cttsmap` (compute a test statistics map)

Inspecting observation definition files
---------------------------------------
After having created an observation definition file, the user might want to
know a summary of its content. Important quantities usually are the number of
observations, the total livetime, or the mean zenith angle. For this purpose,
:ref:`csobsinfo` inspects an observation XML container for its content.

.. code-block:: bash

  $ csobsinfo
  Event list, counts cube, or observation definition file [obs.xml] selected_obs.xml 
  
This dumps general information of the observation container in a logfile ``csobsinfo.log``.
For a direct output on the screen, pass the hidden parameter ``debug=yes`` on the command line:

.. code-block:: bash

  $ csobsinfo debug=yes
  Event list, counts cube, or observation definition file [obs.xml] selected_obs.xml 
  ...
  2016-02-19T15:43:08: +=========+
  2016-02-19T15:43:08: | Summary |
  2016-02-19T15:43:08: +=========+
  2016-02-19T15:43:08:  Unbinned observations .....: 4
  2016-02-19T15:43:08:  Number of events ..........: 6085
  2016-02-19T15:43:08:  Binned observations .......: 0
  2016-02-19T15:43:08:  Number of bins ............: 0
  2016-02-19T15:43:08: 
  2016-02-19T15:43:08: === Pointings ===
  2016-02-19T15:43:08:  Mean zenith angle .........: 46.89
  2016-02-19T15:43:08:  Mean azimuth angle ........: 13.63
  2016-02-19T15:43:08: 
  2016-02-19T15:43:08: === Energy range ===
  2016-02-19T15:43:08:  Emin ......................: 461.48 GeV
  2016-02-19T15:43:08:  Emax ......................: 37.8832 TeV
  2016-02-19T15:43:08: 
  2016-02-19T15:43:08: === Time range ===
  2016-02-19T15:43:08:  Start [MJD] ...............: 53343.9212153
  2016-02-19T15:43:08:  Stop [MJD] ................: 53347.9321065
  2016-02-19T15:43:08:  Total ontime ..............: 6877.00 s = 114.62 min = 1.91 h
  2016-02-19T15:43:08:  Total livetime ............: 6427.52 s = 107.13 min = 1.79 h
  ... 

Further options are described on the reference page :ref:`csobsinfo`.

Inspecting model XML files
--------------------------
Similar to the observation containers, model XML file can be inspected, too.
The tool :ref:`csmodelinfo` gives a summary of the model container. In particular,
the number of free parameters, number of sky and background models and also number of
parameters that are at their limits might be of interest.

.. code-block:: bash

  $ csmodelinfo debug=yes
  Input model definiton file XML file [$CTOOLS/share/models/crab.xml]
  Output DS9 region file [ds9.reg]
  2016-02-19T15:50:40: +=========+
  2016-02-19T15:50:40: | Summary |
  2016-02-19T15:50:40: +=========+
  2016-02-19T15:50:40: === Instrument specific models ===
  2016-02-19T15:50:40:  All .......................: 1
  2016-02-19T15:50:40:  CTA .......................: 1
  2016-02-19T15:50:40: === Model types ===
  2016-02-19T15:50:40:  PointSource ...............: 1
  2016-02-19T15:50:40:  CTAIrfBackground ..........: 1
  2016-02-19T15:50:40: 
  2016-02-19T15:50:40: +=======================+
  2016-02-19T15:50:40: | Parameter information |
  2016-02-19T15:50:40: +=======================+
  2016-02-19T15:50:40:  All parameters ............: 10
  2016-02-19T15:50:40:  Fixed parameters ..........: 6
  2016-02-19T15:50:40:  Free parameters (total) ...: 4
  2016-02-19T15:50:40:  Free background parameters : 2
  2016-02-19T15:50:40:  Free source parameters ....: 2
  2016-02-19T15:50:40:  Free spectral parameters ..: 2
  2016-02-19T15:50:40:  Parameters at limit .......: 0

Similar to :ref:`csobsinfo`, the script also generates a ds9 region file
including the positions of all the sky model components.

Compute upper limit
-------------------
Very often in gamma-ray astronomy sources are at the verge of detection or even not detectable.
In such cases, it is useful to derive an upper limit using :ref:`ctulimit`.

.. code-block:: bash

	$ ctulimit
	Input event list, counts cube or observation definition XML file [selected_obs.xml]
	Input model XML file [crab_models.xml]  
	Source of interest [Crab] 
	
The upper limit will be stored in the log file. To get the limit printed on screen, use the hidden parameter
``debug=yes``.

Compute asymmetric errors
-------------------------
When an analysis approaches its final state, it makes sense to have asymmetric errors on the parameters of the
source of interest. For this purpose, the tool :ref:`cterror` can be used:

.. code-block:: bash

	$ cterror debug=yes
	Input event list, counts cube or observation definition XML file [selected_obs.xml] 
	Input model XML file [crab_models.xml] 
	Source of interest [Crab] 
	Output model XML file [cterror_results.xml] 

The output model does not contain asymmetric errors yet. The positive and negative uncertainties can be read from the
logfile (or from screen if ``debug=yes`` was specified).

Compute spectral points
-----------------------
A very common task in astronomy is to compute spectral data points. To determine a spectral data point, a small energy
range is considered and the model prefactor and its uncertainty is evaluated. The tool :ref:`csspec` works on both, binned
and unbinned data. The hidden parameter ``edisp=yes`` can be specified in both cases to consider the energy migration
matrix in the fit.

.. note::
  
  The resulting spectral points are provided as a function of reconstructed energy. 

Unbinned
^^^^^^^^

.. code-block:: bash

	$ csspec debug=yes
	Input event list, counts cube, or observation definition XML file [events.fits] selected_obs.xml 
	Input model XML file [$CTOOLS/share/models/crab.xml] crab_models.xml 
	Source name [Crab] 
	Algorithm for defining energy bins (FILE|LIN|LOG) [LOG] 
	Lower energy limit for spectral points (TeV) [0.1] 0.5
	Upper energy limit for spectral points (TeV) [100.0] 50.0
	Number of spectral points (1-10000) [20] 10
	Output spectrum file [spectrum.fits] 
	
Binned
^^^^^^

.. code-block:: bash

	$ csspec debug=yes
	Input event list, counts cube, or observation definition XML file [cntcube.fits] 
	Input exposure cube file (only needed for stacked analysis) [expcube.fits] 
	Input PSF cube file (only needed for stacked analysis) [psfcube.fits] 
	Input background cube file (only needed for stacked analysis) [bkgcube.fits] 
	Input model XML file [binned_models.xml] 
	Source name [Crab] 
	Number of spectral points (1-10000) [10] 
	Lower energy limit for spectral points (TeV) [0.5] 
	Upper energy limit for spectral points (TeV) [50.0] 
	Output spectrum file [spectrum.fits]  

Plot spectral points
^^^^^^^^^^^^^^^^^^^^
Instead of supporting a plotting library, simple example scripts are available to visualise data products.
Have a look at ``$CTOOLS/examples/show_spectrum.py`` how a spectrum can be plotted using ``matplotlib``.
To have first glance at the above computed spectrum one can use this script in the following way:

.. code-block:: bash
  
  $ python $CTOOLS/examples/show_spectrum.py spectrum.fits


Compute light curves
--------------------
Euqally common in astronomy are light curves from a time-series analysis. The tool :ref:`cslightcrv` takes as input only
an unbinned observation container (or a single event list). Counts cube are not possible as inout since the time informstion
is lost during the binning procedure. Nevertheless, in case of large time bins, it is possible to require a binned analysis
in each time bin. The tool will accordingly slice the event list and create the data products for the stacked analysis in each
time span.

.. code-block:: bash

	$ cslightcrv
	Input event list, counts cube, or observation definition XML file [selected_obs.xml]
	Input model XML file [crab_models.xml]  
	Source name [Crab] 
	Algorithm for defining time bins (FILE|LIN|GTI) [GTI] 
	Number of energy bins per light curve bin (0=unbinned) [0] 
	Lower energy limit of events (TeV) [0.1] 0.5
	Upper energy limit of events (TeV) [100.0] 50
	Output light curve file [lightcurve.fits]  

Compute a residual map
----------------------
A frequent means to visually inspect the fitted model with respect to the input data is to create residual maps.
Using :ref:`csresmap` such a map can easily be computed. The tool internally bins the data according to user parameters.
(if data is not already provided in a binned state). Taking into account the corresponding IRFs,
:reF:`csresmap` computes a model map (running :ref:`ctmodel`). Subsequently, there are several choices how data
and model should be compared. There are three options:

* ``SUB``: the subtraction of the model from the counts. The resulting map will display differences in absolute counts.
* ``SUBDIV``: the subtraction and division by the model. The resulting map will display relative differences of the data with respect to the model.
* ``SUBDIVSQRT``: the subtraction and division by the square root of the model. The resulting map will display an approximation of a residual significance map. In case of sufficient count statistic per bin, the bin value represents the significance. 

Example unbinned observation container
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

	$ csresmap
	Input event list, counts cube, or observation definition XML file [selected_obs.xml]
	Input model XML file [crab_results.xml]
	First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 
	Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 
	Coordinate System (CEL|GAL) [CEL] 
	Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR] 
	Size of the X axis in pixels [100]
	Size of the Y axis in pixels [100]
	Pixel size (deg/pixel) [0.02]
	Residual map computation algorithm (SUB|SUBDIV|SUBDIVSQRT) [SUBDIV]
	Output residual map file [resmap.fits] 

Example binned/stacked observation container
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The tool behaves differently if the input observation parameter is a count cube. It will query for the additional binned reponse components.

.. code-block:: bash

	$ csresmap 
	Input event list, counts cube, or observation definition XML file [cntcube.fits]
	Input model cube file (generated with ctmodel) [NONE] 
	Input exposure cube file (only needed for stacked analysis) [expcube.fits]
	Input PSF cube file (only needed for stacked analysis) [psfcube.fits]  
	Input background cube file (only needed for stacked analysis) [bkgcube.fits] 
	Input model XML file [binned_results.xml] 
	Residual map computation algorithm (SUB|SUBDIV|SUBDIVSQRT) [SUBDIV] 
	Output residual map file [resmap.fits] 
	
Note that the tool also queries for an input model cube file (which we set to ``NONE`` here). This is very convenient in case the model cube has already been precomputed
using :ref:`ctmodel`.

Example with input model cube
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In case a model cube can be provided, :ref:`csresmap` simply collapses both cube into one skymap and applies the specified algorithm.

.. code-block:: bash

	$ csresmap
	Input event list, counts cube, or observation definition XML file [cntcube.fits]
	Input model cube file (generated with ctmodel) [modcube.fits] 
	Residual map computation algorithm (SUB|SUBDIV|SUBDIVSQRT) [SUB] 
	Output residual map file [resmap.fits] 


Compute a test statistics (TS) map
----------------------------------
The test statistics quantity is a more precise measurement of statistic than the ``SUBDIVSQRT`` option in :ref:`csresmap`. 
To get familiar with the concept of this procedure, :ref:`read more about TS maps <sec_tsmap>`.

.. code-block:: bash

	$ cttsmap
	Input event list, counts cube or observation definition XML file [selected_obs.xml]  
	Input model XML file [crab_results.xml] 
	Test source name [Crab] 
	First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 
	Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 
	Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR] 
	Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] 
	Image scale (in degrees/pixel) [0.05]
	Size of the X axis in pixels [20]
	Size of the Y axis in pixels [20]
	Output Test Statistic map file [tsmap.fits] 

The output file will contain a sky map holding the TS map. In complementary extensions, for each free parameter of the source model, a map displaying the parameter
value is added to the output file. For instance a map of the prefactor can be visualised with `DS9 <http://ds9.si.edu/>`_ as follows:

.. code-block:: bash
  
  $ ds9 tsmap.fits[Prefactor] 

Since :ref:`cttsmap` has to run a parameter optimisation in each skymap bin, it is very time comsuming to
compute a fine-granulated TS sky map. Read here how to split up the :ref:`computation into several jobs <sec_tips>`. 
Similar as many other tools :ref:`cttsmap` can also work on binned observation input:

.. code-block:: bash
  
	$ cttsmap debug=yes
	Input event list, counts cube or observation definition XML file [cntcube.fits]
	Input exposure cube file (only needed for stacked analysis) [expcube.fits] 
	Input PSF cube file (only needed for stacked analysis) [psfcube.fits] 
	Input background cube file (only needed for stacked analysis) [bkgcube.fits] 
	Input model XML file [binned_results.xml] 
	Test source name [Crab] 
	First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 
	Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 
	Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR] 
	Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] 
	Image scale (in degrees/pixel) [0.05] 
	Size of the X axis in pixels [20] 
	Size of the Y axis in pixels [20] 
	Output Test Statistic map file [tsmap.fits] 

	
