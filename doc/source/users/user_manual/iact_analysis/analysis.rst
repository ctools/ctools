.. _sec_analysis:

Analysing IACT data
===================

In order to analyse data from current IACTs, the data have to be retrieved. Information how
to copy IACT data to a local machine can be found here: :ref:`sec_copy`

Check available FITS production
-------------------------------
Before starting an analysis, it is important to know what kind of data need to be analysed.
In the current scheme of `storing IACT data <http://gamma-astro-data-formats.readthedocs.org/en/latest/data_storage/index.html>`_,
every FITS data set is represented with a unique name, e.g. ``"hess-fits-pa-release-1.0-Prod26-MppStd"``. This string
has to be passed to some set of ctools to find and use the data for analysis.

To check the available names of FITS data sets, it is recommended to run :ref:`csiactdata`.

.. code-block:: bash

  $ csiactdata
  Path were data is located [] /path/to/fits/data/
  
This will list names of available FITS data sets on the screen. In the 
following it is assumed the user has selected a data set called ``fits-data-name``.
The query for the path where the data is located can be
omitted by setting the environment variable ``VHEFITS`` to this locations. This might be convenient
since some other tools also query for the same parameter and can use ``VHEFITS`` instead.

Find observations
-----------------
To start an analysis, it is required to assemble a list of observations that should be used. For this purpose,
:ref:`csfindobs` searches for observations according to user criteria.

Search for sky pointings
^^^^^^^^^^^^^^^^^^^^^^^^
The most common search for observations is by pointing position in the sky.

.. code-block:: bash

  $ csfindobs
  Name of FITS production (Run csiactdata to view your options) [fits-data-name] 
  Right ascension [83.6331]
  Declination [22.01] 
  Search radius [2.5]
  Runlist outfile [runlist.lis]
  
  
Specifying additional search criteria
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The search for observations can be further constrained arbitrarily. For this purpose,
:ref:`csfindobs` contains a hidden parameter ``expression``. It is possible to use any
FITS selection expression that is `supported by cfitsio <http://www.isdc.unige.ch/integral/download/osa/doc/10.1/osa_um_intro/node38.html>`_.
Available properties (i.e. column names) for selection, can be found 
`here <http://gamma-astro-data-formats.readthedocs.org/en/latest/data_storage/obs_index/index.html>`_
or by simply browsing the observation index file with a FITS viewer.

Example:

.. code-block:: bash

  $ csfindobs expression="ZEN_PNT<30&&LIVETIME>1800"
  ...
  
This command selects all observations that have a zenith angle less than 
30 degrees and a livetime larger than 1800 seconds (i.e. 30 minutes).


Omitting the pointing selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You can also omit the pointing selection, for example in case where you 
don't know the sky coordinates of your object of interest. In that case, you 
may directly use the object name in the expression and pass ``NONE`` as
coordinate to :ref:`csfindobs`:

.. code-block:: bash

  $ csfindobs expression="OBJECT=='Crab Nebula'"
  Name of FITS production (Run csiactdata to view your options) [fits-data-name] 
  Right ascension [NONE]
  Runlist outfile [runlist.lis]

This will omit the selection of pointing positions.

.. note::

   By default, :ref:`csfindobs` only selects data of highest quality (i.e. QUALITY=0).
   You may overwrite this default by specifying the hidden parameter 
   ``min_qual``. For example, ``min_qual=1`` selects all data with a 
   looser quality criteria.

Create an observation list
--------------------------
The runlist ASCII file containing a list of selected observation IDs must now be converted to an
observation definition XML file. This file contains information about the location of the files that are required
for the analysis. The purpose of the tool :ref:`csiactobs` is to do this conversion.

.. code-block:: bash

  $ csiactobs
  Data storage name (Run csiactdata to view your options) [] fits-data-name
  Runlist file [runlist.lis]
  Number of free parameters per background model [1] 
  Output model XML file [bgmodels.xml] 
  Observation XML outfile [obs.xml]

This tools has now performed various tasks:

* Use the IACT data storage to convert the observation list ``runlist.lis`` into an observation definition XML file ``obs.xml``.
* Create a model XML file that contains the background models. In ctools, each observation has its own independent background model. The number of free parameters per model was set to 1, i.e. the normalisation of the background is left free for each observation. The models were saved in ``bgmodels.xml``.
* The tool uses an internal hierarchy for assigning the background models. The hidden parameter ``bkg_mod_hiera`` (default=irf|aeff|gauss) steers the order how background models should be created. In case the IRF background model is not available, the script automatically falls back to the background model from effective area (GCTAModelAeffBackground). 
* There are some further hidden parameters to steer the start parameters for the Aeff and Gaussian background model. Have a look at :ref:`csiactobs` to view a full list of parameters.
* In ``csiactobs.log`` (or on screen if ``debug=yes``), the script logs the complete energy range of the observation container. These values might be important for later usage (e.g. in binned analysis).

In case a sky model is already prepared, it is possible to also provide the hidden parameter ``inmodel``. The output
model XML file will then contain both, the background model and the input sky model:

.. code-block:: bash

  $ csiactobs inmodel="mymodel.xml"
  
Alternatively, models can be merged at any times using the simple tool :ref:`csmodelmerge`:

.. code-block:: bash

  $ csmodelmerge
  Input model XML files [bgmodels.xml crab.xml]
  Output model file [crab_models.xml]
  
Note that the number of files to merge is not limited to two. Detailled options how the input model XML file can
be passed is given on the reference page of :ref:`csmodelmerge`. It is also important to know that each model in the container
must have a unique name. This implies, for instance, merging the same XML model twice will result in an exception.

A list of available sky models can be found `here <http://gammalib.sourceforge.net/user_manual/modules/model.html>`_.
Of particular help to create sky models for your dataset is the section about :ref:`modelling CTA data <models>`.


Example XML files
-----------------
To get familiar with the XML syntax and format, example files for an observation container and a model container
are shown in the following.

Observation XML file
^^^^^^^^^^^^^^^^^^^^

.. code-block:: xml

	<?xml version="1.0" encoding="UTF-8" standalone="no"?>
	<observation_list title="observation list">
	  <observation name="Crab Nebula" id="11111" instrument="HESS">
	    <parameter name="EventList" file="/path/to/fits/file/events_11111.fits.gz" />
	    <parameter name="EffectiveArea" file="/path/to/fits/file/aeff_11111.fits.gz" />
	    <parameter name="PointSpreadFunction" file="/path/to/fits/file/psf_11111.fits.gz" />
	    <parameter name="EnergyDispersion" file="/path/to/fits/file/edisp_11111.fits.gz" />
	    <parameter name="Background" file="/path/to/fits/file/bgmodel_11111.fits.gz" />
	  </observation>
	  <observation name="Crab Nebula" id="11112" instrument="HESS">
	    <parameter name="EventList" file="/path/to/fits/file/events_11112.fits.gz" />
	    <parameter name="EffectiveArea" file="/path/to/fits/file/aeff_11112.fits.gz" />
	    <parameter name="PointSpreadFunction" file="/path/to/fits/file/psf_11112.fits.gz" />
	    <parameter name="EnergyDispersion" file="/path/to/fits/file/edisp_11112.fits.gz" />
	    <parameter name="Background" file="/path/to/fits/file/bgmodel_11112.fits.gz" />
	  </observation>
	  <observation name="Crab Nebula" id="11113" instrument="HESS">
	    <parameter name="EventList" file="/path/to/fits/file/events_11113.fits.gz" />
	    <parameter name="EffectiveArea" file="/path/to/fits/file/aeff_11113.fits.gz" />
	    <parameter name="PointSpreadFunction" file="/path/to/fits/file/psf_11113.fits.gz" />
	    <parameter name="EnergyDispersion" file="/path/to/fits/file/edisp_11113.fits.gz" />
	    <parameter name="Background" file="/path/to/fits/file/bgmodel_11113.fits.gz" />
	  </observation>
	  <observation name="Crab Nebula" id="11114" instrument="HESS">
	    <parameter name="EventList" file="/path/to/fits/file/events_11114.fits.gz" />
	    <parameter name="EffectiveArea" file="/path/to/fits/file/aeff_11114.fits.gz" />
	    <parameter name="PointSpreadFunction" file="/path/to/fits/file/psf_11114.fits.gz" />
	    <parameter name="EnergyDispersion" file="/path/to/fits/file/edisp_11114.fits.gz" />
	    <parameter name="Background" file="/path/to/fits/file/bgmodel_11114.fits.gz" />
	  </observation>
	</observation_list>

Model XML file
^^^^^^^^^^^^^^

.. code-block:: xml

	<?xml version="1.0" encoding="UTF-8" standalone="no"?>
	<source_library title="source library">
	  <source name="bkg_11111" type="CTAIrfBackground" instrument="HESS" id="11111">
	    <spectrum type="ConstantValue">
	      <parameter name="Value" value="1" error="0" scale="1" min="0.01" max="100" free="1" />
	    </spectrum>
	  </source>
	  <source name="bkg_11112" type="CTAIrfBackground" instrument="HESS" id="11112">
	    <spectrum type="ConstantValue">
	      <parameter name="Value" value="1" error="0" scale="1" min="0.01" max="100" free="1" />
	    </spectrum>
	  </source>
	  <source name="bkg_11113" type="CTAIrfBackground" instrument="HESS" id="11113">
	    <spectrum type="ConstantValue">
	      <parameter name="Value" value="1" error="0" scale="1" min="0.01" max="100" free="1" />
	    </spectrum>
	  </source>
	  <source name="bkg_11114" type="CTAIrfBackground" instrument="HESS" id="11114">
	    <spectrum type="ConstantValue">
	      <parameter name="Value" value="1" error="0" scale="1" min="0.01" max="100" free="1" />
	    </spectrum>
	  </source>
	    <source name="Crab" type="PointSource">
	    <spectrum type="PowerLaw">
	       <parameter name="Prefactor" scale="1e-17" value="3.0"  min="1e-07" max="1000.0" free="1"/>
	       <parameter name="Index"     scale="-1"    value="2.5" min="0.0"   max="+5.0"   free="1"/>
	       <parameter name="Scale"     scale="1e6"   value="1.0"  min="0.01"  max="1000.0" free="0"/>
	    </spectrum>
	    <spatialModel type="SkyDirFunction">
	      <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
	      <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
	    </spatialModel>
	  </source>
	</source_library>

.. note::

   It is important to ensure background models are properly linked to their respective observation.
   Therefore it is required to keep the attributes ``instrument`` and ``id`` the same for the observation
   and the corresponding background model. The tool :ref:`csiactobs` assures this automatically.

Run ctselect
------------
To prepare the data for analysis, cuts have to be applied to the event data. The selection is performed by :ref:`ctselect`.
This tool writes out selected event lists into the local directory. If the observation XML file contains several runs, it is recommended
to first create a separate folder and specify this folder in the hidden ``prefix`` argument.

.. code-block:: bash

  $ mkdir selected
  
.. code-block:: bash

  $ ctselect usethres=DEFAULT usepnt=yes prefix=selected/selected_
  Input event list or observation definition file [events.fits] obs.xml
  Lower energy limit (TeV) [0.1]
  Upper energy limit (TeV) [100.0]
  Radius of ROI (degrees) (0-180) [3.0] 2.5
  Start time (CTA MET in seconds) [0.0]
  End time (CTA MET in seconds) [0.0]
  Output event list or observation definition file [selected_events.fits] selected_obs.xml 
  
For IACT analysis, it is recommended to use the hidden parameter ``usethres="DEFAULT"``. This instructs :ref:`ctselect`
to extract the safe energy range from the instrument response functions and apply them to the data. This safe energy range
is thus superior to the energy limit passed via the user parameters. In addition, to analyse the complete field of view,
the parameter ``usepnt=yes`` uses, for each observation, the pointing position as centre for the selection radius.
The radius parameter is dependent on the intrument, for an instrument with a 5 degree field of view, a radius of 2.5 degrees
seems reasonable. The time selection is not applied in the above example; 
specifying 0 as start and end time skips the time selection.
For time-resolved analysis, it is important to know the MET time that is required to extract.
The result of the selection step is written into the observation XML file ``selected_obs.xml``, which now contains references
to the new selected event FITS files.

Unbinned analysis
-----------------
Once the data is selected, the easiest way to analyse is an unbinned analysis. Note that the input model XML
file must now contain the background models and source components to describe the field of view.

.. code-block:: bash

  $ ctlike
  Event list, counts cube or observation definition file [events.fits] selected_obs.xml
  Source model [$CTOOLS/share/models/crab.xml] crab_models.xml
  Source model output file [crab_results.xml]
  
The result of the fit was stored in ``crab_results.xml``. Note that fitted parameters, ``Prefactors`` in particular,
typically use MeV as energy unit. To monitor the progress of the fit on the screen, one can simply run with the option ``debug=yes``.
Alternatively, the logfile ``ctlike.log`` can be inspected after the fit. 

On default, energy dispersion is not considered in the fit. To switch on the usage of the energy migration matrix,
the hidden parameter ``edisp=yes`` can be provided. Note that this will cause a significant reduction of the computing
speed.

Stacked (binned) analysis
-------------------------
The stacked analysis mode is using a binned analysis where all observations are included and stacked into one event cube.
This analysis mode is much faster than unbinned analysis when having a large dataset (e.g. > 100 hours).
For this type of binned analysis, some intermediate data products have to be produced. The products are a binned data cube,
an exposure cube, a PSF cube, and a background cube. 

Bin events
^^^^^^^^^^

.. code-block:: bash

  $ ctbin
  Event list or observation definition file [events.fits] selected_obs.xml
  First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63]
  Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01]
  Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR]
  Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL]
  Image scale (in degrees/pixel) [0.02]
  Size of the X axis in pixels [200] 
  Size of the Y axis in pixels [200]
  Algorithm for defining energy bins (FILE|LIN|LOG) [LOG]
  Start value for first energy bin in TeV [0.1] 0.5
  Stop value for last energy bin in TeV [100.0] 50
  Number of energy bins [20]
  Output counts cube [cntcube.fits]

Note that bins only get filled if the bin of the cube is fully contained in the energy range and RoI of a considered observation.
It is therefore useful to provide the energy range given by :ref:`csiactobs` above. This ensures a maximum agreement between
observations and binning and reduces the loss of data.

Create exposure cube
^^^^^^^^^^^^^^^^^^^^
After binning the events into a three-dimensional cube, an exposure cube has to be computed.
The exposure is defined as the effective area times the dead-time corrected observation time.
Each observation from the input container gets stacked in the resulting cube. The exposure is stored
in units of :math:`cm^2 s`. The exposure cube does not have to contain the same binning as the event cube
but for simplicity, the event cube can be passed to adopt the binning parameters. Note, however, that the
exposure cube is defined in true sky coordinates and energy while the counts cube is defined in reconstructed
sky coordinates and energy. Consequently, the sky area and energy range covered by the exposure cube should be
slightly larger than that of the counts cube to accommodate for spill over of events due to the point spread function
and energy dispersion. By computing the exposure cube on the same grid as the counts cube, the spill over of events
from sources at the edge of cube will not be handled correctly. In the example below, however, no source at the edge of the
field of view is present. Therefore, for simplicity, the count cube is used as input cube to extract the binning.

This task of computing the exposure cube is is performed by :ref:`ctexpcube`. 

.. code-block:: bash

  $ ctexpcube
  Event list or observation definition file [NONE] selected_obs.xml
  Input counts cube file to extract exposure cube definition [NONE] cntcube.fits
  Output exposure cube file [expcube.fits]
  
Alternatively, the exposure cube can be created with different binning than the event cube:

.. code-block:: bash

  $ ctexpcube
  Input event list or observation definition XML file [NONE] selected_obs.xml 
  Input counts cube file to extract exposure cube definition [NONE] 
  First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 
  Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 
  Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR] 
  Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] 
  Image scale (in degrees/pixel) [0.02] 0.04
  Size of the X axis in pixels [200] 100
  Size of the Y axis in pixels [200] 100
  Lower energy limit (TeV) [0.5] 
  Upper energy limit (TeV) [50.0] 
  Number of energy bins [20] 30
  Output exposure cube file [expcube.fits] 

Create PSF cube
^^^^^^^^^^^^^^^
As a next step for the binned analysis, a cube containing the point spread function (PSF) must be computed.
Since the PSF cannot be stored by one single parameter, the PSF cube computed by :ref:`ctpsfcube` has a fourth dimension.
In each bin of the cube, the PSF is stored as a function of offset from source. The granularity of the
PSF histogram is determined by the hidden parameter ``anumbins`` (default: 200). Therefore, when passing
the event cube to adopt the sky binning for the PSF cube, the resulting FITS file can become quite large due to
the fourth dimension. Usually in IACT analysis, the PSF doesn't change too dramatically across the field of view.
Therefore the user can think about reducing the spatial binning of the PSF cube:

.. code-block:: bash

  $ ctpsfcube
  Input event list or observation definition XML file [NONE] selected_obs.xml 
  Input counts cube file to extract PSF cube definition [NONE] 
  First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 
  Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 
  Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR] 
  Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] 
  Image scale (in degrees/pixel) [1.0] 0.2
  Size of the X axis in pixels [10] 20
  Size of the Y axis in pixels [10] 20
  Lower energy limit (TeV) [0.1] 0.5
  Upper energy limit (TeV) [100.0] 50
  Number of energy bins [20] 
  Output PSF cube file [psfcube.fits]
  
Depending on the required PSF precision, one could reduce the number of offset bins via the hidden parameter
``anumbins``:

.. code-block:: bash

  $ ctpsfcube anumbins=100

Create background cube
^^^^^^^^^^^^^^^^^^^^^^
Last but not least, for a binned IACT analysis a cube containing the background rate in sky coordinates and
reconstructed energy has to be computed. This task is performed by :ref:`ctbkgcube`. The binning here can also differ
from the event cube. For simplicity, however, the example below uses the event cube to adopt the binning.

.. code-block:: bash

  $ ctbkgcube debug=yes
  Input event list or observation definition XML file [NONE] selected_obs.xml 
  Input counts cube file to extract background cube definition [NONE] cntcube.fits 
  Input model XML file [NONE] crab_models.xml 
  Output background cube file [bkgcube.fits] 
  Output model XML file [NONE] binned_models.xml

Note that this tool also requires the parameters of an input and output model. In the model XML file that came out of
:ref:`csiactobs`, one background model per observation is included. This models get merged and averaged in the background
sky cube. In the output model the background models per obseervation will be removed. Instead, a global background model
for the newly created background cube is included. Sky models present in the input model XML file will also be included in the 
new XML file, which subsequently can be used for binned :ref:`ctlike`.

Example for stacked model XML file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The output model of :ref:`ctbkgcube` looks the following:

.. code-block:: xml

	<?xml version="1.0" encoding="UTF-8" standalone="no"?>
	<source_library title="source library">
	  <source name="BackgroundModel" type="CTACubeBackground" instrument="CTA,HESS,MAGIC,VERITAS">
	    <spectrum type="PowerLaw">
	      <parameter name="Prefactor" value="1" error="0" scale="1" min="0.01" max="100" free="1" />
	      <parameter name="Index" value="0" error="0" scale="1" min="-5" max="5" free="1" />
	      <parameter name="Scale" value="1" scale="1e+06" free="0" />
	    </spectrum>
	  </source>
	  <source name="Crab" type="PointSource">
	    <spectrum type="PowerLaw">
	       <parameter name="Prefactor" scale="1e-17" value="3.0"  min="1e-07" max="1000.0" free="1"/>
	       <parameter name="Index"     scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
	       <parameter name="Scale"     scale="1e6"   value="1."  min="0.01"  max="1000.0" free="0"/>
	    </spectrum>
	    <spatialModel type="SkyDirFunction">
	      <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="0"/>
	      <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="0"/>
	    </spatialModel>
	  </source>
	</source_library>
	
The background model with ``type=CTACubeBackground`` is used to scale the background cube stored in the FITS file
created by :ref:`ctbkgcube`.

Run ctlike
^^^^^^^^^^
Having all the intermediate data products ready, a binned analysis can be conducted using :ref:`ctlike`.

.. code-block:: bash

  $ ctlike
  Input event list, counts cube or observation definition XML file [events.fits] cntcube.fits 
  Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits 
  Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits 
  Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits 
  Input model XML file [binned_models.xml]
  Output model XML file [binned_results.xml]
	
Note that when passing an event cube to :ref:`ctlike`, the tool behaves differently than in unbinned mode.
It queries directly for the additional ingredients for the binned analysis. It is important to pass the background model
generated by :ref:`ctbkgcube` here to ensure the proper modelling of the background in the fit.




