.. _sec_analysis:

Analysing IACT data
===================

In order to analyse data from current IACTs, you first have to retrieve the
data. Information on how to do that can be found in the section :ref:`sec_copy`.

Check available FITS production
-------------------------------
Before you start an analysis, it is important that you know what kind of data
need to be analysed. In the current scheme of
`storing IACT data <http://gamma-astro-data-formats.readthedocs.org/en/latest/data_storage/index.html>`_,
every FITS data set is represented by a unique name, e.g.
``"hess-fits-pa-release-1.0-Prod26-MppStd"``.
You have to pass this string to some set of scripts to find and use the data
for analysis.

To check the available names of FITS data sets, you should run the
:ref:`csiactdata` script:

.. code-block:: bash

  $ csiactdata
  Path were data is located [] /path/to/fits/data/
  
This script will list the names of available FITS data sets into the console.
In the following it is assumed that you have selected a data set named
``fits-data-name``.
You can omit the query for the path where the data is located by setting the
environment variable ``VHEFITS`` to this location. Also other scripts will
benefit from having this environment variable set, and this avoids that you
have to type in the path name every time you use one of these scripts.


Find observations
-----------------

Now that you know which data are available you have to assemble a list of
observations (a.k.a. runs) that you want to use in an analysis. You will
use the :ref:`csfindobs` for that task, which allows you to select the
observations according to a number of user criteria.


Select observations according to pointing direction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most common search for observations is by pointing direction on the sky,
and you do this by simply invoking the :ref:`csfindobs` script:

.. code-block:: bash

  $ csfindobs
  Name of FITS production (Run csiactdata to view your options) [fits-data-name] 
  Right ascension of selection region centre (deg) [83.63]
  Declination of selection region centre (deg) [22.01]
  Search radius of selection region (deg) [2.5]
  Output runlist file [runlist.lis]

The script will create an ASCII file named ``runlist.lis`` that contains the
identifiers for all observations with pointing directions located within 2.5°
around the position of the Crab nebula.


Select observations using additional criteria
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You may add further selection criteria using the hidden ``expression`` parameter.
You may use any FITS selection expression that is
`supported by the cfitsio row filtering specification <https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node97.html>`_.
Column names that can be used in the expression can be found
`here <http://gamma-astro-data-formats.readthedocs.org/en/latest/data_storage/obs_index/index.html>`_
or by simply browsing the observation index file with a FITS viewer. For example

.. code-block:: bash

  $ csfindobs expression="ZEN_PNT<30&&LIVETIME>1800"

will select all observations that have a zenith angle less than 30° and a
livetime larger than 1800 seconds.


Select observations independent of pointing direction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also omit the pointing selection, for example in case that you do not
know the sky coordinates of your object of interest. In that case, you may
directly use the object name in the expression and pass ``NONE`` as coordinate
to :ref:`csfindobs`:

.. code-block:: bash

  $ csfindobs expression="OBJECT=='Crab Nebula'"
  Name of FITS production (Run csiactdata to view your options) [fits-data-name] 
  Right ascension of selection region centre [NONE]
  Output runlist file [runlist.lis]

.. note::

   By default, :ref:`csfindobs` only selects data of highest quality (i.e. QUALITY=0).
   You may overwrite this default by specifying the hidden parameter 
   ``min_qual``. For example, ``min_qual=1`` selects all data with a 
   looser quality criteria.


Create an observation list
--------------------------

As next step, you must convert the runlist ASCII file into an observation
definition XML file. The observation definition XML file will contain the
file names of all files that are needed for the analysis. You do this
conversion with the :ref:`csiactobs` script:

.. code-block:: bash

  $ csiactobs
  Data storage name [] fits-data-name
  Input runlist file [runlist.lis]
  Number of free parameters per background model [1] 
  Output model definition XML file [bgmodels.xml]
  Output observation definition XML file [obs.xml]

The :ref:`csiactobs` script will create two output files: the observations
definition XML file ``obs.xml`` and an output model definition XML file
``bgmodels.xml``. To generate ``obs.xml``, :ref:`csiactobs` has used the
IACT data storage and extracted the relevant file names. ``bgmodels.xml`` is
a file that is used for background modeling, where each observation will have
its own independent background model. In the example above, you have set the
number of free parameters per background model to one, hence the normalisation
of the background model for each observation will be a free parameter that is
later adjusted by the maximum likelihood fit.

There are some further hidden parameters to steer the start parameters for the
effective area and the Gaussian background models. You may have a look at
:ref:`csiactobs` to see the full list of available parameters.

You may also note that the :ref:`csiactobs` script has created a ``csiactobs.log``
file that logs the complete energy range of the observations that have been
selected. These values may be important for later, in particular if you
intend to do a stacked analysis.

In case that you have already a model definition XML file that describes a model
of the celestial source distribution (a.k.a. a sky model), you may provide this
sky model to :ref:`csiactobs` using the hidden ``inmodel`` parameter. The
output model definition XML file will then contain the sky model and the
background model, and you can use the model directly in a maximum likelihood
fit.

.. code-block:: bash

  $ csiactobs inmodel="mymodel.xml"
  
Alternatively, you can merged models using the :ref:`csmodelmerge` script:

.. code-block:: bash

  $ csmodelmerge
  Input model definition XML files [mymodel.xml bgmodels.xml]
  Output model definition XML file [combined_model.xml]

.. note::

   The number of files that can be merged with :ref:`csmodelmerge` is not
   limited to two. The list of input file names may be either separated by
   whitespace or semi-colons, can be specified using wildcards, or can be given
   in an ASCII file (see :ref:`csmodelmerge`).

.. warning::

   Each model component in the input model definition XML files that are provided
   to :ref:`csmodelmerge` must have a different and unique name. Merging for
   example the same XML model twice will raise an exception.

A list of available sky models can be found
`here <http://gammalib.sourceforge.net/user_manual/modules/model.html>`_.
If you are not familiar with creating sky models you should read the section
about :ref:`modelling CTA data <models>`.


Example XML files
-----------------
To get familiar with the XML syntax and format, example files for an observation
definition XML file and a model definition XML file are shown in the following.

Observation definition XML file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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


Model definition XML file
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: xml

	<?xml version="1.0" encoding="UTF-8" standalone="no"?>
	<source_library title="source library">
	  <source name="bkg_11111" type="CTAIrfBackground" instrument="HESS" id="11111">
	    <spectrum type="ConstantValue">
	      <parameter name="Normalization" value="1" error="0" scale="1" min="0.01" max="100" free="1" />
	    </spectrum>
	  </source>
	  <source name="bkg_11112" type="CTAIrfBackground" instrument="HESS" id="11112">
	    <spectrum type="ConstantValue">
	      <parameter name="Normalization" value="1" error="0" scale="1" min="0.01" max="100" free="1" />
	    </spectrum>
	  </source>
	  <source name="bkg_11113" type="CTAIrfBackground" instrument="HESS" id="11113">
	    <spectrum type="ConstantValue">
	      <parameter name="Normalization" value="1" error="0" scale="1" min="0.01" max="100" free="1" />
	    </spectrum>
	  </source>
	  <source name="bkg_11114" type="CTAIrfBackground" instrument="HESS" id="11114">
	    <spectrum type="ConstantValue">
	      <parameter name="Normalization" value="1" error="0" scale="1" min="0.01" max="100" free="1" />
	    </spectrum>
	  </source>
	    <source name="Crab" type="PointSource">
	    <spectrum type="PowerLaw">
	       <parameter name="Prefactor"   scale="1e-17" value="3.0"  min="1e-07" max="1000.0" free="1"/>
	       <parameter name="Index"       scale="-1"    value="2.5" min="0.0"   max="+5.0"   free="1"/>
	       <parameter name="PivotEnergy" scale="1e6"   value="1.0"  min="0.01"  max="1000.0" free="0"/>
	    </spectrum>
	    <spatialModel type="PointSource">
	      <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
	      <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
	    </spatialModel>
	  </source>
	</source_library>

.. note::

   It is important to ensure background models are properly linked to their
   respective observation. Therefore it is required to keep the attributes
   ``instrument`` and ``id`` the same for the observation and the corresponding
   background model. The tool :ref:`csiactobs` assures this automatically.


Run ctselect
------------
To prepare the data for analysis, cuts have to be applied to the event data.
The selection is performed by :ref:`ctselect`.
This tool writes out selected event lists into the local directory. If the
observation definition XML file contains several runs, it is recommended
to first create a separate folder and specify this folder in the hidden
``prefix`` argument.

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
  
For IACT analysis, it is recommended to use the hidden parameter
``usethres="DEFAULT"``. This instructs :ref:`ctselect` to extract the safe
energy range from the instrument response functions and apply them to the
data. This safe energy range is thus superior to the energy limit passed via
the user parameters. In addition, to analyse the complete field of view, the
parameter ``usepnt=yes`` uses, for each observation, the pointing position as
centre for the selection radius. The radius parameter is dependent on the
intrument, for an instrument with a 5° field of view, a radius of 2.5°
seems reasonable. The time selection is not applied in the above example; 
specifying 0 as start and end time skips the time selection. For time-resolved
analysis, it is important to know the MET time that is required to extract.
The result of the selection step is written into the observation XML file
``selected_obs.xml``, which now contains references to the new selected event
FITS files.


Unbinned analysis
-----------------
Once the data is selected, the easiest way to analyse is an unbinned analysis.
Note that the input model definition XML file must now contain the background
models and source components to describe the field of view.

.. code-block:: bash

  $ ctlike
  Event list, counts cube or observation definition file [events.fits] selected_obs.xml
  Source model [$CTOOLS/share/models/crab.xml] crab_models.xml
  Source model output file [crab_results.xml]
  
The result of the fit was stored in ``crab_results.xml``. Note that fitted
parameters, ``Prefactors`` in particular, typically use MeV as energy unit.
To monitor the progress of the fit on the screen, one can simply run with the
option ``debug=yes``. Alternatively, the logfile ``ctlike.log`` can be inspected
after the fit.

On default, energy dispersion is not considered in the fit. To switch on the
usage of the energy migration matrix, the hidden parameter ``edisp=yes`` can
be provided. Note that this will cause a significant reduction of the
computing speed.


Stacked analysis
----------------
In a stacked analysis the events of all observations are stacked into a single
counts cube. This analysis mode is much faster than unbinned analysis when
having a large dataset (e.g. > 100 hours). For this type of analysis, some
intermediate data products have to be produced. The products are a binned counts
cube, an exposure cube, a point spead function cube, optionally an energy
dispersion cube, and a background cube.

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


Create exposure cube
^^^^^^^^^^^^^^^^^^^^
After binning the events into a three-dimensional cube, an exposure cube has to
be computed. The exposure is defined as the effective area times the dead-time
corrected observation time. Each observation from the input container gets
stacked in the resulting cube. The exposure is stored in units of :math:`cm^2 s`.
The exposure cube does not have to contain the same binning as the event cube
but for simplicity, the event cube can be passed to adopt the binning
parameters. Note, however, that the exposure cube is defined in true sky
coordinates and energy while the counts cube is defined in reconstructed sky
coordinates and energy. Consequently, the sky area and energy range covered by
the exposure cube should be slightly larger than that of the counts cube to
accommodate for spill over of events due to the point spread function and energy
dispersion. By computing the exposure cube on the same grid as the counts cube,
the spill over of events from sources at the edge of cube will not be handled
correctly. In the example below, however, no source at the edge of the field of
view is present. Therefore, for simplicity, the count cube is used as input cube
to extract the binning.

This task of computing the exposure cube is is performed by :ref:`ctexpcube`. 

.. code-block:: bash

  $ ctexpcube
  Event list or observation definition file [NONE] selected_obs.xml
  Input counts cube file to extract exposure cube definition [NONE] cntcube.fits
  Output exposure cube file [expcube.fits]
  
Alternatively, the exposure cube can be created with different binning than the
counts cube:

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

Create point spread function cube
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As a next step for the stacked analysis, a cube containing the point spread
function (PSF) must be computed. Since the PSF cannot be stored by one single
parameter, the PSF cube computed by :ref:`ctpsfcube` has a fourth dimension.
In each bin of the cube, the PSF is stored as a function of offset from source.
The granularity of the PSF histogram is determined by the hidden parameter
``anumbins`` (default: 200). Therefore, when passing the event cube to adopt
the sky binning for the PSF cube, the resulting FITS file can become quite large
due to the fourth dimension. Usually in IACT analysis, the PSF doesn't change
too dramatically across the field of view. Therefore the user can think about
reducing the spatial binning of the PSF cube:

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
  
Depending on the required PSF precision, one could reduce the number of offset
bins via the hidden parameter ``anumbins``:

.. code-block:: bash

  $ ctpsfcube anumbins=100

Create energy dispersion cube (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In case energy dispersion should be considered used in the analysis, the tool
:ref:`ctedispcube` computes the energy migration response for the given
observations. Analoguous to the PSF cube, the energy dispersion cube is stored
in a four-dimensional map. In each bin of the cube, the migration value
``E_reco/E_true`` is stored. The granularity of the migration histogram can be
steered via the hidden parameter ``migrabins`` (default: 100). In addition, the
maximum migration value can be set via the hidden parameter ``migramax``
(default: 2.0). Similar to the PSF cube, the energy dispersion cube FITS file
can become quite large. Therefore, since the energy dispersion shouldn't vary
too much across the sky region, one could reduce the spatial binning:

.. code-block:: bash

  $ ctedispcube
  Input event list or observation definition XML file [NONE] selected_obs.xml 
  Input counts cube file to extract energy dispersion cube definition [NONE] 
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
  Output energy dispersion cube file [edispcube.fits]

Create background cube
^^^^^^^^^^^^^^^^^^^^^^
Last but not least, for a stacked IACT analysis a cube containing the background
rate in sky coordinates and reconstructed energy has to be computed. This task
is performed by :ref:`ctbkgcube`. The binning here can also differ from the counts
cube. For simplicity, however, the example below uses the counts cube to adopt
the binning.

.. code-block:: bash

  $ ctbkgcube debug=yes
  Input event list or observation definition XML file [NONE] selected_obs.xml 
  Input counts cube file to extract background cube definition [NONE] cntcube.fits 
  Input model XML file [NONE] crab_models.xml 
  Output background cube file [bkgcube.fits] 
  Output model XML file [NONE] binned_models.xml

Note that this tool also requires the parameters of an input and output model.
In the model XML file that came out of :ref:`csiactobs`, one background model
per observation is included. This models get merged and averaged in the background
sky cube. In the output model the background models per obseervation will be
removed. Instead, a global background model for the newly created background
cube is included. Sky models present in the input model XML file will also be
included in the new XML file, which subsequently can be used for stacked
:ref:`ctlike`.

Example for stacked model XML file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The output model of :ref:`ctbkgcube` looks the following:

.. code-block:: xml

	<?xml version="1.0" encoding="UTF-8" standalone="no"?>
	<source_library title="source library">
	  <source name="BackgroundModel" type="CTACubeBackground" instrument="CTA,HESS,MAGIC,VERITAS">
	    <spectrum type="PowerLaw">
	      <parameter name="Prefactor"   value="1" error="0" scale="1" min="0.01" max="100" free="1" />
	      <parameter name="Index"       value="0" error="0" scale="1" min="-5" max="5" free="1" />
	      <parameter name="PivotEnergy" value="1" scale="1e+06" free="0" />
	    </spectrum>
	  </source>
	  <source name="Crab" type="PointSource">
	    <spectrum type="PowerLaw">
	       <parameter name="Prefactor"   scale="1e-17" value="3.0"  min="1e-07" max="1000.0" free="1"/>
	       <parameter name="Index"       scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
	       <parameter name="PivotEnergy" scale="1e6"   value="1."  min="0.01"  max="1000.0" free="0"/>
	    </spectrum>
	    <spatialModel type="PointSource">
	      <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="0"/>
	      <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="0"/>
	    </spatialModel>
	  </source>
	</source_library>
	
The background model with ``type=CTACubeBackground`` is used to scale the
background cube stored in the FITS file created by :ref:`ctbkgcube`.

Run ctlike
^^^^^^^^^^
Having all the intermediate data products ready, a stacked analysis can be
conducted using :ref:`ctlike`.

.. code-block:: bash

  $ ctlike
  Input event list, counts cube or observation definition XML file [events.fits] cntcube.fits 
  Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits 
  Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits 
  Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits 
  Input model XML file [binned_models.xml]
  Output model XML file [binned_results.xml]
	
Note that when passing an counts cube to :ref:`ctlike`, the tool behaves
differently than in unbinned mode. It queries directly for the additional
ingredients for the stacked analysis. It is important to pass the background
model generated by :ref:`ctbkgcube` here to ensure the proper modelling of
the background in the fit.

To consider also the energy dispersion, you have to execute the :ref:`ctlike`
tool with the hidden parameter ``edisp=yes``. The tool will then also query for
the energy dispersion cube:

.. code-block:: bash
  
  $ ctlike debug=yes
  Input event list, counts cube or observation definition XML file [events.fits] cntcube.fits 
  Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits 
  Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits 
  Input energy dispersion cube file (only needed for stacked analysis) [NONE] edispcube.fits 
  Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits 
  Input model XML file [binned_models.xml]
  Output model XML file [binned_results.xml] 
