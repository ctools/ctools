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
  
This will list names of available FITS data sets on the screen. In the following it is assumed the users' choice
of the data set is called ``fits-data-name``. The query for the path where the data is located can be
omitted by setting environment variable ``$VHEFITS`` to this locations. This might be convenient
since some other tools also query for the same parameter and can use ``$VHEFITS`` instead.

Find observations
-----------------
To start an analysis, it is required to assemble a list of observations that should be used. For this purpose,
:ref:`csfindobs` searches for observations according to the user requirement.

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
  
  
Specifying additional search requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The search for observations can be further constrained arbitrarily. For this purpose, :ref:`csfindobs` contains a hidden parameter
``expression``. It is possible to use any FITS selection expression that is `supported by cfitsio <http://www.isdc.unige.ch/integral/download/osa/doc/10.1/osa_um_intro/node38.html>`_.
Available properties (i.e. column names) for selection, can be found `here <http://gamma-astro-data-formats.readthedocs.org/en/latest/data_storage/obs_index/index.html>`_ or by simply browsing the observation index file with a FITS viewer.

Example:

.. code-block:: bash

  $ csfindobs expression="ZEN_PNT<30&&LIVETIME>1800"
  ...
  
This command selects all observations that have a zenith angle less than 30 degrees and a livetime larger than 30 minutes.

Selections that don't include sky coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The sky selection can also be omitted. For instance one does not always know the sky coordinate of a certain object.
Therefore, it is possible to search for the object name instead.
For this purpose one can simply pass an invalid coordinate to :ref:`csfindobs`.

.. code-block:: bash

  $ csfindobs expression="OBJECT=='Crab Nebula'"
  Name of FITS production (Run csiactdata to view your options) [fits-data-name] 
  Right ascension [NONE]
  Runlist outfile [runlist.lis]

The selection by region will simply be omitted.

Note
^^^^
On default, :ref:`csfindobs` only select data of highest quality (i.e. QUALITY=0). Specifying the hidden parameter e.g. ``min_qual=1``
allows to select all data with looser quality criteria.

Create an observation list
--------------------------
The runlist ASCII file containing a list of selected observation IDs must now be converted to an
observation XML file. This file contains information about the location of the files that are required
for the analysis. The tool :ref:`csiactobs` is intended to do this conversion.

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
* In ``csiactobs.log`` (or on screen if ``debug=yes``), the script dumps the complete energy range of the observation container. These values might be important for later usage (e.g. in binned analysis).

In case a sky model is already prepared, it is possible to also provide the hidden parameter ``inmodel``. The output
model XML file will then contain both, the background model and the input sky model:

.. code-block:: bash

  $ csiactobs inmodel="mymodel.xml"
  
A list of available sky models can be found here: `models <http://gammalib.sourceforge.net/user_manual/modules/model.html>`_

Alternatively, models can be merged at any times using the simple tool :ref:`csmodelmerge`:

.. code-block:: bash

  $ csmodelmerge
  Input model XML files [bgmodels.xml crab.xml]
  Output model file [crab_models.xml]
  
Note that the number of files to merge is not limited to two. Detailled options how the input model XML file can
be passed is given on the reference page of :ref:`csmodelmerge`.

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
seems reasonable. The time selection is not applied in the above example since an invalid time range was provided.
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
  
The result of the fit was stored in 'crab_results.xml'. Note that fitted parameters, ``Prefactors`` in particular,
typically use MeV as energy unit. To monitor the progress of the fit on the screen, one can simply run with the option ``debug=yes``.
Alternatively, the logfile ``ctlike.log`` can be inspected after the fit. 

On default, energy dispersion is not considered in the fit. To switch on the usage of the energy migration matrix,
the hidden parameter ``edisp=yes`` can be provided. Note that this will cause a significant reduction of the computing
speed.

Binned analysis
---------------
For binned analysis, some intermediated data products have to be produced. The products are a binned data cube,
an exposure cube, a psf cube, and a background cube. 

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

Note that bins only get filled if bin of the cube is fully contained in the energy range and RoI of a considered observation.
It is therefore useful to provide the energy range given by :ref:`csiactobs` above. This ensures a maximum agreement between
observations and binning and reduces the loss of data.

Create exposure cube
^^^^^^^^^^^^^^^^^^^^

Create psf cube
^^^^^^^^^^^^^^^

Create background cube
^^^^^^^^^^^^^^^^^^^^^^

Run ctlike
^^^^^^^^^^

Compute spectral points
-----------------------

