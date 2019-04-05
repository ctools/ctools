.. _sec_iact_selection:

Selecting IACT observations
===========================

  .. admonition:: What you will learn

     You will learn how to select observations from a dataset.

As next step you have to select observations from a dataset and to build
an observation definition XML file for the analysis.


Selecting observations
----------------------

First run :ref:`csfindobs` to select the observations of interest.
The most common search for observations is by pointing direction on the sky,
and you do this by typing

.. code-block:: bash

   $ csfindobs
   Path were data are located (NONE uses $VHEFITS environment variable) [NONE] data
   Name of FITS production (Run csiactdata to view your options) [prod-name] hess_dl3_dr1
   Right ascension of selection region centre (deg) [83.63]
   Declination of selection region centre (deg) [22.51]
   Search radius of selection region (deg) [2.5]
   Output runlist file [runlist.lis]

This creates an output ASCII file named ``runlist.lis`` that contains the
identifiers for all observations with pointing directions located within 2.5 deg
around the position of the Crab nebula. The content of the file is

.. code-block:: none

   23523
   23526
   23559
   23592

hence four observations were found in the production. An excerpt of the log
file ``csfindobs.log`` of the run is shown below:

.. code-block:: none

   2019-04-05T13:22:08: +===================+
   2019-04-05T13:22:08: | Find observations |
   2019-04-05T13:22:08: +===================+
   2019-04-05T13:22:08:  Expression ................: ANGSEP(83.63,22.51,RA_PNT,DEC_PNT)<=2.5&&QUALITY<=0
   2019-04-05T13:22:08:  Observations ..............: 4
   2019-04-05T13:22:08:  Observation 1 .............: 23523
   2019-04-05T13:22:08:  Observation 2 .............: 23526
   2019-04-05T13:22:08:  Observation 3 .............: 23559
   2019-04-05T13:22:08:  Observation 4 .............: 23592

.. note::
   You may add further selection criteria using the hidden ``expression``
   parameter.
   You may use any FITS selection expression that is
   `supported by the cfitsio row filtering specification <https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node97.html>`_.
   Column names that can be used in the expression can be found
   `here <http://gamma-astro-data-formats.readthedocs.org/en/latest/data_storage/obs_index/index.html>`_
   or by simply browsing the observation index file with a FITS viewer.
   For example

   .. code-block:: bash

      $ csfindobs expression="ZEN_PNT < 30 && LIVETIME > 1600"
      Path were data are located (NONE uses $VHEFITS environment variable) [data] 
      Name of FITS production (Run csiactdata to view your options) [hess_dl3_dr1]
      Right ascension of selection region centre (deg) [83.63] NONE
      Output runlist file [runlist.lis]

   will select all observations that have a zenith angle less than 30 deg and a
   livetime larger than 1600 seconds. Note that pointing selection was omitted
   by specifying ``NONE`` for the Right Ascension.

.. note::
   By default, :ref:`csfindobs` only selects data of highest quality
   (i.e. ``QUALITY=0``).
   You may overwrite this default by specifying the hidden parameter 
   ``min_qual``. For example, ``min_qual=1`` selects all data with a 
   looser quality criteria.


Create an observation list
--------------------------

Now you can convert the runlist into an observation definition XML file.
You do this conversion with the :ref:`csiactobs` script:

.. code-block:: bash

   $ csiactobs
   Path were data are located (NONE uses $VHEFITS environment variable) [NONE] data
   Data storage name [fits-prod-name] hess_dl3_dr1
   Input runlist file [runlist.lis]
   Number of free parameters per background model [1] 2
   Output model definition XML file [bkgmodels.xml]
   Output observation definition XML file [obs.xml]

The :ref:`csiactobs` script will create two output files: the observation
definition XML file ``obs.xml`` and an output model definition XML file
``bkgmodels.xml``. To generate ``obs.xml``, :ref:`csiactobs` has used the
IACT data storage and extracted the relevant file names. ``bgmodels.xml`` is
a file that is used for background modeling, where each observation will have
its own independent background model. In the example above, you have set the
number of free parameters per background model to one, hence the normalisation
of the background model for each observation will be a free parameter that is
later adjusted by the maximum likelihood fit.

Here the observation definition file that was created

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <observation_list title="observation list">
     <observation name="Crab Nebula" id="23523" instrument="HESS">
       <parameter name="EventList" file="data/data/hess_dl3_dr1_obs_id_023523.fits.gz[events]" />
       <parameter name="EffectiveArea" file="data/data/hess_dl3_dr1_obs_id_023523.fits.gz[aeff]" />
       <parameter name="PointSpreadFunction" file="data/data/hess_dl3_dr1_obs_id_023523.fits.gz[psf]" />
       <parameter name="EnergyDispersion" file="data/data/hess_dl3_dr1_obs_id_023523.fits.gz[edisp]" />
       <parameter name="Background" file="" />
     </observation>
     <observation name="Crab Nebula" id="23526" instrument="HESS">
       <parameter name="EventList" file="data/data/hess_dl3_dr1_obs_id_023526.fits.gz[events]" />
       <parameter name="EffectiveArea" file="data/data/hess_dl3_dr1_obs_id_023526.fits.gz[aeff]" />
       <parameter name="PointSpreadFunction" file="data/data/hess_dl3_dr1_obs_id_023526.fits.gz[psf]" />
       <parameter name="EnergyDispersion" file="data/data/hess_dl3_dr1_obs_id_023526.fits.gz[edisp]" />
       <parameter name="Background" file="" />
     </observation>
     <observation name="Crab Nebula" id="23559" instrument="HESS">
       <parameter name="EventList" file="data/data/hess_dl3_dr1_obs_id_023559.fits.gz[events]" />
       <parameter name="EffectiveArea" file="data/data/hess_dl3_dr1_obs_id_023559.fits.gz[aeff]" />
       <parameter name="PointSpreadFunction" file="data/data/hess_dl3_dr1_obs_id_023559.fits.gz[psf]" />
       <parameter name="EnergyDispersion" file="data/data/hess_dl3_dr1_obs_id_023559.fits.gz[edisp]" />
       <parameter name="Background" file="" />
     </observation>
     <observation name="Crab Nebula" id="23592" instrument="HESS">
       <parameter name="EventList" file="data/data/hess_dl3_dr1_obs_id_023592.fits.gz[events]" />
       <parameter name="EffectiveArea" file="data/data/hess_dl3_dr1_obs_id_023592.fits.gz[aeff]" />
       <parameter name="PointSpreadFunction" file="data/data/hess_dl3_dr1_obs_id_023592.fits.gz[psf]" />
       <parameter name="EnergyDispersion" file="data/data/hess_dl3_dr1_obs_id_023592.fits.gz[edisp]" />
       <parameter name="Background" file="" />
     </observation>
   </observation_list>
   <observation_list title="observation list" />

and here the model definition file

.. code-block:: xml

  <?xml version="1.0" encoding="UTF-8" standalone="no"?>
  <source_library title="source library">
    <source name="bkg_23523" type="CTAAeffBackground" instrument="HESS" id="23523" tscalc="0">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor" value="1" error="0" scale="1e-14" min="0.01" max="100" free="1" />
        <parameter name="Index" value="-2" scale="1" min="-5" max="5" free="1" />
        <parameter name="PivotEnergy" value="1" scale="1000000" free="0" />
      </spectrum>
    </source>
    <source name="bkg_23526" type="CTAAeffBackground" instrument="HESS" id="23526" tscalc="0">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor" value="1" error="0" scale="1e-14" min="0.01" max="100" free="1" />
        <parameter name="Index" value="-2" scale="1" min="-5" max="5" free="1" />
        <parameter name="PivotEnergy" value="1" scale="1000000" free="0" />
      </spectrum>
    </source>
    <source name="bkg_23559" type="CTAAeffBackground" instrument="HESS" id="23559" tscalc="0">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor" value="1" error="0" scale="1e-14" min="0.01" max="100" free="1" />
        <parameter name="Index" value="-2" scale="1" min="-5" max="5" free="1" />
        <parameter name="PivotEnergy" value="1" scale="1000000" free="0" />
      </spectrum>
    </source>
    <source name="bkg_23592" type="CTAAeffBackground" instrument="HESS" id="23592" tscalc="0">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor" value="1" error="0" scale="1e-14" min="0.01" max="100" free="1" />
        <parameter name="Index" value="-2" scale="1" min="-5" max="5" free="1" />
        <parameter name="PivotEnergy" value="1" scale="1000000" free="0" />
      </spectrum>
    </source>
  </source_library>

.. note::
   In case that you have already a model definition XML file that describes a
   model of the celestial source distribution (a.k.a. a sky model), you may
   provide this sky model to :ref:`csiactobs` using the hidden ``inmodel``
   parameter. The output model definition XML file will then contain the sky
   model and the background model, and you can use the model directly in a
   maximum likelihood fit.

   .. code-block:: bash

      $ csiactobs inmodel=$CTOOLS/share/models/crab.xml
      Path were data are located (NONE uses $VHEFITS environment variable) [data]
      Data storage name [hess_dl3_dr1]
      Input runlist file [runlist.lis]
      Number of free parameters per background model [2]
      Output model definition XML file [bkgmodels.xml] models.xml
      Output observation definition XML file [obs.xml]


Select events
-------------

The final step needed is the selection of events, to define in particular
the region of interest and energy range for the analysis.
You do this by typing

.. code-block:: bash

   $ ctselect
   Input event list or observation definition XML file [events.fits] obs.xml
   Radius of ROI around pointing or specified RA/DEC (degrees) (0-180) [3.0] 2.0
   Start time (UTC string, JD, MJD or MET in seconds) [NONE]
   Lower energy limit (TeV) [0.1] 0.5
   Upper energy limit (TeV) [100.0] 10.0
   Output event list or observation definition XML file [selected_events.fits] obs_selected.xml

which will select all events within a radius of 2 degrees around the pointing
direction with energies between 500 GeV and 10 TeV.

