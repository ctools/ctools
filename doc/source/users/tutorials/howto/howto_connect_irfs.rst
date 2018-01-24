.. _howto_connect_irfs:

How to connect observations to dedicated instrument response functions?
-----------------------------------------------------------------------

  .. admonition:: What you will learn

     You will learn how to **connect individual observations to dedicated
     Instrument Response Functions**.


So far you used the same instrument response function for all observations.
You always specified the calibration database ``prod2`` and the instrument
response function ``South_0.5h`` that are shipped with ctools in all
examples.
When combining multiple observations in a joint analysis you may however
need to specify dedicated instrument response functions for the
different observations.

To illustrate how you can do this, let's simulate two observations of the
Crab nebula, one done with the northern array and the other with the southern
array:

.. code-block:: bash

   $ ctobssim
   RA of pointing (degrees) (0-360) [83.63]
   Dec of pointing (degrees) (-90-90) [22.51] 21.51
   Radius of FOV (degrees) (0-180) [5.0]
   Start time (UTC string, JD, MJD or MET in seconds) [2020-01-01T00:00:00]
   Stop time (UTC string, JD, MJD or MET in seconds) [2020-01-01T00:30:00]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Calibration database [prod2]
   Instrument response function [South_0.5h] North_0.5h
   Input model definition XML file [$CTOOLS/share/models/crab.xml]
   Output event data file or observation definition XML file [events.fits] north.fits

.. code-block:: bash

   $ ctobssim
   RA of pointing (degrees) (0-360) [83.63]
   Dec of pointing (degrees) (-90-90) [21.51] 22.51
   Radius of FOV (degrees) (0-180) [5.0]
   Start time (UTC string, JD, MJD or MET in seconds) [2020-01-01T00:00:00] 2020-01-01T00:30:00
   Stop time (UTC string, JD, MJD or MET in seconds) [2020-01-01T00:30:00] 2020-01-01T01:00:00
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Calibration database [prod2]
   Instrument response function [North_0.5h] South_0.5h
   Input model definition XML file [$CTOOLS/share/models/crab.xml]
   Output event data file or observation definition XML file [north.fits] south.fits

You now need to combine the information about these two observations in an
observation definition XML file:

.. code-block:: xml

   <?xml version="1.0" standalone="no"?>
   <observation_list title="observation library">
     <observation name="Crab" id="00001" instrument="CTA">
       <parameter name="EventList" file="north.fits"/>
       <parameter name="Calibration" database="prod2" response="North_0.5h"/>
     </observation>
     <observation name="Crab" id="00002" instrument="CTA">
       <parameter name="EventList" file="south.fits"/>
       <parameter name="Calibration" database="prod2" response="South_0.5h"/>
     </observation>
  </observation_list>

Each observation now has a ``Calibration`` parameter with the attributes
``database`` and ``response``.
To see which calibration databases and response functions are available on
your system you can run the :ref:`cscaldb` script (use the ``debug=yes``
option to get the output of the script on the console):

.. code-block:: bash

   $ cscaldb debug=yes
   2018-01-24T16:12:35: +============+
   2018-01-24T16:12:35: | Parameters |
   2018-01-24T16:12:35: +============+
   2018-01-24T16:12:35:  chatter ...................: 2
   2018-01-24T16:12:35:  clobber ...................: yes
   2018-01-24T16:12:35:  debug .....................: yes
   2018-01-24T16:12:35:  mode ......................: ql
   2018-01-24T16:12:35:  logfile ...................: cscaldb.log
   2018-01-24T16:12:35:
   2018-01-24T16:12:35: +==============+
   2018-01-24T16:12:35: | Mission: cta |
   2018-01-24T16:12:35: +==============+
   2018-01-24T16:12:35: === Response functions in database "prod2" ===
   2018-01-24T16:12:35: North_0.5h
   2018-01-24T16:12:35: North_50h
   2018-01-24T16:12:35: North_5h
   2018-01-24T16:12:35: South_0.5h
   2018-01-24T16:12:35: South_50h
   2018-01-24T16:12:35: South_5h

You will see that ctools ships with one database for the CTA observatory.
This is the ``prod2`` database.
Within this database there are six response functions:
``North_0.5h``, ``North_5h``, ``North_50h``,
``South_0.5h``, ``South_5h``, and ``South_50h``.

You now jointly fit both observations using :ref:`ctlike`:

.. code-block:: bash

   $ ctlike chatter=3
   Input event list, counts cube or observation definition XML file [events.fits] obs_irf.xml
   Input model definition XML file [$CTOOLS/share/models/crab.xml]
   Output model definition XML file [crab_results.xml]

To see the usage of the individual response functions you may inspect the
log file (you need to set ``chatter=3`` to see the details of the 
observations in the log file).
You will notice that each observation now has a specific response
function and that the filename of the response information differs for
both observations.

.. code-block:: none

   2018-01-24T16:13:12: +====================+
   2018-01-24T16:13:12: | Input observations |
   2018-01-24T16:13:12: +====================+
   2018-01-24T16:13:12: === GObservations ===
   2018-01-24T16:13:12:  Number of observations ....: 2
   2018-01-24T16:13:12:  Number of models ..........: 2
   2018-01-24T16:13:12:  Number of observed events .: 32097
   2018-01-24T16:13:12:  Number of predicted events : 0
   2018-01-24T16:13:12: === GCTAObservation ===
   2018-01-24T16:13:12:  Name ......................: Crab
   2018-01-24T16:13:12:  Identifier ................: 00001
   ...
   2018-01-24T16:13:12: === GCTAResponseIrf ===
   2018-01-24T16:13:12:  Caldb mission .............: cta
   2018-01-24T16:13:12:  Caldb instrument ..........: prod2
   2018-01-24T16:13:12:  Response name .............: North_0.5h
   2018-01-24T16:13:12:  Energy dispersion .........: Not used
   2018-01-24T16:13:12:  Save energy range .........: undefined
   ...
   2018-01-24T16:13:12: === GCTAObservation ===
   2018-01-24T16:13:12:  Name ......................: Crab
   2018-01-24T16:13:12:  Identifier ................: 00002
   ...
   2018-01-24T16:13:12: === GCTAResponseIrf ===
   2018-01-24T16:13:12:  Caldb mission .............: cta
   2018-01-24T16:13:12:  Caldb instrument ..........: prod2
   2018-01-24T16:13:12:  Response name .............: South_0.5h
   2018-01-24T16:13:12:  Energy dispersion .........: Not used
   2018-01-24T16:13:12:  Save energy range .........: undefined
   ...
   2018-01-24T16:13:12: +=================================+
   2018-01-24T16:13:12: | Maximum likelihood optimisation |
   2018-01-24T16:13:12: +=================================+
   2018-01-24T16:13:12:  >Iteration   0: -logL=220368.120, Lambda=1.0e-03
   2018-01-24T16:13:12:  >Iteration   1: -logL=220367.708, Lambda=1.0e-03, delta=0.412, step=1.0e+00, max(|grad|)=-0.544370 [Index:3]
   2018-01-24T16:13:12:  >Iteration   2: -logL=220367.708, Lambda=1.0e-04, delta=0.000, step=1.0e+00, max(|grad|)=0.013339 [Index:3]
   ...
   2018-01-24T16:13:12: +=========================================+
   2018-01-24T16:13:12: | Maximum likelihood optimisation results |
   2018-01-24T16:13:12: +=========================================+
   2018-01-24T16:13:12: === GOptimizerLM ===
   2018-01-24T16:13:12:  Optimized function value ..: 220367.708
   2018-01-24T16:13:12:  Absolute precision ........: 0.005
   2018-01-24T16:13:12:  Acceptable value decrease .: 2
   2018-01-24T16:13:12:  Optimization status .......: converged
   2018-01-24T16:13:12:  Number of parameters ......: 10
   2018-01-24T16:13:12:  Number of free parameters .: 4
   2018-01-24T16:13:12:  Number of iterations ......: 2
   2018-01-24T16:13:12:  Lambda ....................: 1e-05
   2018-01-24T16:13:12:  Maximum log likelihood ....: -220367.708
   2018-01-24T16:13:12:  Observed events  (Nobs) ...: 32097.000
   2018-01-24T16:13:12:  Predicted events (Npred) ..: 32097.000 (Nobs - Npred = 0.000298515627946472)
   2018-01-24T16:13:12: === GModels ===
   2018-01-24T16:13:12:  Number of models ..........: 2
   2018-01-24T16:13:12:  Number of parameters ......: 10
   2018-01-24T16:13:12: === GModelSky ===
   2018-01-24T16:13:12:  Name ......................: Crab
   2018-01-24T16:13:12:  Instruments ...............: all
   2018-01-24T16:13:12:  Instrument scale factors ..: unity
   2018-01-24T16:13:12:  Observation identifiers ...: all
   2018-01-24T16:13:12:  Model type ................: PointSource
   2018-01-24T16:13:12:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
   2018-01-24T16:13:12:  Number of parameters ......: 6
   2018-01-24T16:13:12:  Number of spatial par's ...: 2
   2018-01-24T16:13:12:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
   2018-01-24T16:13:12:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
   2018-01-24T16:13:12:  Number of spectral par's ..: 3
   2018-01-24T16:13:12:   Prefactor ................: 5.69752531012742e-16 +/- 8.48707814807091e-18 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
   2018-01-24T16:13:12:   Index ....................: -2.46938547601713 +/- 0.0131947771405262 [-0,-5]  (free,scale=-1,gradient)
   2018-01-24T16:13:12:   PivotEnergy ..............: 300000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
   2018-01-24T16:13:12:  Number of temporal par's ..: 1
   2018-01-24T16:13:12:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2018-01-24T16:13:12: === GCTAModelIrfBackground ===
   2018-01-24T16:13:12:  Name ......................: CTABackgroundModel
   2018-01-24T16:13:12:  Instruments ...............: CTA
   2018-01-24T16:13:12:  Instrument scale factors ..: unity
   2018-01-24T16:13:12:  Observation identifiers ...: all
   2018-01-24T16:13:12:  Model type ................: "PowerLaw" * "Constant"
   2018-01-24T16:13:12:  Number of parameters ......: 4
   2018-01-24T16:13:12:  Number of spectral par's ..: 3
   2018-01-24T16:13:12:   Prefactor ................: 1.00359239305231 +/- 0.0102298598204261 [0.001,1000] ph/cm2/s/MeV (free,scale=1,gradient)
   2018-01-24T16:13:12:   Index ....................: 0.00205626275645075 +/- 0.0063475481805347 [-5,5]  (free,scale=1,gradient)
   2018-01-24T16:13:12:   PivotEnergy ..............: 1000000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
   2018-01-24T16:13:12:  Number of temporal par's ..: 1
   2018-01-24T16:13:12:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

You can have a more fine grained control over the response function by
specifying individual filenames for the various response components.
An example for an observation definition XML file is shown below.
This is definitely expert mode, to be used with utmost care.

.. code-block:: xml

   <?xml version="1.0" standalone="no"?>
   <observation_list title="observation library">
     <observation name="Crab" id="00001" instrument="CTA">
       <parameter name="EventList"           file="north.fits"/>
       <parameter name="EffectiveArea"       file="$(CALDB)/data/cta/prod2/bcf/North_0.5h/irf_file.fits.gz"/>
       <parameter name="PointSpreadFunction" file="$(CALDB)/data/cta/prod2/bcf/North_0.5h/irf_file.fits.gz"/>
       <parameter name="EnergyDispersion"    file="$(CALDB)/data/cta/prod2/bcf/North_0.5h/irf_file.fits.gz"/>
       <parameter name="Background"          file="$(CALDB)/data/cta/prod2/bcf/North_0.5h/irf_file.fits.gz"/>
     </observation>
     <observation name="Crab" id="00002" instrument="CTA">
       <parameter name="EventList"           file="south.fits"/>
       <parameter name="EffectiveArea"       file="$(CALDB)/data/cta/prod2/bcf/South_0.5h/irf_file.fits.gz"/>
       <parameter name="PointSpreadFunction" file="$(CALDB)/data/cta/prod2/bcf/South_0.5h/irf_file.fits.gz"/>
       <parameter name="EnergyDispersion"    file="$(CALDB)/data/cta/prod2/bcf/South_0.5h/irf_file.fits.gz"/>
       <parameter name="Background"          file="$(CALDB)/data/cta/prod2/bcf/South_0.5h/irf_file.fits.gz"/>
     </observation>
   </observation_list>

Finally, response information may also be provided to combine stacked
observations. An example for the syntax of the observation definition XML 
file is given below:

.. code-block:: xml

   <?xml version="1.0" standalone="no"?>
   <observation_list title="observation library">
     <observation name="Crab" id="00001" instrument="CTA">
       <parameter name="CountsCube"   file="cntcube1.fits"/>
       <parameter name="ExposureCube" file="expcube1.fits"/>
       <parameter name="PsfCube"      file="psfcube1.fits"/>
       <parameter name="EdispCube"    file="edispcube1.fits"/>
       <parameter name="BkgCube"      file="bkgcube1.fits"/>
     </observation>
     <observation name="Crab" id="00002" instrument="CTA">
       <parameter name="CountsCube"   file="cntcube2.fits"/>
       <parameter name="ExposureCube" file="expcube2.fits"/>
       <parameter name="PsfCube"      file="psfcube2.fits"/>
       <parameter name="EdispCube"    file="edispcube2.fits"/>
       <parameter name="BkgCube"      file="bkgcube2.fits"/>
     </observation>
   </observation_list>
