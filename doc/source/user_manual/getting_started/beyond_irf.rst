.. _sec_connecting_irf:

Connecting observations to dedicated instrument response functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

So far the same instrument response function has been used for all
observations in our analysis examples.
You may have recognised that we always specified the calibration database
``prod2`` and the instrument response function ``South_50h`` that
are shipped with ctools.
When combining multiple observations in a joint analysis we may however
need to specify dedicated instrument response functions for the
different observations.

Let's simulate two observations of the Crab nebula, one done with
the northern array and the other with the southern array:

.. code-block:: bash

  $ ctobssim
  RA of pointing (degrees) (0-360) [83.63] 
  Dec of pointing (degrees) (-90-90) [22.01] 21.51
  Radius of FOV (degrees) (0-180) [5.0] 
  Start time (MET in s) [0.0] 
  End time (MET in s) [1800.0] 
  Lower energy limit (TeV) [0.1] 
  Upper energy limit (TeV) [100.0] 
  Calibration database [prod2] 
  Instrument response function [South_50h] North_50h
  Model [$CTOOLS/share/models/crab.xml] 
  Output event data file or observation definition file [events.fits] events1.fits

.. code-block:: bash

  $ ctobssim
  RA of pointing (degrees) (0-360) [83.63] 
  Dec of pointing (degrees) (-90-90) [22.01] 22.51
  Radius of FOV (degrees) (0-180) [5.0] 
  Start time (MET in s) [0.0] 
  End time (MET in s) [1800.0] 
  Lower energy limit (TeV) [0.1] 
  Upper energy limit (TeV) [100.0] 
  Calibration database [prod2] 
  Instrument response function [North_50h] South_50h
  Model [$CTOOLS/share/models/crab.xml] 
  Output event data file or observation definition file [events1.fits] events2.fits

Now we combine the information about these two observations in an
observation definition XML file:

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <observation_list title="observation library">
    <observation name="Crab" id="00001" instrument="CTA">
      <parameter name="EventList" file="events1.fits"/>
      <parameter name="Calibration" database="prod2" response="North_50h"/>
    </observation>
    <observation name="Crab" id="00002" instrument="CTA">
      <parameter name="EventList" file="events2.fits"/>
      <parameter name="Calibration" database="prod2" response="South_50h"/>
    </observation>
  </observation_list>

Each observation now has a ``Calibration`` parameter with the attributes
``database`` and ``response``.
To see which calibration databases and response functions are available on
your system you can run the :ref:`cscaldb` script (use the ``debug=yes``
option to get the output of the script on the console):

.. code-block:: bash

   $ cscaldb debug=yes
   Parfile cscaldb.par not found. Create default parfile.
   2015-05-22T21:35:39: +============+
   2015-05-22T21:35:39: | Parameters |
   2015-05-22T21:35:39: +============+
   2015-05-22T21:35:39:  chatter ...................: 2
   2015-05-22T21:35:39:  clobber ...................: yes
   2015-05-22T21:35:39:  debug .....................: yes
   2015-05-22T21:35:39:  mode ......................: ql
   2015-05-22T21:35:39: 
   2015-05-22T21:35:39: +==============+
   2015-05-22T21:35:39: | Mission: cta |
   2015-05-22T21:35:39: +==============+
   2015-05-22T21:35:40: === Response functions in database "prod2" ===
   2015-05-22T21:35:40: North_0.5h
   2015-05-22T21:35:40: North_50h
   2015-05-22T21:35:40: North_5h
   2015-05-22T21:35:40: South_0.5h
   2015-05-22T21:35:40: South_50h
   2015-05-22T21:35:40: South_5h

You will see that ctools ships with one database for the CTA observatory.
This is the ``prod2`` database.
Within this database there are six response functions:
``North_0.5h``, ``North_5h``, ``North_50h``,
``South_0.5h``, ``South_5h``, and ``South_50h``.

We now jointly fit both observations using :ref:`ctlike`:

.. code-block:: bash

   $ ctlike
   Event list, counts cube or observation definition file [obs.xml] obs_irf.xml
   Source model [$CTOOLS/share/models/crab.xml] 
   Source model output file [crab_results.xml] 

To see the usage of the individual response functions you may inspect the
log file. You will notice that each observation now has a specific response
function and that the filename of the response information differs for
both observations.

.. code-block:: xml

  2015-05-22T22:06:43: +==============+
  2015-05-22T22:06:43: | Observations |
  2015-05-22T22:06:43: +==============+
  2015-05-22T22:06:43: === GObservations ===
  2015-05-22T22:06:43:  Number of observations ....: 2
  2015-05-22T22:06:43:  Number of predicted events : 0
  2015-05-22T22:06:43: === GCTAObservation ===
  2015-05-22T22:06:43:  Name ......................: 
  2015-05-22T22:06:43:  Identifier ................: 00001
  2015-05-22T22:06:43:  Instrument ................: CTA
  2015-05-22T22:06:43:  Event file ................: events1.fits
  2015-05-22T22:06:43:  Event type ................: EventList
  2015-05-22T22:06:43:  Statistics ................: Poisson
  2015-05-22T22:06:43:  Ontime ....................: 1800 s
  2015-05-22T22:06:43:  Livetime ..................: 1710 s
  2015-05-22T22:06:43:  Deadtime correction .......: 0.95
  2015-05-22T22:06:43:  User energy range .........: undefined
  2015-05-22T22:06:43: === GCTAPointing ===
  2015-05-22T22:06:43:  Pointing direction ........: (RA,Dec)=(83.63,21.51)
  2015-05-22T22:06:43: === GCTAResponseIrf ===
  2015-05-22T22:06:43:  Caldb mission .............: cta
  2015-05-22T22:06:43:  Caldb instrument ..........: prod2
  2015-05-22T22:06:43:  Response name .............: North_50h
  2015-05-22T22:06:43:  Energy dispersion .........: Not used
  2015-05-22T22:06:43:  Save energy range .........: undefined
  ...
  2015-05-22T22:06:43: === GCTAObservation ===
  2015-05-22T22:06:43:  Name ......................: 
  2015-05-22T22:06:43:  Identifier ................: 00002
  2015-05-22T22:06:43:  Instrument ................: CTA
  2015-05-22T22:06:43:  Event file ................: events2.fits
  2015-05-22T22:06:43:  Event type ................: EventList
  2015-05-22T22:06:43:  Statistics ................: Poisson
  2015-05-22T22:06:43:  Ontime ....................: 1800 s
  2015-05-22T22:06:43:  Livetime ..................: 1710 s
  2015-05-22T22:06:43:  Deadtime correction .......: 0.95
  2015-05-22T22:06:43:  User energy range .........: undefined
  2015-05-22T22:06:43: === GCTAPointing ===
  2015-05-22T22:06:43:  Pointing direction ........: (RA,Dec)=(83.63,22.51)
  2015-05-22T22:06:43: === GCTAResponseIrf ===
  2015-05-22T22:06:43:  Caldb mission .............: cta
  2015-05-22T22:06:43:  Caldb instrument ..........: prod2
  2015-05-22T22:06:43:  Response name .............: South_50h
  2015-05-22T22:06:43:  Energy dispersion .........: Not used
  2015-05-22T22:06:43:  Save energy range .........: undefined
  ...
  2015-05-22T22:06:43: +=================================+
  2015-05-22T22:06:43: | Maximum likelihood optimisation |
  2015-05-22T22:06:43: +=================================+
  2015-05-22T22:06:43:  >Iteration   0: -logL=138039.472, Lambda=1.0e-03
  2015-05-22T22:06:43:  >Iteration   1: -logL=138035.501, Lambda=1.0e-03, delta=3.971, max(|grad|)=16.203761 [Index:7]
  2015-05-22T22:06:43:  >Iteration   2: -logL=138035.496, Lambda=1.0e-04, delta=0.005, max(|grad|)=-0.053284 [Index:7]
  2015-05-22T22:06:43:  >Iteration   3: -logL=138035.496, Lambda=1.0e-05, delta=0.000, max(|grad|)=0.001386 [Index:7]
  2015-05-22T22:06:44: 
  2015-05-22T22:06:44: +=========================================+
  2015-05-22T22:06:44: | Maximum likelihood optimization results |
  2015-05-22T22:06:44: +=========================================+
  2015-05-22T22:06:44: === GOptimizerLM ===
  2015-05-22T22:06:44:  Optimized function value ..: 138035.496
  2015-05-22T22:06:44:  Absolute precision ........: 0.005
  2015-05-22T22:06:44:  Acceptable value decrease .: 2
  2015-05-22T22:06:44:  Optimization status .......: converged
  2015-05-22T22:06:44:  Number of parameters ......: 10
  2015-05-22T22:06:44:  Number of free parameters .: 4
  2015-05-22T22:06:44:  Number of iterations ......: 3
  2015-05-22T22:06:44:  Lambda ....................: 1e-06
  2015-05-22T22:06:44:  Maximum log likelihood ....: -138035.496
  2015-05-22T22:06:44:  Observed events  (Nobs) ...: 20814.000
  2015-05-22T22:06:44:  Predicted events (Npred) ..: 20814.000 (Nobs - Npred = 2.67698e-06)
  2015-05-22T22:06:44: === GModels ===
  2015-05-22T22:06:44:  Number of models ..........: 2
  2015-05-22T22:06:44:  Number of parameters ......: 10
  2015-05-22T22:06:44: === GModelSky ===
  2015-05-22T22:06:44:  Name ......................: Crab
  2015-05-22T22:06:44:  Instruments ...............: all
  2015-05-22T22:06:44:  Instrument scale factors ..: unity
  2015-05-22T22:06:44:  Observation identifiers ...: all
  2015-05-22T22:06:44:  Model type ................: PointSource
  2015-05-22T22:06:44:  Model components ..........: "SkyDirFunction" * "PowerLaw" * "Constant"
  2015-05-22T22:06:44:  Number of parameters ......: 6
  2015-05-22T22:06:44:  Number of spatial par's ...: 2
  2015-05-22T22:06:44:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
  2015-05-22T22:06:44:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
  2015-05-22T22:06:44:  Number of spectral par's ..: 3
  2015-05-22T22:06:44:   Prefactor ................: 5.9184e-16 +/- 9.68924e-18 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
  2015-05-22T22:06:44:   Index ....................: -2.49569 +/- 0.0145376 [-0,-5]  (free,scale=-1,gradient)
  2015-05-22T22:06:44:   PivotEnergy ..............: 300000 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2015-05-22T22:06:44:  Number of temporal par's ..: 1
  2015-05-22T22:06:44:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
  2015-05-22T22:06:44: === GCTAModelIrfBackground ===
  2015-05-22T22:06:44:  Name ......................: CTABackgroundModel
  2015-05-22T22:06:44:  Instruments ...............: CTA
  2015-05-22T22:06:44:  Instrument scale factors ..: unity
  2015-05-22T22:06:44:  Observation identifiers ...: all
  2015-05-22T22:06:44:  Model type ................: "PowerLaw" * "Constant"
  2015-05-22T22:06:44:  Number of parameters ......: 4
  2015-05-22T22:06:44:  Number of spectral par's ..: 3
  2015-05-22T22:06:44:   Prefactor ................: 1.02524 +/- 0.0182054 [0.001,1000] ph/cm2/s/MeV (free,scale=1,gradient)
  2015-05-22T22:06:44:   Index ....................: 0.0229728 +/- 0.0104713 [-5,5]  (free,scale=1,gradient)
  2015-05-22T22:06:44:   PivotEnergy ..............: 1e+06 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2015-05-22T22:06:44:  Number of temporal par's ..: 1
  2015-05-22T22:06:44:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)



A more fine grained control over the response function can be achieved by 
specifying individual filenames for the various response components.
An example for an observation definition XML file is shown below.
This is definitely expert mode, to be used with utmost care.

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <observation_list title="observation library">
    <observation name="Crab" id="00001" instrument="CTA">
      <parameter name="EventList"           file="events1.fits"/>
      <parameter name="EffectiveArea"       file="$(CALDB)/data/cta/prod2/bcf/North_50h/irf_file.fits.gz"/>
      <parameter name="PointSpreadFunction" file="$(CALDB)/data/cta/prod2/bcf/North_50h/irf_file.fits.gz"/>
      <parameter name="EnergyDispersion"    file="$(CALDB)/data/cta/prod2/bcf/North_50h/irf_file.fits.gz"/>
      <parameter name="Background"          file="$(CALDB)/data/cta/prod2/bcf/North_50h/irf_file.fits.gz"/>
    </observation>
    <observation name="Crab" id="00002" instrument="CTA">
      <parameter name="EventList"           file="events2.fits"/>
      <parameter name="EffectiveArea"       file="$(CALDB)/data/cta/prod2/bcf/South_50h/irf_file.fits.gz"/>
      <parameter name="PointSpreadFunction" file="$(CALDB)/data/cta/prod2/bcf/South_50h/irf_file.fits.gz"/>
      <parameter name="EnergyDispersion"    file="$(CALDB)/data/cta/prod2/bcf/South_50h/irf_file.fits.gz"/>
      <parameter name="Background"          file="$(CALDB)/data/cta/prod2/bcf/South_50h/irf_file.fits.gz"/>
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
      <parameter name="BkgCube"      file="bkgcube1.fits"/>
    </observation>
    <observation name="Crab" id="00002" instrument="CTA">
      <parameter name="CountsCube"   file="cntcube2.fits"/>
      <parameter name="ExposureCube" file="expcube2.fits"/>
      <parameter name="PsfCube"      file="psfcube2.fits"/>
      <parameter name="BkgCube"      file="bkgcube2.fits"/>
    </observation>
  </observation_list>


