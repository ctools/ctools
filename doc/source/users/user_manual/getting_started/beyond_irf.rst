.. _sec_connecting_irf:

Connecting observations to dedicated instrument response functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

So far the same instrument response function has been used for all
observations in our analysis examples.
You may have recognised that we always specified the calibration database
``prod2`` and the instrument response function ``South_0.5h`` that
are shipped with ctools.
When combining multiple observations in a joint analysis you may however
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
  Instrument response function [South_0.5h] North_0.5h
  Input model XML file [$CTOOLS/share/models/crab.xml] 
  Output event data file or observation definition XML file [events.fits] north.fits

.. code-block:: bash

  $ ctobssim
  RA of pointing (degrees) (0-360) [83.63] 
  Dec of pointing (degrees) (-90-90) [21.51] 22.51
  Radius of FOV (degrees) (0-180) [5.0] 
  Start time (MET in s) [0.0] 
  End time (MET in s) [1800.0] 
  Lower energy limit (TeV) [0.1] 
  Upper energy limit (TeV) [100.0] 
  Calibration database [prod2] 
  Instrument response function [North_0.5h] South_0.5h
  Input model XML file [$CTOOLS/share/models/crab.xml] 
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

You now jointly fit both observations using :ref:`ctlike`:

.. code-block:: bash

   $ ctlike chatter=3
  Input event list, counts cube or observation definition XML file [events.fits] obs_irf.xml
  Input model XML file [$CTOOLS/share/models/crab.xml] 
  Output model XML file [crab_results.xml] 

To see the usage of the individual response functions you may inspect the
log file (you need to set ``chatter=3`` to see the details of the 
observations in the log file).
You will notice that each observation now has a specific response
function and that the filename of the response information differs for
both observations.

.. code-block:: xml

  2015-12-07T22:22:18: +==============+
  2015-12-07T22:22:18: | Observations |
  2015-12-07T22:22:18: +==============+
  2015-12-07T22:22:18: === GObservations ===
  2015-12-07T22:22:18:  Number of observations ....: 2
  2015-12-07T22:22:18:  Number of models ..........: 2
  2015-12-07T22:22:18:  Number of predicted events : 0
  2015-12-07T22:22:18: === GCTAObservation ===
  2015-12-07T22:22:18:  Name ......................: Crab
  2015-12-07T22:22:18:  Identifier ................: 00001
  2015-12-07T22:22:18:  Instrument ................: CTA
  2015-12-07T22:22:18:  Event file ................: north.fits
  2015-12-07T22:22:18:  Event type ................: EventList
  2015-12-07T22:22:18:  Statistics ................: Poisson
  2015-12-07T22:22:18:  Ontime ....................: 1800 s
  2015-12-07T22:22:18:  Livetime ..................: 1710 s
  2015-12-07T22:22:18:  Deadtime correction .......: 0.95
  2015-12-07T22:22:18:  User energy range .........: undefined
  2015-12-07T22:22:18: === GCTAPointing ===
  2015-12-07T22:22:18:  Pointing direction ........: (RA,Dec)=(83.63,21.51)
  2015-12-07T22:22:18: === GCTAResponseIrf ===
  2015-12-07T22:22:18:  Caldb mission .............: cta
  2015-12-07T22:22:18:  Caldb instrument ..........: prod2
  2015-12-07T22:22:18:  Response name .............: North_0.5h
  2015-12-07T22:22:18:  Energy dispersion .........: Not used
  2015-12-07T22:22:18:  Save energy range .........: undefined
  2015-12-07T22:22:18: === GCTAEventList ===
  2015-12-07T22:22:18:  Number of events ..........: 8469
  2015-12-07T22:22:18:  Time interval .............: 51544.5 - 51544.5 days
  2015-12-07T22:22:18:  Energy interval ...........: 0.1 - 100 TeV
  2015-12-07T22:22:18:  Region of interest ........: RA=83.63, DEC=21.51 [0,0] Radius=5 deg
  2015-12-07T22:22:18: === GCTAObservation ===
  2015-12-07T22:22:18:  Name ......................: Crab
  2015-12-07T22:22:18:  Identifier ................: 00002
  2015-12-07T22:22:18:  Instrument ................: CTA
  2015-12-07T22:22:18:  Event file ................: south.fits
  2015-12-07T22:22:18:  Event type ................: EventList
  2015-12-07T22:22:18:  Statistics ................: Poisson
  2015-12-07T22:22:18:  Ontime ....................: 1800 s
  2015-12-07T22:22:18:  Livetime ..................: 1710 s
  2015-12-07T22:22:18:  Deadtime correction .......: 0.95
  2015-12-07T22:22:18:  User energy range .........: undefined
  2015-12-07T22:22:18: === GCTAPointing ===
  2015-12-07T22:22:18:  Pointing direction ........: (RA,Dec)=(83.63,22.51)
  2015-12-07T22:22:18: === GCTAResponseIrf ===
  2015-12-07T22:22:18:  Caldb mission .............: cta
  2015-12-07T22:22:18:  Caldb instrument ..........: prod2
  2015-12-07T22:22:18:  Response name .............: South_0.5h
  2015-12-07T22:22:18:  Energy dispersion .........: Not used
  2015-12-07T22:22:18:  Save energy range .........: undefined
  2015-12-07T22:22:18: === GCTAEventList ===
  2015-12-07T22:22:18:  Number of events ..........: 22529
  2015-12-07T22:22:18:  Time interval .............: 51544.5 - 51544.5 days
  2015-12-07T22:22:18:  Energy interval ...........: 0.1 - 100 TeV
  2015-12-07T22:22:18:  Region of interest ........: RA=83.63, DEC=22.51 [0,0] Radius=5 deg
  ...
  2015-12-07T22:22:18: +=================================+
  2015-12-07T22:22:18: | Maximum likelihood optimisation |
  2015-12-07T22:22:18: +=================================+
  2015-12-07T22:22:18:  >Iteration   0: -logL=214691.468, Lambda=1.0e-03
  2015-12-07T22:22:18:  >Iteration   1: -logL=214679.380, Lambda=1.0e-03, delta=12.088, max(|grad|)=9.713100 [Index:3]
  2015-12-07T22:22:18:  >Iteration   2: -logL=214679.365, Lambda=1.0e-04, delta=0.015, max(|grad|)=0.092975 [Index:3]
  2015-12-07T22:22:18:  >Iteration   3: -logL=214679.365, Lambda=1.0e-05, delta=0.000, max(|grad|)=0.000749 [Index:3]
  ...
  2015-12-07T22:22:18: +=========================================+
  2015-12-07T22:22:18: | Maximum likelihood optimisation results |
  2015-12-07T22:22:18: +=========================================+
  2015-12-07T22:22:18: === GOptimizerLM ===
  2015-12-07T22:22:18:  Optimized function value ..: 214679.365
  2015-12-07T22:22:18:  Absolute precision ........: 0.005
  2015-12-07T22:22:18:  Acceptable value decrease .: 2
  2015-12-07T22:22:18:  Optimization status .......: converged
  2015-12-07T22:22:18:  Number of parameters ......: 10
  2015-12-07T22:22:18:  Number of free parameters .: 4
  2015-12-07T22:22:18:  Number of iterations ......: 3
  2015-12-07T22:22:18:  Lambda ....................: 1e-06
  2015-12-07T22:22:18:  Maximum log likelihood ....: -214679.365
  2015-12-07T22:22:18:  Observed events  (Nobs) ...: 30998.000
  2015-12-07T22:22:18:  Predicted events (Npred) ..: 30998.000 (Nobs - Npred = 1.02755e-06)
  2015-12-07T22:22:18: === GModels ===
  2015-12-07T22:22:18:  Number of models ..........: 2
  2015-12-07T22:22:18:  Number of parameters ......: 10
  2015-12-07T22:22:18: === GModelSky ===
  2015-12-07T22:22:18:  Name ......................: Crab
  2015-12-07T22:22:18:  Instruments ...............: all
  2015-12-07T22:22:18:  Instrument scale factors ..: unity
  2015-12-07T22:22:18:  Observation identifiers ...: all
  2015-12-07T22:22:18:  Model type ................: PointSource
  2015-12-07T22:22:18:  Model components ..........: "SkyDirFunction" * "PowerLaw" * "Constant"
  2015-12-07T22:22:18:  Number of parameters ......: 6
  2015-12-07T22:22:18:  Number of spatial par's ...: 2
  2015-12-07T22:22:18:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
  2015-12-07T22:22:18:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
  2015-12-07T22:22:18:  Number of spectral par's ..: 3
  2015-12-07T22:22:18:   Prefactor ................: 5.81237e-16 +/- 8.63534e-18 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
  2015-12-07T22:22:18:   Index ....................: -2.53954 +/- 0.0137432 [-0,-5]  (free,scale=-1,gradient)
  2015-12-07T22:22:18:   PivotEnergy ..............: 300000 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2015-12-07T22:22:18:  Number of temporal par's ..: 1
  2015-12-07T22:22:18:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
  2015-12-07T22:22:18: === GCTAModelIrfBackground ===
  2015-12-07T22:22:18:  Name ......................: CTABackgroundModel
  2015-12-07T22:22:18:  Instruments ...............: CTA
  2015-12-07T22:22:18:  Instrument scale factors ..: unity
  2015-12-07T22:22:18:  Observation identifiers ...: all
  2015-12-07T22:22:18:  Model type ................: "PowerLaw" * "Constant"
  2015-12-07T22:22:18:  Number of parameters ......: 4
  2015-12-07T22:22:18:  Number of spectral par's ..: 3
  2015-12-07T22:22:18:   Prefactor ................: 1.01067 +/- 0.0103997 [0.001,1000] ph/cm2/s/MeV (free,scale=1,gradient)
  2015-12-07T22:22:18:   Index ....................: 0.0135414 +/- 0.0064452 [-5,5]  (free,scale=1,gradient)
  2015-12-07T22:22:18:   PivotEnergy ..............: 1e+06 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2015-12-07T22:22:18:  Number of temporal par's ..: 1
  2015-12-07T22:22:18:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

A more fine grained control over the response function can be achieved by 
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
      <parameter name="BkgCube"      file="bkgcube1.fits"/>
    </observation>
    <observation name="Crab" id="00002" instrument="CTA">
      <parameter name="CountsCube"   file="cntcube2.fits"/>
      <parameter name="ExposureCube" file="expcube2.fits"/>
      <parameter name="PsfCube"      file="psfcube2.fits"/>
      <parameter name="BkgCube"      file="bkgcube2.fits"/>
    </observation>
  </observation_list>
