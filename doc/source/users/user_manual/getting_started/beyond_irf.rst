.. _sec_connecting_irf:

Connecting observations to dedicated instrument response functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  2016-06-29T20:09:44: +============+
  2016-06-29T20:09:44: | Parameters |
  2016-06-29T20:09:44: +============+
  2016-06-29T20:09:44:  chatter ...................: 2
  2016-06-29T20:09:44:  clobber ...................: yes
  2016-06-29T20:09:44:  debug .....................: yes
  2016-06-29T20:09:44:  mode ......................: ql
  2016-06-29T20:09:44:  logfile ...................: cscaldb.log
  2016-06-29T20:09:44:
  2016-06-29T20:09:44: +==============+
  2016-06-29T20:09:44: | Mission: cta |
  2016-06-29T20:09:44: +==============+
  2016-06-29T20:09:44: === Response functions in database "prod2" ===
  2016-06-29T20:09:44: North_0.5h
  2016-06-29T20:09:44: North_50h
  2016-06-29T20:09:44: North_5h
  2016-06-29T20:09:44: South_0.5h
  2016-06-29T20:09:44: South_50h
  2016-06-29T20:09:44: South_5h

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

.. code-block:: none

  2016-06-29T20:10:29: +==============+
  2016-06-29T20:10:29: | Observations |
  2016-06-29T20:10:29: +==============+
  2016-06-29T20:10:29: === GObservations ===
  2016-06-29T20:10:29:  Number of observations ....: 2
  2016-06-29T20:10:29:  Number of models ..........: 2
  2016-06-29T20:10:29:  Number of observed events .: 31270
  2016-06-29T20:10:29:  Number of predicted events : 0
  2016-06-29T20:10:29: === GCTAObservation ===
  2016-06-29T20:10:29:  Name ......................: Crab
  2016-06-29T20:10:29:  Identifier ................: 00001
  ...
  2016-06-29T20:10:29: === GCTAResponseIrf ===
  2016-06-29T20:10:29:  Caldb mission .............: cta
  2016-06-29T20:10:29:  Caldb instrument ..........: prod2
  2016-06-29T20:10:29:  Response name .............: North_0.5h
  2016-06-29T20:10:29:  Energy dispersion .........: Not used
  2016-06-29T20:10:29:  Save energy range .........: undefined
  ...
  2016-06-29T20:10:29: === GCTAObservation ===
  2016-06-29T20:10:29:  Name ......................: Crab
  2016-06-29T20:10:29:  Identifier ................: 00002
  ...
  2016-06-29T20:10:29: === GCTAResponseIrf ===
  2016-06-29T20:10:29:  Caldb mission .............: cta
  2016-06-29T20:10:29:  Caldb instrument ..........: prod2
  2016-06-29T20:10:29:  Response name .............: South_0.5h
  2016-06-29T20:10:29:  Energy dispersion .........: Not used
  2016-06-29T20:10:29:  Save energy range .........: undefined
  ...
  2016-06-29T20:10:29: +=================================+
  2016-06-29T20:10:29: | Maximum likelihood optimisation |
  2016-06-29T20:10:29: +=================================+
  2016-06-29T20:10:30:  >Iteration   0: -logL=216208.923, Lambda=1.0e-03
  2016-06-29T20:10:30:  >Iteration   1: -logL=216207.021, Lambda=1.0e-03, delta=1.902, max(|grad|)=4.436064 [Index:7]
  2016-06-29T20:10:30:  >Iteration   2: -logL=216207.021, Lambda=1.0e-04, delta=0.001, max(|grad|)=0.012146 [Index:7]
  ...
  2016-06-29T20:10:30: +=========================================+
  2016-06-29T20:10:30: | Maximum likelihood optimisation results |
  2016-06-29T20:10:30: +=========================================+
  2016-06-29T20:10:30: === GOptimizerLM ===
  2016-06-29T20:10:30:  Optimized function value ..: 216207.021
  2016-06-29T20:10:30:  Absolute precision ........: 0.005
  2016-06-29T20:10:30:  Acceptable value decrease .: 2
  2016-06-29T20:10:30:  Optimization status .......: converged
  2016-06-29T20:10:30:  Number of parameters ......: 10
  2016-06-29T20:10:30:  Number of free parameters .: 4
  2016-06-29T20:10:30:  Number of iterations ......: 2
  2016-06-29T20:10:30:  Lambda ....................: 1e-05
  2016-06-29T20:10:30:  Maximum log likelihood ....: -216207.021
  2016-06-29T20:10:30:  Observed events  (Nobs) ...: 31270.000
  2016-06-29T20:10:30:  Predicted events (Npred) ..: 31269.998 (Nobs - Npred = 0.00190254)
  2016-06-29T20:10:30: === GModels ===
  2016-06-29T20:10:30:  Number of models ..........: 2
  2016-06-29T20:10:30:  Number of parameters ......: 10
  2016-06-29T20:10:30: === GModelSky ===
  2016-06-29T20:10:30:  Name ......................: Crab
  2016-06-29T20:10:30:  Instruments ...............: all
  2016-06-29T20:10:30:  Instrument scale factors ..: unity
  2016-06-29T20:10:30:  Observation identifiers ...: all
  2016-06-29T20:10:30:  Model type ................: PointSource
  2016-06-29T20:10:30:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
  2016-06-29T20:10:30:  Number of parameters ......: 6
  2016-06-29T20:10:30:  Number of spatial par's ...: 2
  2016-06-29T20:10:30:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
  2016-06-29T20:10:30:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
  2016-06-29T20:10:30:  Number of spectral par's ..: 3
  2016-06-29T20:10:30:   Prefactor ................: 5.7247e-16 +/- 8.64521e-18 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
  2016-06-29T20:10:30:   Index ....................: -2.46519 +/- 0.013273 [-0,-5]  (free,scale=-1,gradient)
  2016-06-29T20:10:30:   PivotEnergy ..............: 300000 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2016-06-29T20:10:30:  Number of temporal par's ..: 1
  2016-06-29T20:10:30:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
  2016-06-29T20:10:30: === GCTAModelIrfBackground ===
  2016-06-29T20:10:30:  Name ......................: CTABackgroundModel
  2016-06-29T20:10:30:  Instruments ...............: CTA
  2016-06-29T20:10:30:  Instrument scale factors ..: unity
  2016-06-29T20:10:30:  Observation identifiers ...: all
  2016-06-29T20:10:30:  Model type ................: "PowerLaw" * "Constant"
  2016-06-29T20:10:30:  Number of parameters ......: 4
  2016-06-29T20:10:30:  Number of spectral par's ..: 3
  2016-06-29T20:10:30:   Prefactor ................: 1.01522 +/- 0.0104442 [0.001,1000] ph/cm2/s/MeV (free,scale=1,gradient)
  2016-06-29T20:10:30:   Index ....................: 0.00802507 +/- 0.00643391 [-5,5]  (free,scale=1,gradient)
  2016-06-29T20:10:30:   PivotEnergy ..............: 1e+06 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2016-06-29T20:10:30:  Number of temporal par's ..: 1
  2016-06-29T20:10:30:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

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
