.. _sec_combining_obs:

Combining observations
~~~~~~~~~~~~~~~~~~~~~~

Generally, the CTA data you may want to analyse will not only be composed of
a single observation (a.k.a. run) but of a list of observations that should
be combined in a joint analysis.
ctools has the capability to collect individual observations in a list and
to perform for example a joint maximum likelihood fit of all observations 
in a single shot.
Here is an example that illustrates how to do that.

Let's start with the simulation of two 30 min long observations of the Crab
nebula, each offset by 0.5° from the nebula in opposite directions:

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
  Instrument response function [South_0.5h] 
  Input model XML file [$CTOOLS/share/models/crab.xml] 
  Output event data file or observation definition XML file [events.fits] events1.fits

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
  Instrument response function [South_0.5h] 
  Input model XML file [$CTOOLS/share/models/crab.xml] 
  Output event data file or observation definition XML file [events1.fits] events2.fits

This will produce the two event files ``events1.fits`` and ``events2.fits``
on disk.

As next step you have to create an **observation definition XML file** that
collects both observations in a list.
Here is how that file looks like:

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <observation_list title="observation library">
    <observation name="Crab" id="00001" instrument="CTA">
      <parameter name="EventList" file="events1.fits"/>
    </observation>
    <observation name="Crab" id="00002" instrument="CTA">
      <parameter name="EventList" file="events2.fits"/>
    </observation>
  </observation_list>

The file contains a single ``<observation_list>`` tag that contains
two ``<observation>`` tags that each define an observation.
Each observation has a ``name``, an ``id`` and an ``instrument``
attribute.
The ``name`` attribute can have any arbitrary value, and may be the
same for all observations.
**However, the** ``id`` **attribute needs to be a unique character
string for any given instrument.**
The ``instrument`` attribute is a case-sensitive string that identifies
the instrument with which the observation was taken.
**Please make sure that the instrument string is set correctly so that
ctools knows which instrument specific functions need to be called.**

.. note::

   The instrument string for a CTA observation is obviously ``CTA``.
   In case that you want to analyse data from an existing Imaging Air
   Cherenkov Telescope you can also set the instrument string to ``HESS``,
   ``MAGIC``, or ``VERITAS``. You may also combine observations from different
   telescopes for a joint analysis in an observation definition file.
   **Please recall that instrument strings are case sensitive.**

Now you are ready to do a joint maximum likelihood analysis using
:ref:`ctlike`. Type:

.. code-block:: bash

  $ ctlike chatter=3
  Input event list, counts cube or observation definition XML file [events.fits] obs.xml
  Calibration database [prod2] 
  Instrument response function [South_0.5h] 
  Input model XML file [$CTOOLS/share/models/crab.xml] 
  Output model XML file [crab_results.xml] 

Instead of providing an event list or a counts cube, you now provide the
filename of the observation definition XML file (here ``obs.xml``) as input
parameter. :ref:`ctlike` recognises this format and automatically performs a
joint maximum likelihood analysis. To see the result of this you needed to
specify the ``chatter=3`` parameter so that the chattiness is increased in
the log file:

.. code-block:: none

  2016-06-29T19:14:38: +==============+
  2016-06-29T19:14:38: | Observations |
  2016-06-29T19:14:38: +==============+
  2016-06-29T19:14:38: === GObservations ===
  2016-06-29T19:14:38:  Number of observations ....: 2
  2016-06-29T19:14:38:  Number of models ..........: 2
  2016-06-29T19:14:38:  Number of observed events .: 46028
  2016-06-29T19:14:38:  Number of predicted events : 0
  2016-06-29T19:14:38: === GCTAObservation ===
  2016-06-29T19:14:38:  Name ......................: Crab
  2016-06-29T19:14:38:  Identifier ................: 00001
  2016-06-29T19:14:38:  Instrument ................: CTA
  2016-06-29T19:14:38:  Event file ................: events1.fits
  2016-06-29T19:14:38:  Event type ................: EventList
  2016-06-29T19:14:38:  Statistics ................: Poisson
  2016-06-29T19:14:38:  Ontime ....................: 1800 s
  2016-06-29T19:14:38:  Livetime ..................: 1710 s
  2016-06-29T19:14:38:  Deadtime correction .......: 0.95
  2016-06-29T19:14:38:  User energy range .........: undefined
  2016-06-29T19:14:38: === GCTAPointing ===
  2016-06-29T19:14:38:  Pointing direction ........: (RA,Dec)=(83.63,21.51)
  2016-06-29T19:14:38: === GCTAResponseIrf ===
  2016-06-29T19:14:38:  Caldb mission .............: cta
  2016-06-29T19:14:38:  Caldb instrument ..........: prod2
  2016-06-29T19:14:38:  Response name .............: South_0.5h
  2016-06-29T19:14:38:  Energy dispersion .........: Not used
  2016-06-29T19:14:38:  Save energy range .........: undefined
  2016-06-29T19:14:38: === GCTAEventList ===
  2016-06-29T19:14:38:  Number of events ..........: 23014 (disposed in "events1.fits")
  2016-06-29T19:14:38:  Time interval .............: 51544.5 - 51544.5 days
  2016-06-29T19:14:38:  Energy interval ...........: 0.1 - 100 TeV
  2016-06-29T19:14:38:  Region of interest ........: RA=83.63, DEC=21.51 [0,0] Radius=5 deg
  2016-06-29T19:14:38: === GCTAObservation ===
  2016-06-29T19:14:38:  Name ......................: Crab
  2016-06-29T19:14:38:  Identifier ................: 00002
  2016-06-29T19:14:38:  Instrument ................: CTA
  2016-06-29T19:14:38:  Event file ................: events2.fits
  2016-06-29T19:14:38:  Event type ................: EventList
  2016-06-29T19:14:38:  Statistics ................: Poisson
  2016-06-29T19:14:38:  Ontime ....................: 1800 s
  2016-06-29T19:14:38:  Livetime ..................: 1710 s
  2016-06-29T19:14:38:  Deadtime correction .......: 0.95
  2016-06-29T19:14:38:  User energy range .........: undefined
  2016-06-29T19:14:38: === GCTAPointing ===
  2016-06-29T19:14:38:  Pointing direction ........: (RA,Dec)=(83.63,21.51)
  2016-06-29T19:14:38: === GCTAResponseIrf ===
  2016-06-29T19:14:38:  Caldb mission .............: cta
  2016-06-29T19:14:38:  Caldb instrument ..........: prod2
  2016-06-29T19:14:38:  Response name .............: South_0.5h
  2016-06-29T19:14:38:  Energy dispersion .........: Not used
  2016-06-29T19:14:38:  Save energy range .........: undefined
  2016-06-29T19:14:38: === GCTAEventList ===
  2016-06-29T19:14:38:  Number of events ..........: 23014 (disposed in "events2.fits")
  2016-06-29T19:14:38:  Time interval .............: 51544.5 - 51544.5 days
  2016-06-29T19:14:38:  Energy interval ...........: 0.1 - 100 TeV
  2016-06-29T19:14:38:  Region of interest ........: RA=83.63, DEC=21.51 [0,0] Radius=5 deg
  2016-06-29T19:14:38:
  2016-06-29T19:14:38: +=================================+
  2016-06-29T19:14:38: | Maximum likelihood optimisation |
  2016-06-29T19:14:38: +=================================+
  2016-06-29T19:14:38:  >Iteration   0: -logL=309355.758, Lambda=1.0e-03
  2016-06-29T19:14:38:  >Iteration   1: -logL=309352.158, Lambda=1.0e-03, delta=3.600, max(|grad|)=7.976081 [Index:7]
  2016-06-29T19:14:38:  >Iteration   2: -logL=309352.157, Lambda=1.0e-04, delta=0.001, max(|grad|)=0.033525 [Index:3]
  ...
  2016-06-29T19:14:38: +=========================================+
  2016-06-29T19:14:38: | Maximum likelihood optimisation results |
  2016-06-29T19:14:38: +=========================================+
  2016-06-29T19:14:38: === GOptimizerLM ===
  2016-06-29T19:14:38:  Optimized function value ..: 309352.157
  2016-06-29T19:14:38:  Absolute precision ........: 0.005
  2016-06-29T19:14:38:  Acceptable value decrease .: 2
  2016-06-29T19:14:38:  Optimization status .......: converged
  2016-06-29T19:14:38:  Number of parameters ......: 10
  2016-06-29T19:14:38:  Number of free parameters .: 4
  2016-06-29T19:14:38:  Number of iterations ......: 2
  2016-06-29T19:14:38:  Lambda ....................: 1e-05
  2016-06-29T19:14:38:  Maximum log likelihood ....: -309352.157
  2016-06-29T19:14:38:  Observed events  (Nobs) ...: 46028.000
  2016-06-29T19:14:38:  Predicted events (Npred) ..: 46027.996 (Nobs - Npred = 0.00423575)
  2016-06-29T19:14:38: === GModels ===
  2016-06-29T19:14:38:  Number of models ..........: 2
  2016-06-29T19:14:38:  Number of parameters ......: 10
  2016-06-29T19:14:38: === GModelSky ===
  2016-06-29T19:14:38:  Name ......................: Crab
  2016-06-29T19:14:38:  Instruments ...............: all
  2016-06-29T19:14:38:  Instrument scale factors ..: unity
  2016-06-29T19:14:38:  Observation identifiers ...: all
  2016-06-29T19:14:38:  Model type ................: PointSource
  2016-06-29T19:14:38:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
  2016-06-29T19:14:38:  Number of parameters ......: 6
  2016-06-29T19:14:38:  Number of spatial par's ...: 2
  2016-06-29T19:14:38:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
  2016-06-29T19:14:38:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
  2016-06-29T19:14:38:  Number of spectral par's ..: 3
  2016-06-29T19:14:38:   Prefactor ................: 5.7909e-16 +/- 7.3315e-18 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
  2016-06-29T19:14:38:   Index ....................: -2.47181 +/- 0.0109973 [-0,-5]  (free,scale=-1,gradient)
  2016-06-29T19:14:38:   PivotEnergy ..............: 300000 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2016-06-29T19:14:38:  Number of temporal par's ..: 1
  2016-06-29T19:14:38:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
  2016-06-29T19:14:38: === GCTAModelIrfBackground ===
  2016-06-29T19:14:38:  Name ......................: CTABackgroundModel
  2016-06-29T19:14:38:  Instruments ...............: CTA
  2016-06-29T19:14:38:  Instrument scale factors ..: unity
  2016-06-29T19:14:38:  Observation identifiers ...: all
  2016-06-29T19:14:38:  Model type ................: "PowerLaw" * "Constant"
  2016-06-29T19:14:38:  Number of parameters ......: 4
  2016-06-29T19:14:38:  Number of spectral par's ..: 3
  2016-06-29T19:14:38:   Prefactor ................: 1.01608 +/- 0.0084739 [0.001,1000] ph/cm2/s/MeV (free,scale=1,gradient)
  2016-06-29T19:14:38:   Index ....................: 0.00524546 +/- 0.00516392 [-5,5]  (free,scale=1,gradient)
  2016-06-29T19:14:38:   PivotEnergy ..............: 1e+06 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2016-06-29T19:14:38:  Number of temporal par's ..: 1
  2016-06-29T19:14:38:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

The log file indicates that the fit converged quickly, the spectral
parameters of the Crab nebula have now been constrained using the events
from both observations.
The computation time increases roughly linearly with the number of
observations that are combined, although ctools implements parallel 
multi-core processing which will spread the likelihood computation for 
the different observations over all CPU cores that are available. 
**Doing a joint unbinned analysis is thus an efficient solution if
data from multiple observations should be combined.**

Combining observations is not limited to unbinned data (i.e. event lists)
but may also be applied to binned data (i.e. counts cubes).
Using :ref:`ctbin` you can create counts cubes from both event lists which
may then be combined in an observation definition XML file:

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <observation_list title="observation library">
    <observation name="Crab" id="00001" instrument="CTA">
      <parameter name="CountsCube" file="cntcube1.fits"/>
    </observation>
    <observation name="Crab" id="00002" instrument="CTA">
      <parameter name="CountsCube" file="cntcube2.fits"/>
    </observation>
  </observation_list>

Feeding the observation definition XML file to :ref:`ctlike` will then
lead to a joint binned analysis.
In the joint binned analysis, the events of individual observations are
not combined, but are kept separate in distinct counts cubes.
This is not very efficient, as generally counts cubes for short duration
observations are only sparsly populated and the likelihood computation 
has to loop over a hugh number of data space bins (though also here
:ref:`ctlike` benefits from multi-core parallel processing).
**Though possible, a joint binned analysis is thus not the recommended
method for combining observations.**
An alternative is to stack the events of all observations into a single
counts cube.
The :ref:`following section <sec_stacked>` describes how such a stacked
analysis is done with ctools.

.. note::

  Given that logic, unbinned and binned observations may also be combined
  in a joint analysis, although this Use Case may be a bit academic:

  .. code-block:: xml

    <?xml version="1.0" standalone="no"?>
    <observation_list title="observation library">
      <observation name="Crab" id="00001" instrument="CTA">
        <parameter name="EventList" file="events1.fits"/>
      </observation>
      <observation name="Crab" id="00002" instrument="CTA">
        <parameter name="CountsCube" file="cntcube2.fits"/>
      </observation>
    </observation_list>
