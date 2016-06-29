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
nebula, each offset by 0.5 deg from the nebula in opposite directions:

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
   In case that data from any existing Imaging Air Cherenkov Telescope
   are provided in the correct format, the instrument strings
   ``HESS``, ``MAGIC`` and ``VERITAS`` are also recognised. This opens
   the obvious possibility to combine data from different
   Imaging Air Cherenkov Telescope in a joint analysis.
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

Instead of providing an event list or a counts cube, you now provided the
filename of the observation definition XML file (here ``obs.xml``) as input
parameter. :ref:`ctlike` recognises this format and automatically performs a
joint analysis. To see the result of this you needed to specify the
``chatter=3`` parameter so that the chattiness is increased in the log 
file:

.. code-block:: xml

  2015-12-07T21:23:07: +==============+
  2015-12-07T21:23:07: | Observations |
  2015-12-07T21:23:07: +==============+
  2015-12-07T21:23:07: === GObservations ===
  2015-12-07T21:23:07:  Number of observations ....: 2
  2015-12-07T21:23:07:  Number of models ..........: 2
  2015-12-07T21:23:07:  Number of predicted events : 0
  2015-12-07T21:23:07: === GCTAObservation ===
  2015-12-07T21:23:07:  Name ......................: Crab
  2015-12-07T21:23:07:  Identifier ................: 00001
  2015-12-07T21:23:07:  Instrument ................: CTA
  2015-12-07T21:23:07:  Event file ................: events1.fits
  2015-12-07T21:23:07:  Event type ................: EventList
  2015-12-07T21:23:07:  Statistics ................: Poisson
  2015-12-07T21:23:07:  Ontime ....................: 1800 s
  2015-12-07T21:23:07:  Livetime ..................: 1710 s
  2015-12-07T21:23:07:  Deadtime correction .......: 0.95
  2015-12-07T21:23:07:  User energy range .........: undefined
  2015-12-07T21:23:07: === GCTAPointing ===
  2015-12-07T21:23:07:  Pointing direction ........: (RA,Dec)=(83.63,21.51)
  2015-12-07T21:23:07: === GCTAResponseIrf ===
  2015-12-07T21:23:07:  Caldb mission .............: cta
  2015-12-07T21:23:07:  Caldb instrument ..........: prod2
  2015-12-07T21:23:07:  Response name .............: South_0.5h
  2015-12-07T21:23:07:  Energy dispersion .........: Not used
  2015-12-07T21:23:07:  Save energy range .........: undefined
  2015-12-07T21:23:07: === GCTAEventList ===
  2015-12-07T21:23:07:  Number of events ..........: 22593
  2015-12-07T21:23:07:  Time interval .............: 51544.5 - 51544.5 days
  2015-12-07T21:23:07:  Energy interval ...........: 0.1 - 100 TeV
  2015-12-07T21:23:07:  Region of interest ........: RA=83.63, DEC=21.51 [0,0] Radius=5 deg
  2015-12-07T21:23:07: === GCTAObservation ===
  2015-12-07T21:23:07:  Name ......................: Crab
  2015-12-07T21:23:07:  Identifier ................: 00002
  2015-12-07T21:23:07:  Instrument ................: CTA
  2015-12-07T21:23:07:  Event file ................: events2.fits
  2015-12-07T21:23:07:  Event type ................: EventList
  2015-12-07T21:23:07:  Statistics ................: Poisson
  2015-12-07T21:23:07:  Ontime ....................: 1800 s
  2015-12-07T21:23:07:  Livetime ..................: 1710 s
  2015-12-07T21:23:07:  Deadtime correction .......: 0.95
  2015-12-07T21:23:07:  User energy range .........: undefined
  2015-12-07T21:23:07: === GCTAPointing ===
  2015-12-07T21:23:07:  Pointing direction ........: (RA,Dec)=(83.63,22.51)
  2015-12-07T21:23:07: === GCTAResponseIrf ===
  2015-12-07T21:23:07:  Caldb mission .............: cta
  2015-12-07T21:23:07:  Caldb instrument ..........: prod2
  2015-12-07T21:23:07:  Response name .............: South_0.5h
  2015-12-07T21:23:07:  Energy dispersion .........: Not used
  2015-12-07T21:23:07:  Save energy range .........: undefined
  2015-12-07T21:23:07: === GCTAEventList ===
  2015-12-07T21:23:07:  Number of events ..........: 22529
  2015-12-07T21:23:07:  Time interval .............: 51544.5 - 51544.5 days
  2015-12-07T21:23:07:  Energy interval ...........: 0.1 - 100 TeV
  2015-12-07T21:23:07:  Region of interest ........: RA=83.63, DEC=22.51 [0,0] Radius=5 deg
  ...
  2015-12-07T21:23:07: +=================================+
  2015-12-07T21:23:07: | Maximum likelihood optimisation |
  2015-12-07T21:23:07: +=================================+
  2015-12-07T21:23:07:  >Iteration   0: -logL=304176.886, Lambda=1.0e-03
  2015-12-07T21:23:08:  >Iteration   1: -logL=304161.799, Lambda=1.0e-03, delta=15.088, max(|grad|)=12.922804 [Index:3]
  2015-12-07T21:23:08:  >Iteration   2: -logL=304161.783, Lambda=1.0e-04, delta=0.015, max(|grad|)=0.098155 [Index:3]
  2015-12-07T21:23:08:  >Iteration   3: -logL=304161.783, Lambda=1.0e-05, delta=0.000, max(|grad|)=0.000598 [Index:3]
  ...

  2015-12-07T21:23:08: +=========================================+
  2015-12-07T21:23:08: | Maximum likelihood optimisation results |
  2015-12-07T21:23:08: +=========================================+
  2015-12-07T21:23:08: === GOptimizerLM ===
  2015-12-07T21:23:08:  Optimized function value ..: 304161.783
  2015-12-07T21:23:08:  Absolute precision ........: 0.005
  2015-12-07T21:23:08:  Acceptable value decrease .: 2
  2015-12-07T21:23:08:  Optimization status .......: converged
  2015-12-07T21:23:08:  Number of parameters ......: 10
  2015-12-07T21:23:08:  Number of free parameters .: 4
  2015-12-07T21:23:08:  Number of iterations ......: 3
  2015-12-07T21:23:08:  Lambda ....................: 1e-06
  2015-12-07T21:23:08:  Maximum log likelihood ....: -304161.783
  2015-12-07T21:23:08:  Observed events  (Nobs) ...: 45122.000
  2015-12-07T21:23:08:  Predicted events (Npred) ..: 45122.000 (Nobs - Npred = 1.04131e-06)
  2015-12-07T21:23:08: === GModels ===
  2015-12-07T21:23:08:  Number of models ..........: 2
  2015-12-07T21:23:08:  Number of parameters ......: 10
  2015-12-07T21:23:08: === GModelSky ===
  2015-12-07T21:23:08:  Name ......................: Crab
  2015-12-07T21:23:08:  Instruments ...............: all
  2015-12-07T21:23:08:  Instrument scale factors ..: unity
  2015-12-07T21:23:08:  Observation identifiers ...: all
  2015-12-07T21:23:08:  Model type ................: PointSource
  2015-12-07T21:23:08:  Model components ..........: "SkyDirFunction" * "PowerLaw" * "Constant"
  2015-12-07T21:23:08:  Number of parameters ......: 6
  2015-12-07T21:23:08:  Number of spatial par's ...: 2
  2015-12-07T21:23:08:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
  2015-12-07T21:23:08:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
  2015-12-07T21:23:08:  Number of spectral par's ..: 3
  2015-12-07T21:23:08:   Prefactor ................: 5.78382e-16 +/- 7.27218e-18 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
  2015-12-07T21:23:08:   Index ....................: -2.5314 +/- 0.0112546 [-0,-5]  (free,scale=-1,gradient)
  2015-12-07T21:23:08:   PivotEnergy ..............: 300000 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2015-12-07T21:23:08:  Number of temporal par's ..: 1
  2015-12-07T21:23:08:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
  2015-12-07T21:23:08: === GCTAModelIrfBackground ===
  2015-12-07T21:23:08:  Name ......................: CTABackgroundModel
  2015-12-07T21:23:08:  Instruments ...............: CTA
  2015-12-07T21:23:08:  Instrument scale factors ..: unity
  2015-12-07T21:23:08:  Observation identifiers ...: all
  2015-12-07T21:23:08:  Model type ................: "PowerLaw" * "Constant"
  2015-12-07T21:23:08:  Number of parameters ......: 4
  2015-12-07T21:23:08:  Number of spectral par's ..: 3
  2015-12-07T21:23:08:   Prefactor ................: 1.0008 +/- 0.00839192 [0.001,1000] ph/cm2/s/MeV (free,scale=1,gradient)
  2015-12-07T21:23:08:   Index ....................: 0.0102998 +/- 0.0051965 [-5,5]  (free,scale=1,gradient)
  2015-12-07T21:23:08:   PivotEnergy ..............: 1e+06 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2015-12-07T21:23:08:  Number of temporal par's ..: 1
  2015-12-07T21:23:08:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

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
Using :ref:`ctbin` we can create counts cubes from both event lists which
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
