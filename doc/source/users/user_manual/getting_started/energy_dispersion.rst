.. _sec_edisp_cta:

Taking into account energy dispersion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may not have recognised it, but the examples you have exercised so far
have neglected the impact of energy dispersion on the analysis. In reality,
however, the reconstructed event energy will differ from the true photon energy,
and this effect will become particularily important at low energies. There are
therefore cases where you want to take the effects of energy dispersion into
account, and you can do this in all ctools that make use of the instrument
response functions by setting the hidden ``edisp`` parameter to ``yes``.

.. warning::

   The maximum likelihood computations including energy dispersion are
   relatively time consuming, and in many situations the impact of the
   energy dispersion on the analysis results will be very small. For this
   reason, energy dispersion is by default disabled in ctools.

To simulate events taking into account the energy dispersion you should
run the :ref:`ctobssim` tool with ``edisp=yes`` parameter:

.. code-block:: bash

  $ ctobssim edisp=yes
  RA of pointing (degrees) (0-360) [83.63] 
  Dec of pointing (degrees) (-90-90) [22.01] 
  Radius of FOV (degrees) (0-180) [5.0] 
  Start time (MET in s) [0.0] 
  End time (MET in s) [1800.0] 
  Lower energy limit (TeV) [0.1] 
  Upper energy limit (TeV) [100.0] 
  Calibration database [prod2] 
  Instrument response function [South_0.5h] 
  Input model XML file [$CTOOLS/share/models/crab.xml] 
  Output event data file or observation definition XML file [events.fits]

You may then select events as before using the :ref:`ctselect` tool
(since this tool does not need the instrument response function it has no
hidden ``edisp`` parameter):

.. code-block:: bash

  $ ctselect
  Input event list or observation definition XML file [events.fits]
  RA for ROI centre (degrees) (0-360) [83.63]
  Dec for ROI centre (degrees) (-90-90) [22.01]
  Radius of ROI (degrees) (0-180) [3.0]
  Start time (CTA MET in seconds) [0.0]
  End time (CTA MET in seconds) [0.0]
  Lower energy limit (TeV) [0.1]
  Upper energy limit (TeV) [100.0]
  Output event list or observation definition XML file [selected_events.fits]

You then perform a maximum likelihood fit using the :ref:`ctlike` tool
with ``edisp=yes`` so that energy dispersion is correctly taken into account
in the likelihood computation:

.. code-block:: bash

  $ ctlike edisp=yes
  Input event list, counts cube or observation definition XML file [selected_events.fits]
  Calibration database [prod2]
  Instrument response function [South_0.5h]
  Input model XML file [$CTOOLS/share/models/crab.xml]
  Output model XML file [crab_results.xml]

Below you see the output from the ``ctlike.log`` file.

.. code-block:: none

  2016-06-29T20:33:36: +=================================+
  2016-06-29T20:33:36: | Maximum likelihood optimisation |
  2016-06-29T20:33:36: +=================================+
  2016-06-29T20:33:39:  >Iteration   0: -logL=133777.933, Lambda=1.0e-03
  2016-06-29T20:33:42:  >Iteration   1: -logL=133772.125, Lambda=1.0e-03, delta=5.808, max(|grad|)=23.221525 [Index:7]
  2016-06-29T20:33:45:  >Iteration   2: -logL=133772.111, Lambda=1.0e-04, delta=0.013, max(|grad|)=-0.042654 [Prefactor:6]
  2016-06-29T20:33:48:  >Iteration   3: -logL=133772.111, Lambda=1.0e-05, delta=0.000, max(|grad|)=-0.000611 [Index:3]
  ...
  2016-06-29T20:33:50: +=========================================+
  2016-06-29T20:33:50: | Maximum likelihood optimisation results |
  2016-06-29T20:33:50: +=========================================+
  2016-06-29T20:33:50: === GOptimizerLM ===
  2016-06-29T20:33:50:  Optimized function value ..: 133772.111
  2016-06-29T20:33:50:  Absolute precision ........: 0.005
  2016-06-29T20:33:50:  Acceptable value decrease .: 2
  2016-06-29T20:33:50:  Optimization status .......: converged
  2016-06-29T20:33:50:  Number of parameters ......: 10
  2016-06-29T20:33:50:  Number of free parameters .: 4
  2016-06-29T20:33:50:  Number of iterations ......: 3
  2016-06-29T20:33:50:  Lambda ....................: 1e-06
  2016-06-29T20:33:50:  Maximum log likelihood ....: -133772.111
  2016-06-29T20:33:50:  Observed events  (Nobs) ...: 20952.000
  2016-06-29T20:33:50:  Predicted events (Npred) ..: 20952.000 (Nobs - Npred = 1.01277e-06)
  2016-06-29T20:33:50: === GModels ===
  2016-06-29T20:33:50:  Number of models ..........: 2
  2016-06-29T20:33:50:  Number of parameters ......: 10
  2016-06-29T20:33:50: === GModelSky ===
  2016-06-29T20:33:50:  Name ......................: Crab
  2016-06-29T20:33:50:  Instruments ...............: all
  2016-06-29T20:33:50:  Instrument scale factors ..: unity
  2016-06-29T20:33:50:  Observation identifiers ...: all
  2016-06-29T20:33:50:  Model type ................: PointSource
  2016-06-29T20:33:50:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
  2016-06-29T20:33:50:  Number of parameters ......: 6
  2016-06-29T20:33:50:  Number of spatial par's ...: 2
  2016-06-29T20:33:50:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
  2016-06-29T20:33:50:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
  2016-06-29T20:33:50:  Number of spectral par's ..: 3
  2016-06-29T20:33:50:   Prefactor ................: 5.74029e-16 +/- 1.04025e-17 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
  2016-06-29T20:33:50:   Index ....................: -2.49899 +/- 0.015621 [-0,-5]  (free,scale=-1,gradient)
  2016-06-29T20:33:50:   PivotEnergy ..............: 300000 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2016-06-29T20:33:50:  Number of temporal par's ..: 1
  2016-06-29T20:33:50:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
  2016-06-29T20:33:50: === GCTAModelIrfBackground ===
  2016-06-29T20:33:50:  Name ......................: CTABackgroundModel
  2016-06-29T20:33:50:  Instruments ...............: CTA
  2016-06-29T20:33:50:  Instrument scale factors ..: unity
  2016-06-29T20:33:50:  Observation identifiers ...: all
  2016-06-29T20:33:50:  Model type ................: "PowerLaw" * "Constant"
  2016-06-29T20:33:50:  Number of parameters ......: 4
  2016-06-29T20:33:50:  Number of spectral par's ..: 3
  2016-06-29T20:33:50:   Prefactor ................: 0.959074 +/- 0.0134496 [0.001,1000] ph/cm2/s/MeV (free,scale=1,gradient)
  2016-06-29T20:33:50:   Index ....................: -0.0162199 +/- 0.00850751 [-5,5]  (free,scale=1,gradient)
  2016-06-29T20:33:50:   PivotEnergy ..............: 1e+06 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2016-06-29T20:33:50:  Number of temporal par's ..: 1
  2016-06-29T20:33:50:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

