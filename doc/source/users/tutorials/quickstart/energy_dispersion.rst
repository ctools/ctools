.. _start_edisp:

Taking the energy dispersion into account
-----------------------------------------

  .. admonition:: What you will learn

     You will learn how to **take the energy dispersion (or energy redistribution)
     into account**.

     Although the effect of the energy dispersion is often neglegible there
     may be cases where you want to consider energy dispersion in an analysis,
     for example if you are analysing the data down to very low energies.

  .. warning::
     Energy dispersion is fully implemented in ctools but is computationally
     intensive. So be aware that the **tools and scripts will take a substantial
     amount of computing time if energy dispersion is considered**.

You may not have recognised, but the examples you have exercised so far
have neglected the impact of the energy dispersion on the analysis. In reality,
however, the reconstructed event energy will differ from the true photon energy,
and this effect will become particularily important at low energies. There are
therefore cases where you want to take the energy dispersion into account.

To simulate events taking the energy dispersion into account you run the
:ref:`ctobssim` tool with ``edisp=yes`` parameter:

.. code-block:: bash

   $ ctobssim edisp=yes
   RA of pointing (degrees) (0-360) [83.63]
   Dec of pointing (degrees) (-90-90) [22.51]
   Radius of FOV (degrees) (0-180) [5.0]
   Start time (UTC string, JD, MJD or MET in seconds) [2020-01-01T00:00:00]
   Stop time (UTC string, JD, MJD or MET in seconds) [2020-01-01T01:00:00]
   Lower energy limit (TeV) [0.03]
   Upper energy limit (TeV) [150.0]
   Calibration database [prod2]
   Instrument response function [South_0.5h]
   Input model definition XML file [$CTOOLS/share/models/crab.xml]
   Output event data file or observation definition XML file [events.fits] events_edisp.fits

You then select events from the simulated data as before using the
:ref:`ctselect` tool:

.. code-block:: bash

   $ ctselect
   Input event list or observation definition XML file [events.fits] events_edisp.fits
   Radius of ROI around pointing or specified RA/DEC (degrees) (0-180) [3.0]
   Start time (UTC string, JD, MJD or MET in seconds) [2020-01-01T00:10:00]
   Stop time (UTC string, JD, MJD or MET in seconds) [2020-01-01T00:40:00]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Output event list or observation definition XML file [selected_events.fits] selected_events_edisp.fits

Then you bin as before the selected events into a counts cube using the
:ref:`ctbin` tool:

.. code-block:: bash

   $ ctbin
   Input event list or observation definition XML file [selected_events.fits] selected_events_edisp.fits
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL]
   Projection method (AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN) [CAR]
   First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63]
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.51]
   Image scale (in degrees/pixel) [0.02]
   Size of the X axis in pixels [200]
   Size of the Y axis in pixels [200]
   Algorithm for defining energy bins (FILE|LIN|LOG|POW) [LOG]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins (1-200) [20]
   Output counts cube file or observation definition XML file [cntcube.fits] cntcube_edisp.fits

As next step you need to compute the
:ref:`energy dispersion cube <glossary_edispcube>`
using the :ref:`ctexpcube` tool. You run the tool as follows:

.. code-block:: bash

   $ ctedispcube
   Input event list or observation definition XML file [NONE] selected_events_edisp.fits
   Calibration database [prod2]
   Instrument response function [South_0.5h]
   Input counts cube file to extract energy dispersion cube definition [NONE]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL]
   Projection method (AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN) [CAR]
   First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63]
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.51]
   Image scale (in degrees/pixel) [1.0]
   Size of the X axis in pixels [10]
   Size of the Y axis in pixels [10]
   Algorithm for defining energy bins (FILE|LIN|LOG|POW) [LOG]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins (1-1000) [20]
   Output energy dispersion cube file [edispcube.fits]

Now you are ready to perform a binned maximum likelihood analysis taking the
energy dispersion into account. You do this by running the :ref:`ctlike` tool
with the ``edisp=yes`` parameter. The :ref:`ctlike` tool will now query for the
energy dispersion cube:

.. code-block:: bash

   $ ctlike edisp=yes
   Input event list, counts cube or observation definition XML file [selected_events.fits] cntcube_edisp.fits
   Input exposure cube file [expcube.fits]
   Input PSF cube file [psfcube.fits]
   Input background cube file [bkgcube.fits]
   Input energy dispersion cube file [NONE] edispcube.fits
   Input model definition XML file [$CTOOLS/share/models/crab.xml] models.xml
   Output model definition XML file [crab_results.xml] crab_results_edisp.xml

And here is the output in the log file:

.. code-block:: none

   2019-04-02T14:20:38: +=================================+
   2019-04-02T14:20:38: | Maximum likelihood optimisation |
   2019-04-02T14:20:38: +=================================+
   2019-04-02T14:21:49:  >Iteration   0: -logL=57889.195, Lambda=1.0e-03
   2019-04-02T14:22:51:  >Iteration   1: -logL=57887.419, Lambda=1.0e-03, delta=1.776, step=1.0e+00, max(|grad|)=1.111374 [Index:3]
   2019-04-02T14:23:54:  >Iteration   2: -logL=57887.417, Lambda=1.0e-04, delta=0.002, step=1.0e+00, max(|grad|)=-0.005923 [Index:7]
   2019-04-02T14:24:57:
   2019-04-02T14:24:57: +=========================================+
   2019-04-02T14:24:57: | Maximum likelihood optimisation results |
   2019-04-02T14:24:57: +=========================================+
   2019-04-02T14:24:57: === GOptimizerLM ===
   2019-04-02T14:24:57:  Optimized function value ..: 57887.417
   2019-04-02T14:24:57:  Absolute precision ........: 0.005
   2019-04-02T14:24:57:  Acceptable value decrease .: 2
   2019-04-02T14:24:57:  Optimization status .......: converged
   2019-04-02T14:24:57:  Number of parameters ......: 10
   2019-04-02T14:24:57:  Number of free parameters .: 4
   2019-04-02T14:24:57:  Number of iterations ......: 2
   2019-04-02T14:24:57:  Lambda ....................: 1e-05
   2019-04-02T14:24:57:  Maximum log likelihood ....: -57887.417
   2019-04-02T14:24:57:  Observed events  (Nobs) ...: 19137.000
   2019-04-02T14:24:57:  Predicted events (Npred) ..: 19136.996 (Nobs - Npred = 0.00354148293627077)
   2019-04-02T14:24:57: === GModels ===
   2019-04-02T14:24:57:  Number of models ..........: 2
   2019-04-02T14:24:57:  Number of parameters ......: 10
   2019-04-02T14:24:57: === GModelSky ===
   2019-04-02T14:24:57:  Name ......................: Crab
   2019-04-02T14:24:57:  Instruments ...............: all
   2019-04-02T14:24:57:  Observation identifiers ...: all
   2019-04-02T14:24:57:  Model type ................: PointSource
   2019-04-02T14:24:57:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
   2019-04-02T14:24:57:  Number of parameters ......: 6
   2019-04-02T14:24:57:  Number of spatial par's ...: 2
   2019-04-02T14:24:57:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
   2019-04-02T14:24:57:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
   2019-04-02T14:24:57:  Number of spectral par's ..: 3
   2019-04-02T14:24:57:   Prefactor ................: 5.52559284054621e-16 +/- 9.88229994960437e-18 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
   2019-04-02T14:24:57:   Index ....................: -2.48163444213634 +/- 0.015305403980771 [-0,-5]  (free,scale=-1,gradient)
   2019-04-02T14:24:57:   PivotEnergy ..............: 300000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
   2019-04-02T14:24:57:  Number of temporal par's ..: 1
   2019-04-02T14:24:57:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2019-04-02T14:24:57:  Number of scale par's .....: 0
   2019-04-02T14:24:57: === GCTAModelCubeBackground ===
   2019-04-02T14:24:57:  Name ......................: BackgroundModel
   2019-04-02T14:24:57:  Instruments ...............: CTA, HESS, MAGIC, VERITAS
   2019-04-02T14:24:57:  Observation identifiers ...: all
   2019-04-02T14:24:57:  Model type ................: "PowerLaw" * "Constant"
   2019-04-02T14:24:57:  Number of parameters ......: 4
   2019-04-02T14:24:57:  Number of spectral par's ..: 3
   2019-04-02T14:24:57:   Prefactor ................: 1.00540991217377 +/- 0.0157241034891596 [0.01,100] ph/cm2/s/MeV (free,scale=1,gradient)
   2019-04-02T14:24:57:   Index ....................: 0.00380886384530723 +/- 0.00942814666809632 [-5,5]  (free,scale=1,gradient)
   2019-04-02T14:24:57:   PivotEnergy ..............: 1000000 MeV (fixed,scale=1000000,gradient)
   2019-04-02T14:24:57:  Number of temporal par's ..: 1
   2019-04-02T14:24:57:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

You can also perform an unbinned maximum likelihood analysis taking the energy
dispersion into account. In that case the energy dispersion information will be
directly determined from the
:ref:`instrument response functions <glossary_irf>`
and no energy dispersion cube is required:

.. code-block:: bash

   $ ctlike edisp=yes
   Input event list, counts cube or observation definition XML file [cntcube_edisp.fits] selected_events_edisp.fits
   Calibration database [prod2]
   Instrument response function [South_0.5h]
   Input model definition XML file [models.xml] $CTOOLS/share/models/crab.xml
   Output model definition XML file [crab_results_edisp.xml] 

Here the output in the log file:

.. code-block:: none

   2019-04-02T14:26:58: +=================================+
   2019-04-02T14:26:58: | Maximum likelihood optimisation |
   2019-04-02T14:26:58: +=================================+
   2019-04-02T14:27:00:  >Iteration   0: -logL=143060.165, Lambda=1.0e-03
   2019-04-02T14:27:02:  >Iteration   1: -logL=143059.529, Lambda=1.0e-03, delta=0.636, step=1.0e+00, max(|grad|)=-0.938341 [Prefactor:6]
   2019-04-02T14:27:04:  >Iteration   2: -logL=143059.529, Lambda=1.0e-04, delta=0.000, step=1.0e+00, max(|grad|)=-0.002461 [Index:7]
   2019-04-02T14:27:05:
   2019-04-02T14:27:05: +=========================================+
   2019-04-02T14:27:05: | Maximum likelihood optimisation results |
   2019-04-02T14:27:05: +=========================================+
   2019-04-02T14:27:05: === GOptimizerLM ===
   2019-04-02T14:27:05:  Optimized function value ..: 143059.529
   2019-04-02T14:27:05:  Absolute precision ........: 0.005
   2019-04-02T14:27:05:  Acceptable value decrease .: 2
   2019-04-02T14:27:05:  Optimization status .......: converged
   2019-04-02T14:27:05:  Number of parameters ......: 10
   2019-04-02T14:27:05:  Number of free parameters .: 4
   2019-04-02T14:27:05:  Number of iterations ......: 2
   2019-04-02T14:27:05:  Lambda ....................: 1e-05
   2019-04-02T14:27:05:  Maximum log likelihood ....: -143059.529
   2019-04-02T14:27:05:  Observed events  (Nobs) ...: 22407.000
   2019-04-02T14:27:05:  Predicted events (Npred) ..: 22406.999 (Nobs - Npred = 0.000615512442891486)
   2019-04-02T14:27:05: === GModels ===
   2019-04-02T14:27:05:  Number of models ..........: 2
   2019-04-02T14:27:05:  Number of parameters ......: 10
   2019-04-02T14:27:05: === GModelSky ===
   2019-04-02T14:27:05:  Name ......................: Crab
   2019-04-02T14:27:05:  Instruments ...............: all
   2019-04-02T14:27:05:  Observation identifiers ...: all
   2019-04-02T14:27:05:  Model type ................: PointSource
   2019-04-02T14:27:05:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
   2019-04-02T14:27:05:  Number of parameters ......: 6
   2019-04-02T14:27:05:  Number of spatial par's ...: 2
   2019-04-02T14:27:05:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
   2019-04-02T14:27:05:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
   2019-04-02T14:27:05:  Number of spectral par's ..: 3
   2019-04-02T14:27:05:   Prefactor ................: 5.61701723486666e-16 +/- 1.00237442767986e-17 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
   2019-04-02T14:27:05:   Index ....................: -2.48356027289941 +/- 0.0152626158555975 [-0,-5]  (free,scale=-1,gradient)
   2019-04-02T14:27:05:   PivotEnergy ..............: 300000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
   2019-04-02T14:27:05:  Number of temporal par's ..: 1
   2019-04-02T14:27:05:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2019-04-02T14:27:05:  Number of scale par's .....: 0
   2019-04-02T14:27:05: === GCTAModelIrfBackground ===
   2019-04-02T14:27:05:  Name ......................: CTABackgroundModel
   2019-04-02T14:27:05:  Instruments ...............: CTA
   2019-04-02T14:27:05:  Observation identifiers ...: all
   2019-04-02T14:27:05:  Model type ................: "PowerLaw" * "Constant"
   2019-04-02T14:27:05:  Number of parameters ......: 4
   2019-04-02T14:27:05:  Number of spectral par's ..: 3
   2019-04-02T14:27:05:   Prefactor ................: 1.0079566950951 +/- 0.0133706835965654 [0.001,1000] ph/cm2/s/MeV (free,scale=1,gradient)
   2019-04-02T14:27:05:   Index ....................: 0.00303753293882809 +/- 0.00807470190154737 [-5,5]  (free,scale=1,gradient)
   2019-04-02T14:27:05:   PivotEnergy ..............: 1000000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
   2019-04-02T14:27:05:  Number of temporal par's ..: 1
   2019-04-02T14:27:05:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)


