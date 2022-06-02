.. _comptel_fitting:

Fitting the data
----------------

  .. admonition:: What you will learn

     You will learn how you fit a model to COMPTEL data.


Now you are prepared to fit the model to COMPTEL data. You do this using
the :ref:`comlixfit` script. In the example below the script is run with the
``debug=yes`` option so that the fit results are directly displayed in the
terminal. :ref:`comlixfit` performs an iterative background fitting using the
SRCLIX approach; in particuliar, the BGDLIXE method is used for background
model computation.

The fit took 3 SRCLIX iterations to converge. The Crab is detected with
a test statistic value of 1074.78. The fit determined the position of the
gamma-ray source and the spectral power law parameters, which are the intensity
at a pivot energy of 1 MeV and the spectral index. As a quality check, the
fit predicted 517813.013 events in the model which is close to the observed
number of 517813 events.

.. code-block:: bash

   $ comlixfit debug=yes
   Input observation definition file [obs.xml] obs_binned.xml
   Input model definition file [models.xml]
   Method for background computation (PHINOR|BGDLIXA|BGDLIXE) [BGDLIXE]
   Output observation definition file [outobs.xml]
   Output model definition file [results.xml]
   ...
   2022-03-17T13:20:40: +====================+
   2022-03-17T13:20:40: | Input observations |
   2022-03-17T13:20:40: +====================+
   2022-03-17T13:20:40: === GObservations ===
   2022-03-17T13:20:40:  Number of observations ....: 4
   2022-03-17T13:20:40:  Number of models ..........: 5
   2022-03-17T13:20:40:  Number of observed events .: 517813
   2022-03-17T13:20:40:  Number of predicted events : 0
   2022-03-17T13:20:40:
   2022-03-17T13:20:40: +============================================+
   2022-03-17T13:20:40: | Iterative maximum likelihood model fitting |
   2022-03-17T13:20:40: +============================================+
   2022-03-17T13:20:48:  logL after iteration 1 ....: -28465.83123
   2022-03-17T13:20:55:  logL after iteration 2 ....: -28466.41211 (0.58088)
   2022-03-17T13:21:00:  logL after iteration 3 ....: -28466.29714 (-0.11498)
   2022-03-17T13:21:09:  logL after final iteration : -28466.41290 (0.00079)
   2022-03-17T13:21:09:
   2022-03-17T13:21:09: +=========================================+
   2022-03-17T13:21:09: | Maximum likelihood optimisation results |
   2022-03-17T13:21:09: +=========================================+
   2022-03-17T13:21:09: === GOptimizerLM ===
   2022-03-17T13:21:09:  Optimized function value ..: -28466.413
   2022-03-17T13:21:09:  Absolute precision ........: 0.005
   2022-03-17T13:21:09:  Acceptable value decrease .: 0
   2022-03-17T13:21:09:  Optimization status .......: converged
   2022-03-17T13:21:09:  Number of parameters ......: 106
   2022-03-17T13:21:09:  Number of free parameters .: 99
   2022-03-17T13:21:09:  Number of iterations ......: 4
   2022-03-17T13:21:09:  Lambda ....................: 0.1
   2022-03-17T13:21:09:  Total number of iterations : 14
   2022-03-17T13:21:09:  Maximum log likelihood ....: 28466.4129009
   2022-03-17T13:21:09:  Observed events  (Nobs) ...: 517813.0
   2022-03-17T13:21:09:  Predicted events (Npred) ..: 517813.012556
   2022-03-17T13:21:09:  Nobs - Npred ..............: -0.012555927271
   2022-03-17T13:21:09: === GModels ===
   2022-03-17T13:21:09:  Number of models ..........: 5
   2022-03-17T13:21:09:  Number of parameters ......: 106
   2022-03-17T13:21:09: === GModelSky ===
   2022-03-17T13:21:09:  Name ......................: Crab
   2022-03-17T13:21:09:  Instruments ...............: all
   2022-03-17T13:21:09:  Test Statistic ............: 1074.78358831939
   2022-03-17T13:21:09:  Observation identifiers ...: all
   2022-03-17T13:21:09:  Model type ................: PointSource
   2022-03-17T13:21:09:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
   2022-03-17T13:21:09:  Number of parameters ......: 6
   2022-03-17T13:21:09:  Number of spatial par's ...: 2
   2022-03-17T13:21:09:   RA .......................: 83.7688330639726 +/- 0.111131701366959 deg (free,scale=1)
   2022-03-17T13:21:09:   DEC ......................: 21.6073511262309 +/- 0.101155131386771 deg (free,scale=1)
   2022-03-17T13:21:09:  Number of spectral par's ..: 3
   2022-03-17T13:21:09:   Prefactor ................: 0.00213237276990501 +/- 9.60239430063358e-05 [1e-25,infty[ ph/cm2/s/MeV (free,scale=0.002,gradient)
   2022-03-17T13:21:09:   Index ....................: -2.1249932516377 +/- 0.0369030627304123 [-10,10]  (free,scale=-2,gradient)
   2022-03-17T13:21:09:   PivotEnergy ..............: 1 MeV (fixed,scale=1,gradient)
   2022-03-17T13:21:09:  Number of temporal par's ..: 1
   2022-03-17T13:21:09:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2022-03-17T13:21:09:  Number of scale par's .....: 0

.. note::

   To shorten the output, omitted lines were replaced by three dots.
