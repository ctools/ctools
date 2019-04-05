.. _sec_iact_analysis:

Performing a maximum likelihood fit
===================================

  .. admonition:: What you will learn

     You will learn how to perform a maximum likelihood fit of IACT data.

Following a :ref:`csiactobs` run with the ``inmodel=$CTOOLS/share/models/crab.xml``
you have a model that is ready for maximum likelihood fitting.
You do such a fit using

.. code-block:: bash

   $ ctlike
   Input event list, counts cube or observation definition XML file [obs_selected.xml]
   Input model definition XML file [models.xml]
   Output model definition XML file [crab_results.xml]

which will produce the following output in the log file ``ctlike.log``:

.. code-block:: none

   2019-04-05T14:18:17: +=================================+
   2019-04-05T14:18:17: | Maximum likelihood optimisation |
   2019-04-05T14:18:17: +=================================+
   2019-04-05T14:18:17:    Parameter "Prefactor" has zero curvature. Fix parameter.
   2019-04-05T14:18:17:    Parameter "Index" has zero curvature. Fix parameter.
   2019-04-05T14:18:17:  >Iteration   0: -logL=203790.772, Lambda=1.0e-03
   2019-04-05T14:18:17:  >Iteration   1: -logL=192855.522, Lambda=1.0e-03, delta=10935.250, step=1.0e+00, max(|grad|)=-2040.622705 [Prefactor:0]
   2019-04-05T14:18:17:  >Iteration   2: -logL=182999.110, Lambda=1.0e-04, delta=9856.412, step=1.0e+00, max(|grad|)=-981.683624 [Prefactor:0]
   2019-04-05T14:18:17:  >Iteration   3: -logL=174691.964, Lambda=1.0e-05, delta=8307.146, step=1.0e+00, max(|grad|)=620.662474 [Index:9]
   2019-04-05T14:18:17:  >Iteration   4: -logL=168496.493, Lambda=1.0e-06, delta=6195.471, step=1.0e+00, max(|grad|)=587.098865 [Index:9]
   2019-04-05T14:18:17:  >Iteration   5: -logL=164898.487, Lambda=1.0e-07, delta=3598.006, step=1.0e+00, max(|grad|)=423.800306 [Index:5]
   2019-04-05T14:18:18:  >Iteration   6: -logL=163621.925, Lambda=1.0e-08, delta=1276.562, step=1.0e+00, max(|grad|)=237.546032 [Index:1]
   2019-04-05T14:18:18:  >Iteration   7: -logL=163423.407, Lambda=1.0e-09, delta=198.518, step=1.0e+00, max(|grad|)=102.646722 [Index:1]
   2019-04-05T14:18:18:  >Iteration   8: -logL=163409.717, Lambda=1.0e-10, delta=13.690, step=1.0e+00, max(|grad|)=35.371246 [Index:1]
   2019-04-05T14:18:18:  >Iteration   9: -logL=163408.542, Lambda=1.0e-11, delta=1.175, step=1.0e+00, max(|grad|)=11.708402 [Index:1]
   2019-04-05T14:18:18:  >Iteration  10: -logL=163408.416, Lambda=1.0e-12, delta=0.126, step=1.0e+00, max(|grad|)=3.907383 [Index:1]
   2019-04-05T14:18:18:  >Iteration  11: -logL=163408.402, Lambda=1.0e-13, delta=0.014, step=1.0e+00, max(|grad|)=1.309476 [Index:1]
   2019-04-05T14:18:18:  >Iteration  12: -logL=163408.400, Lambda=1.0e-14, delta=0.002, step=1.0e+00, max(|grad|)=0.439462 [Index:1]
   2019-04-05T14:18:18:    Free parameter "Prefactor" after convergence was reached with frozen parameter.
   2019-04-05T14:18:18:    Free parameter "Index" after convergence was reached with frozen parameter.
   2019-04-05T14:18:18:
   2019-04-05T14:18:18: +=========================================+
   2019-04-05T14:18:18: | Maximum likelihood optimisation results |
   2019-04-05T14:18:18: +=========================================+
   2019-04-05T14:18:18: === GOptimizerLM ===
   2019-04-05T14:18:18:  Optimized function value ..: 163408.400
   2019-04-05T14:18:18:  Absolute precision ........: 0.005
   2019-04-05T14:18:18:  Acceptable value decrease .: 2
   2019-04-05T14:18:18:  Optimization status .......: converged
   2019-04-05T14:18:18:  Number of parameters ......: 26
   2019-04-05T14:18:18:  Number of free parameters .: 12
   2019-04-05T14:18:18:  Number of iterations ......: 12
   2019-04-05T14:18:18:  Lambda ....................: 1e-15
   2019-04-05T14:18:18:  Maximum log likelihood ....: -163408.400
   2019-04-05T14:18:18:  Observed events  (Nobs) ...: 18641.000
   2019-04-05T14:18:18:  Predicted events (Npred) ..: 18640.998 (Nobs - Npred = 0.00164345248049358)
   2019-04-05T14:18:18: === GModels ===
   2019-04-05T14:18:18:  Number of models ..........: 6
   2019-04-05T14:18:18:  Number of parameters ......: 26
   2019-04-05T14:18:18: === GCTAModelAeffBackground ===
   2019-04-05T14:18:18:  Name ......................: bkg_23523
   2019-04-05T14:18:18:  Instruments ...............: HESS
   2019-04-05T14:18:18:  Instrument scale factors ..: unity
   2019-04-05T14:18:18:  Observation identifiers ...: 23523
   2019-04-05T14:18:18:  Model type ................: "PowerLaw" * "Constant"
   2019-04-05T14:18:18:  Number of parameters ......: 4
   2019-04-05T14:18:18:  Number of spectral par's ..: 3
   2019-04-05T14:18:18:   Prefactor ................: 2.98790819794319e-13 +/- 4.43400659819813e-15 [1e-16,1e-12] ph/cm2/s/MeV (free,scale=1e-14,gradient)
   2019-04-05T14:18:18:   Index ....................: -3.22702359177827 +/- 0.0221521174108544 [-5,5]  (free,scale=1,gradient)
   2019-04-05T14:18:18:   PivotEnergy ..............: 1000000 MeV (fixed,scale=1000000,gradient)
   2019-04-05T14:18:18:  Number of temporal par's ..: 1
   2019-04-05T14:18:18:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2019-04-05T14:18:18: === GCTAModelAeffBackground ===
   2019-04-05T14:18:18:  Name ......................: bkg_23526
   2019-04-05T14:18:18:  Instruments ...............: HESS
   2019-04-05T14:18:18:  Instrument scale factors ..: unity
   2019-04-05T14:18:18:  Observation identifiers ...: 23526
   2019-04-05T14:18:18:  Model type ................: "PowerLaw" * "Constant"
   2019-04-05T14:18:18:  Number of parameters ......: 4
   2019-04-05T14:18:18:  Number of spectral par's ..: 3
   2019-04-05T14:18:18:   Prefactor ................: 2.44289795381698e-13 +/- 3.75318255643637e-15 [1e-16,1e-12] ph/cm2/s/MeV (free,scale=1e-14,gradient)
   2019-04-05T14:18:18:   Index ....................: -3.21741632072522 +/- 0.0232805723886611 [-5,5]  (free,scale=1,gradient)
   2019-04-05T14:18:18:   PivotEnergy ..............: 1000000 MeV (fixed,scale=1000000,gradient)
   2019-04-05T14:18:18:  Number of temporal par's ..: 1
   2019-04-05T14:18:18:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2019-04-05T14:18:18: === GCTAModelAeffBackground ===
   2019-04-05T14:18:18:  Name ......................: bkg_23559
   2019-04-05T14:18:18:  Instruments ...............: HESS
   2019-04-05T14:18:18:  Instrument scale factors ..: unity
   2019-04-05T14:18:18:  Observation identifiers ...: 23559
   2019-04-05T14:18:18:  Model type ................: "PowerLaw" * "Constant"
   2019-04-05T14:18:18:  Number of parameters ......: 4
   2019-04-05T14:18:18:  Number of spectral par's ..: 3
   2019-04-05T14:18:18:   Prefactor ................: 2.26357632227583e-13 +/- 3.56983854227473e-15 [1e-16,1e-12] ph/cm2/s/MeV (free,scale=1e-14,gradient)
   2019-04-05T14:18:18:   Index ....................: -3.25328282552122 +/- 0.0239266148104449 [-5,5]  (free,scale=1,gradient)
   2019-04-05T14:18:18:   PivotEnergy ..............: 1000000 MeV (fixed,scale=1000000,gradient)
   2019-04-05T14:18:18:  Number of temporal par's ..: 1
   2019-04-05T14:18:18:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2019-04-05T14:18:18: === GCTAModelAeffBackground ===
   2019-04-05T14:18:18:  Name ......................: bkg_23592
   2019-04-05T14:18:18:  Instruments ...............: HESS
   2019-04-05T14:18:18:  Instrument scale factors ..: unity
   2019-04-05T14:18:18:  Observation identifiers ...: 23592
   2019-04-05T14:18:18:  Model type ................: "PowerLaw" * "Constant"
   2019-04-05T14:18:18:  Number of parameters ......: 4
   2019-04-05T14:18:18:  Number of spectral par's ..: 3
   2019-04-05T14:18:18:   Prefactor ................: 2.87721724170504e-13 +/- 4.31204359584396e-15 [1e-16,1e-12] ph/cm2/s/MeV (free,scale=1e-14,gradient)
   2019-04-05T14:18:18:   Index ....................: -3.25393054001656 +/- 0.0227802144778143 [-5,5]  (free,scale=1,gradient)
   2019-04-05T14:18:18:   PivotEnergy ..............: 1000000 MeV (fixed,scale=1000000,gradient)
   2019-04-05T14:18:18:  Number of temporal par's ..: 1
   2019-04-05T14:18:18:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2019-04-05T14:18:18: === GModelSky ===
   2019-04-05T14:18:18:  Name ......................: Crab
   2019-04-05T14:18:18:  Instruments ...............: all
   2019-04-05T14:18:18:  Instrument scale factors ..: unity
   2019-04-05T14:18:18:  Observation identifiers ...: all
   2019-04-05T14:18:18:  Model type ................: PointSource
   2019-04-05T14:18:18:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
   2019-04-05T14:18:18:  Number of parameters ......: 6
   2019-04-05T14:18:18:  Number of spatial par's ...: 2
   2019-04-05T14:18:18:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
   2019-04-05T14:18:18:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
   2019-04-05T14:18:18:  Number of spectral par's ..: 3
   2019-04-05T14:18:18:   Prefactor ................: 1.205595290581e-15 +/- 1.02399444013702e-16 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
   2019-04-05T14:18:18:   Index ....................: -2.64310065863216 +/- 0.0532926238637032 [-0,-5]  (free,scale=-1,gradient)
   2019-04-05T14:18:18:   PivotEnergy ..............: 300000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
   2019-04-05T14:18:18:  Number of temporal par's ..: 1
   2019-04-05T14:18:18:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2019-04-05T14:18:18: === GCTAModelIrfBackground ===
   2019-04-05T14:18:18:  Name ......................: CTABackgroundModel
   2019-04-05T14:18:18:  Instruments ...............: CTA
   2019-04-05T14:18:18:  Instrument scale factors ..: unity
   2019-04-05T14:18:18:  Observation identifiers ...: all
   2019-04-05T14:18:18:  Model type ................: "PowerLaw" * "Constant"
   2019-04-05T14:18:18:  Number of parameters ......: 4
   2019-04-05T14:18:18:  Number of spectral par's ..: 3
   2019-04-05T14:18:18:   Prefactor ................: 1 +/- 0 [0.001,1000] ph/cm2/s/MeV (free,scale=1,gradient)
   2019-04-05T14:18:18:   Index ....................: 0 +/- 0 [-5,5]  (free,scale=1,gradient)
   2019-04-05T14:18:18:   PivotEnergy ..............: 1000000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
   2019-04-05T14:18:18:  Number of temporal par's ..: 1
   2019-04-05T14:18:18:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

.. note::
   Using ``inmodel=$CTOOLS/share/models/crab.xml`` will also append the CTA
   background model to the ``models.xml`` file, yet since this model only
   applies to CTA data it is ignored in the fit. This is signalled by

   .. code-block:: none

      2019-04-05T14:18:17:    Parameter "Prefactor" has zero curvature. Fix parameter.
      2019-04-05T14:18:17:    Parameter "Index" has zero curvature. Fix parameter.

   during the fit, and the corresponding model parameter have an error of
   zero, although they are signalled as free.

.. warning::
   The effective area background model generated by :ref:`csiactobs` is very
   simplistic, and should not be used for a serious data analysis of H.E.S.S.
   data. Please follow the procedure :ref:`hess_dr1` for the generation of
   a more reliable background model.
