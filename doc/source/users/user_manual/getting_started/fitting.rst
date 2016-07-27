.. _sec_fitting_cta:

Fitting CTA data
~~~~~~~~~~~~~~~~

Now you are ready to fit the simulated data with a model. For simplicity
the same model is used in the example that you used to simulate the data
with :ref:`ctobssim`. Model fitting is done using the :ref:`ctlike` tool,
and you do the fit by typing:

.. code-block:: bash

  $ ctlike
  Input event list, counts cube or observation definition XML file [events.fits] cntcube.fits
  Input exposure cube file (only needed for stacked analysis) [NONE] 
  Input PSF cube file (only needed for stacked analysis) [NONE] 
  Input background cube file (only needed for stacked analysis) [NONE] 
  Calibration database [prod2] 
  Instrument response function [South_0.5h] 
  Input model XML file [$CTOOLS/share/models/crab.xml] 
  Output model XML file [crab_results.xml] 

The data are fitted in *binned* mode, which means that the events
have been binned into a counts cube and the fit computes the log-likelihood
function by summing over all 200 x 200 x 20 bins of the counts cube. There is
an alternative method, the so called *unbinned* mode, where the events are
not binned into a counts cube and the log-likelihood is computed directly by
summing over all events. We will explore the *unbinned* mode later.

One of the parameters of :ref:`ctlike` is a source model output file
(we specified ``crab_results.xml`` in the example), and this file will be
a copy of the source model input XML file where the parameter values have
been replaced by the fit results. In addition, the statistical uncertainties
are added for each fitted parameter using the attribute ``error``.
Below we show the XML result file that has been produced by the run:

.. code-block:: xml

  <?xml version="1.0" encoding="UTF-8" standalone="no"?>
  <source_library title="source library">
    <source name="Crab" type="PointSource">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor" value="5.79064" error="0.101828" scale="1e-16" min="1e-07" max="1000" free="1" />
        <parameter name="Index" value="2.47177" error="0.0154935" scale="-1" min="0" max="5" free="1" />
        <parameter name="PivotEnergy" value="0.3" scale="1e+06" min="0.01" max="1000" free="0" />
      </spectrum>
      <spatialModel type="SkyDirFunction">
        <parameter name="RA" value="83.6331" scale="1" min="-360" max="360" free="0" />
        <parameter name="DEC" value="22.0145" scale="1" min="-90" max="90" free="0" />
      </spatialModel>
    </source>
    <source name="CTABackgroundModel" type="CTAIrfBackground" instrument="CTA">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor" value="1.01642" error="0.0162226" scale="1" min="0.001" max="1000" free="1" />
        <parameter name="Index" value="0.00462945" error="0.00964737" scale="1" min="-5" max="5" free="1" />
        <parameter name="PivotEnergy" value="1" scale="1e+06" min="0.01" max="1000" free="0" />
      </spectrum>
    </source>
  </source_library>

In this example, the ``Prefactor`` and ``Index`` of the spectral model for the
Crab nebula as well as the ``Prefactor`` and ``Index`` of the background spectral
model have been fitted (all parameters having the attribute ``free="1"`` are
fitted). Now you see also the effect of having multiplied the background
model with a power law. In that way, the amplitude of the background as
well as it's spectral slope is adjusted by the fit. Obviously, in this
example the adjustment compensates only for the statistical fluctuations
of the background, but with real data, the adjustment may account also for
some of the systematic uncertainties.

.. warning::

   As good practice, the amplitude of the background model should always be
   left as a free parameter of the fit. Otherwise, any uncertainty in the
   background rate will immediately propagate into the flux estimate of the 
   source. 

To get more details about the model fitting you can inspect the log file.
Below the last lines of the ctlike.log log file that has been produced by
this run:

.. code-block:: none

  2016-06-29T16:11:05: +=================================+
  2016-06-29T16:11:05: | Maximum likelihood optimisation |
  2016-06-29T16:11:05: +=================================+
  2016-06-29T16:11:06:  >Iteration   0: -logL=55574.937, Lambda=1.0e-03
  2016-06-29T16:11:06:  >Iteration   1: -logL=55573.383, Lambda=1.0e-03, delta=1.555, max(|grad|)=3.170837 [Index:7]
  2016-06-29T16:11:07:  >Iteration   2: -logL=55573.382, Lambda=1.0e-04, delta=0.001, max(|grad|)=0.056280 [Index:3]
  ...
  2016-06-29T16:11:08: +=========================================+
  2016-06-29T16:11:08: | Maximum likelihood optimisation results |
  2016-06-29T16:11:08: +=========================================+
  2016-06-29T16:11:08: === GOptimizerLM ===
  2016-06-29T16:11:08:  Optimized function value ..: 55573.382
  2016-06-29T16:11:08:  Absolute precision ........: 0.005
  2016-06-29T16:11:08:  Acceptable value decrease .: 2
  2016-06-29T16:11:08:  Optimization status .......: converged
  2016-06-29T16:11:08:  Number of parameters ......: 10
  2016-06-29T16:11:08:  Number of free parameters .: 4
  2016-06-29T16:11:08:  Number of iterations ......: 2
  2016-06-29T16:11:08:  Lambda ....................: 1e-05
  2016-06-29T16:11:08:  Maximum log likelihood ....: -55573.382
  2016-06-29T16:11:08:  Observed events  (Nobs) ...: 18532.000
  2016-06-29T16:11:08:  Predicted events (Npred) ..: 18531.998 (Nobs - Npred = 0.00243959)
  2016-06-29T16:11:08: === GModels ===
  2016-06-29T16:11:08:  Number of models ..........: 2
  2016-06-29T16:11:08:  Number of parameters ......: 10
  2016-06-29T16:11:08: === GModelSky ===
  2016-06-29T16:11:08:  Name ......................: Crab
  2016-06-29T16:11:08:  Instruments ...............: all
  2016-06-29T16:11:08:  Instrument scale factors ..: unity
  2016-06-29T16:11:08:  Observation identifiers ...: all
  2016-06-29T16:11:08:  Model type ................: PointSource
  2016-06-29T16:11:08:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
  2016-06-29T16:11:08:  Number of parameters ......: 6
  2016-06-29T16:11:08:  Number of spatial par's ...: 2
  2016-06-29T16:11:08:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
  2016-06-29T16:11:08:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
  2016-06-29T16:11:08:  Number of spectral par's ..: 3
  2016-06-29T16:11:08:   Prefactor ................: 5.79064e-16 +/- 1.01828e-17 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
  2016-06-29T16:11:08:   Index ....................: -2.47177 +/- 0.0154935 [-0,-5]  (free,scale=-1,gradient)
  2016-06-29T16:11:08:   PivotEnergy ..............: 300000 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2016-06-29T16:11:08:  Number of temporal par's ..: 1
  2016-06-29T16:11:08:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
  2016-06-29T16:11:08: === GCTAModelIrfBackground ===
  2016-06-29T16:11:08:  Name ......................: CTABackgroundModel
  2016-06-29T16:11:08:  Instruments ...............: CTA
  2016-06-29T16:11:08:  Instrument scale factors ..: unity
  2016-06-29T16:11:08:  Observation identifiers ...: all
  2016-06-29T16:11:08:  Model type ................: "PowerLaw" * "Constant"
  2016-06-29T16:11:08:  Number of parameters ......: 4
  2016-06-29T16:11:08:  Number of spectral par's ..: 3
  2016-06-29T16:11:08:   Prefactor ................: 1.01642 +/- 0.0162226 [0.001,1000] ph/cm2/s/MeV (free,scale=1,gradient)
  2016-06-29T16:11:08:   Index ....................: 0.00462945 +/- 0.00964737 [-5,5]  (free,scale=1,gradient)
  2016-06-29T16:11:08:   PivotEnergy ..............: 1e+06 [10000,1e+09] MeV (fixed,scale=1e+06,gradient)
  2016-06-29T16:11:08:  Number of temporal par's ..: 1
  2016-06-29T16:11:08:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

The maximum likelihood optimizer required 2 iterations to converge. This
is pretty fast, but recall that you used the same model file for the simulation
and for fitting, hence the initial parameter values were already very close
to the best fitting values. To see the impact of the initial parameters on
the fit result, you may re-run :ref:`ctlike` using another copy of the model
XML file where you change the value attributes of the parameters that should be 
fitted. You will see that the optimizer requires a couple of more iterations,
but it should converge to the same solution (provided that the initial values
are not too far from the best fitting values).

.. note::

   As sanity check you should verify that the predicted number of events
   (Npred) is equal to the observed number of events (Nobs). To facilitate
   this comparison, :ref:`ctlike` provides the difference Nobs - Npred in 
   the log file. In real life situations, this difference may not always be
   small, in particular if the source model is too constrained. You may 
   then free some of the model parameters so that the fit can correctly
   describe the data.

.. note::

   The :ref:`ctlike` tool has the ability to estimate the detection 
   significance for sources in the XML model. This is done by computing
   the Test Statistic value which is defined as twice the log-likelihood
   difference between fitting a source at a given position on top of a 
   (background) model or fitting no source. Roughly speaken, the square
   root of the Test Statistic value gives the source detection significance
   in Gaussian sigmas, although the exact relation depends somewhat on
   the formulation of the statistical problem.

   To instruct :ref:`ctlike` to compute the Test Statistic value for a
   given source you need to add the attribute ``tscalc="1"`` to the XML
   file:

   .. code-block:: xml

      <source name="Crab" type="PointSource" tscalc="1">

   :ref:`ctlike` will then compute the Test Statistic value for that
   source and dump the result in the log file:

   .. code-block:: none

      2016-06-29T16:17:38: === GModelSky ===
      2016-06-29T16:17:38:  Name ......................: Crab
      2016-06-29T16:17:38:  Instruments ...............: all
      2016-06-29T16:17:38:  Test Statistic ............: 21778.8

   The Test Statistic value will also be added as new attribute
   ``ts`` to the XML result file:

   .. code-block:: xml

      <source name="Crab" type="PointSource" ts="21778.780" tscalc="1">
