.. _start_fitting:

Fitting binned data
-------------------

  .. admonition:: What you will learn

     You will learn how to **adjust a parametrised model to the counts cube
     by means of a maximum likelihood fit**.

Now you are ready to fit the data with a model. For simplicity the same model
is used in the example that you used to simulate the data with :ref:`ctobssim`.
Model fitting is done using the :ref:`ctlike` tool, and you do the fit by
typing:

.. code-block:: bash

   $ ctlike
   Input event list, counts cube or observation definition XML file [events.fits] cntcube.fits
   Input exposure cube file [NONE] expcube.fits
   Input PSF cube file [NONE] psfcube.fits
   Input background cube file [NONE] bkgcube.fits
   Input model definition XML file [$CTOOLS/share/models/crab.xml] models.xml
   Output model definition XML file [crab_results.xml]

The data are fitted in *binned* mode, which means that the events
have been binned into a counts cube and the fit computes the log-likelihood
function by summing over all 200 x 200 x 20 bins of the counts cube. There is
an alternative method, the so called *unbinned* mode, where the events are
not binned into a counts cube and the log-likelihood is computed directly by
summing over all events.
:ref:`We will explore the unbinned mode later <start_unbinned>`.

The :ref:`ctlike` tool produced a
:ref:`model definition XML output file <glossary_moddef>`
``crab_results.xml`` that contains the fit results.
The file is a copy of the
:ref:`model definition XML input file <glossary_moddef>`
``models.xml`` where the parameter values were replaced by the fit results.
In addition, the statistical uncertainties were added for each fitted parameter
using the attribute ``error``.
Below the
:ref:`model definition XML output file <glossary_moddef>`
that was produced by the run:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Crab" type="PointSource">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor" value="5.61368377731832" error="0.0982975197293023" scale="1e-16" min="1e-07" max="1000" free="1" />
         <parameter name="Index" value="2.45569249143075" error="0.0145471260625612" scale="-1" min="0" max="5" free="1" />
         <parameter name="PivotEnergy" value="0.3" scale="1000000" min="0.01" max="1000" free="0" />
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA" value="83.6331" scale="1" min="-360" max="360" free="0" />
         <parameter name="DEC" value="22.0145" scale="1" min="-90" max="90" free="0" />
       </spatialModel>
     </source>
     <source name="BackgroundModel" type="CTACubeBackground" instrument="CTA,HESS,MAGIC,VERITAS">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor" value="0.945994938701748" error="0.0148673466828289" scale="1" min="0.01" max="100" free="1" />
         <parameter name="Index" value="0.0198073832868029" error="0.00948925246912034" scale="1" min="-5" max="5" free="1" />
         <parameter name="PivotEnergy" value="1" scale="1000000" free="0" />
       </spectrum>
     </source>
  </source_library>

In this example, the ``Prefactor`` and ``Index`` of the spectral model for the
Crab nebula as well as the ``Prefactor`` and ``Index`` of the background spectral
model were fitted (these parameters had the attribute ``free="1"``). Now you
see also the effect of having multiplied the background model with a power law.
In that way, the amplitude of the background as well as it's spectral slope is
adjusted by the fit. Obviously, in this example the adjustment compensates only
for the statistical fluctuations of the background, but with real data, the
adjustment may account also for some of the systematic uncertainties.

  .. warning::

     As good practice, the amplitude of the background model should always be
     left as a free parameter of the fit. Otherwise, any uncertainty in the
     background rate will immediately propagate into the flux estimate of the
     source.

  .. warning::

     You may have recognized the ``scale`` and ``value`` attributes in the
     :ref:`model definition XML file <glossary_moddef>`. The value of each
     parameter is obtained by multiplying ``value`` with ``scale``. This allows
     for a pre-scaling of the parameters, and **you should make use of this
     capability to have the value attributes of all parameters that are fitted
     of about the same order, typically 1**. This is necessary to assure a
     proper convergence of the fitting algorithm.

To get more details about the model fitting you can inspect the log file.
Below the last lines of the log file that was produced by this run:

.. code-block:: none

   2017-08-08T20:24:06: +=================================+
   2017-08-08T20:24:06: | Maximum likelihood optimisation |
   2017-08-08T20:24:06: +=================================+
   2017-08-08T20:24:07:  >Iteration   0: -logL=56488.484, Lambda=1.0e-03
   2017-08-08T20:24:07:  >Iteration   1: -logL=56431.025, Lambda=1.0e-03, delta=57.459, step=1.0e+00, max(|grad|)=128.772103 [Index:7]
   2017-08-08T20:24:08:  >Iteration   2: -logL=56430.794, Lambda=1.0e-04, delta=0.232, step=1.0e+00, max(|grad|)=0.505741 [Index:7]
   2017-08-08T20:24:09:  >Iteration   3: -logL=56430.794, Lambda=1.0e-05, delta=0.000, step=1.0e+00, max(|grad|)=0.017423 [Index:3]
   2017-08-08T20:24:09:
   2017-08-08T20:24:09: +=========================================+
   2017-08-08T20:24:09: | Maximum likelihood optimisation results |
   2017-08-08T20:24:09: +=========================================+
   2017-08-08T20:24:09: === GOptimizerLM ===
   2017-08-08T20:24:09:  Optimized function value ..: 56430.794
   2017-08-08T20:24:09:  Absolute precision ........: 0.005
   2017-08-08T20:24:09:  Acceptable value decrease .: 2
   2017-08-08T20:24:09:  Optimization status .......: converged
   2017-08-08T20:24:09:  Number of parameters ......: 10
   2017-08-08T20:24:09:  Number of free parameters .: 4
   2017-08-08T20:24:09:  Number of iterations ......: 3
   2017-08-08T20:24:09:  Lambda ....................: 1e-06
   2017-08-08T20:24:09:  Maximum log likelihood ....: -56430.794
   2017-08-08T20:24:09:  Observed events  (Nobs) ...: 18827.000
   2017-08-08T20:24:09:  Predicted events (Npred) ..: 18827.000 (Nobs - Npred = 1.92705083463807e-05)
   2017-08-08T20:24:09: === GModels ===
   2017-08-08T20:24:09:  Number of models ..........: 2
   2017-08-08T20:24:09:  Number of parameters ......: 10
   2017-08-08T20:24:09: === GModelSky ===
   2017-08-08T20:24:09:  Name ......................: Crab
   2017-08-08T20:24:09:  Instruments ...............: all
   2017-08-08T20:24:09:  Instrument scale factors ..: unity
   2017-08-08T20:24:09:  Observation identifiers ...: all
   2017-08-08T20:24:09:  Model type ................: PointSource
   2017-08-08T20:24:09:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
   2017-08-08T20:24:09:  Number of parameters ......: 6
   2017-08-08T20:24:09:  Number of spatial par's ...: 2
   2017-08-08T20:24:09:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
   2017-08-08T20:24:09:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
   2017-08-08T20:24:09:  Number of spectral par's ..: 3
   2017-08-08T20:24:09:   Prefactor ................: 5.61368377731832e-16 +/- 9.82975197293023e-18 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
   2017-08-08T20:24:09:   Index ....................: -2.45569249143075 +/- 0.0145471260625612 [-0,-5]  (free,scale=-1,gradient)
   2017-08-08T20:24:09:   PivotEnergy ..............: 300000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
   2017-08-08T20:24:09:  Number of temporal par's ..: 1
   2017-08-08T20:24:09:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2017-08-08T20:24:09: === GCTAModelCubeBackground ===
   2017-08-08T20:24:09:  Name ......................: BackgroundModel
   2017-08-08T20:24:09:  Instruments ...............: CTA, HESS, MAGIC, VERITAS
   2017-08-08T20:24:09:  Instrument scale factors ..: unity
   2017-08-08T20:24:09:  Observation identifiers ...: all
   2017-08-08T20:24:09:  Model type ................: "PowerLaw" * "Constant"
   2017-08-08T20:24:09:  Number of parameters ......: 4
   2017-08-08T20:24:09:  Number of spectral par's ..: 3
   2017-08-08T20:24:09:   Prefactor ................: 0.945994938701748 +/- 0.0148673466828289 [0.01,100] ph/cm2/s/MeV (free,scale=1,gradient)
   2017-08-08T20:24:09:   Index ....................: 0.0198073832868029 +/- 0.00948925246912034 [-5,5]  (free,scale=1,gradient)
   2017-08-08T20:24:09:   PivotEnergy ..............: 1000000 MeV (fixed,scale=1000000,gradient)
   2017-08-08T20:24:09:  Number of temporal par's ..: 1
   2017-08-08T20:24:09:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

The maximum likelihood optimizer required 3 iterations to converge. This
is pretty fast, but recall that you used the same model file for the simulation
and for fitting, hence the initial parameter values were already very close
to the best fitting values. To see the impact of the initial parameters on
the fit result, you may re-run :ref:`ctlike` using another copy of the
:ref:`model definition XML input file <glossary_moddef>`
where you change the value attributes of the parameters that should be
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

        2017-08-08T20:34:29: === GModelSky ===
        2017-08-08T20:34:29:  Name ......................: Crab
        2017-08-08T20:34:29:  Instruments ...............: all
        2017-08-08T20:34:29:  Test Statistic ............: 21873.9558591403

     The Test Statistic value will also be added as new attribute
     ``ts`` to the XML result file:

     .. code-block:: xml

        <source name="Crab" type="PointSource" ts="21873.956" tscalc="1">
