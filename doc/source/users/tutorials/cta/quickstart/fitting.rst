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
have been binned into a 3D counts cube and the fit computes the log-likelihood
function by summing over all 200 x 200 x 20 bins of the counts cube. There is
an alternative method, the so called *unbinned* mode, where the events are
not binned into a counts cube and the log-likelihood is computed directly by
summing over all events.
:ref:`We will explore the unbinned mode later <start_unbinned>`. Yet a
third alternative method (*On/Off* mode) consists in binning the events only in energy
integrating over arrival direction within a Source (On) and one or
more Background (Off) regions to minimise the dependency on the IRF
background model. :ref:`We will explore the On/Off mode later <start_onoff>`

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
         <parameter name="Prefactor" value="5.79620732342312" error="0.101207859506902" scale="1e-16" min="1e-07" max="1000" free="1" />
         <parameter name="Index" value="2.48767829376908" error="0.0150183262185825" scale="-1" min="0" max="5" free="1" />
         <parameter name="PivotEnergy" value="0.3" scale="1000000" min="0.01" max="1000" free="0" />
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA" value="83.6331" scale="1" min="-360" max="360" free="0" />
         <parameter name="DEC" value="22.0145" scale="1" min="-90" max="90" free="0" />
       </spatialModel>
     </source>
     <source name="BackgroundModel" type="CTACubeBackground" instrument="CTA,HESS,MAGIC,VERITAS">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor" value="0.997287677101311" error="0.0156090791023617" scale="1" min="0.01" max="100" free="1" />
         <parameter name="Index" value="-0.00905917021700502" error="0.00939297053758322" scale="1" min="-5" max="5" free="1" />
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

   2019-04-02T13:55:43: +=================================+
   2019-04-02T13:55:43: | Maximum likelihood optimisation |
   2019-04-02T13:55:43: +=================================+
   2019-04-02T13:55:44:  >Iteration   0: -logL=57895.557, Lambda=1.0e-03
   2019-04-02T13:55:46:  >Iteration   1: -logL=57893.741, Lambda=1.0e-03, delta=1.815, step=1.0e+00, max(|grad|)=5.091294 [Index:7]
   2019-04-02T13:55:47:  >Iteration   2: -logL=57893.740, Lambda=1.0e-04, delta=0.001, step=1.0e+00, max(|grad|)=0.076392 [Index:7]
   2019-04-02T13:55:48:
   2019-04-02T13:55:48: +=========================================+
   2019-04-02T13:55:48: | Maximum likelihood optimisation results |
   2019-04-02T13:55:48: +=========================================+
   2019-04-02T13:55:48: === GOptimizerLM ===
   2019-04-02T13:55:48:  Optimized function value ..: 57893.740
   2019-04-02T13:55:48:  Absolute precision ........: 0.005
   2019-04-02T13:55:48:  Acceptable value decrease .: 2
   2019-04-02T13:55:48:  Optimization status .......: converged
   2019-04-02T13:55:48:  Number of parameters ......: 10
   2019-04-02T13:55:48:  Number of free parameters .: 4
   2019-04-02T13:55:48:  Number of iterations ......: 2
   2019-04-02T13:55:48:  Lambda ....................: 1e-05
   2019-04-02T13:55:48:  Maximum log likelihood ....: -57893.740
   2019-04-02T13:55:48:  Observed events  (Nobs) ...: 19452.000
   2019-04-02T13:55:48:  Predicted events (Npred) ..: 19451.999 (Nobs - Npred = 0.00144846074545057)
   2019-04-02T13:55:48: === GModels ===
   2019-04-02T13:55:48:  Number of models ..........: 2
   2019-04-02T13:55:48:  Number of parameters ......: 10
   2019-04-02T13:55:48: === GModelSky ===
   2019-04-02T13:55:48:  Name ......................: Crab
   2019-04-02T13:55:48:  Instruments ...............: all
   2019-04-02T13:55:48:  Observation identifiers ...: all
   2019-04-02T13:55:48:  Model type ................: PointSource
   2019-04-02T13:55:48:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
   2019-04-02T13:55:48:  Number of parameters ......: 6
   2019-04-02T13:55:48:  Number of spatial par's ...: 2
   2019-04-02T13:55:48:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
   2019-04-02T13:55:48:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
   2019-04-02T13:55:48:  Number of spectral par's ..: 3
   2019-04-02T13:55:48:   Prefactor ................: 5.79620732342312e-16 +/- 1.01207859506902e-17 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
   2019-04-02T13:55:48:   Index ....................: -2.48767829376908 +/- 0.0150183262185825 [-0,-5]  (free,scale=-1,gradient)
   2019-04-02T13:55:48:   PivotEnergy ..............: 300000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
   2019-04-02T13:55:48:  Number of temporal par's ..: 1
   2019-04-02T13:55:48:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2019-04-02T13:55:48:  Number of scale par's .....: 0
   2019-04-02T13:55:48: === GCTAModelCubeBackground ===
   2019-04-02T13:55:48:  Name ......................: BackgroundModel
   2019-04-02T13:55:48:  Instruments ...............: CTA, HESS, MAGIC, VERITAS
   2019-04-02T13:55:48:  Observation identifiers ...: all
   2019-04-02T13:55:48:  Model type ................: "PowerLaw" * "Constant"
   2019-04-02T13:55:48:  Number of parameters ......: 4
   2019-04-02T13:55:48:  Number of spectral par's ..: 3
   2019-04-02T13:55:48:   Prefactor ................: 0.997287677101311 +/- 0.0156090791023617 [0.01,100] ph/cm2/s/MeV (free,scale=1,gradient)
   2019-04-02T13:55:48:   Index ....................: -0.00905917021700502 +/- 0.00939297053758322 [-5,5]  (free,scale=1,gradient)
   2019-04-02T13:55:48:   PivotEnergy ..............: 1000000 MeV (fixed,scale=1000000,gradient)
   2019-04-02T13:55:48:  Number of temporal par's ..: 1
   2019-04-02T13:55:48:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

The maximum likelihood optimizer required 2 iterations to converge. This
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
   (background) model or fitting no source. As a rule of thumb, the square
   root of the Test Statistic value gives the source detection significance
   in Gaussian sigmas, although the actual conversion depends somewhat on
   the formulation of the statistical problem and the number of
   degrees of freedom associated with the source.

   To instruct :ref:`ctlike` to compute the Test Statistic value for a
   given source you need to add the attribute ``tscalc="1"`` to the XML
   file:

   .. code-block:: xml

      <source name="Crab" type="PointSource" tscalc="1">

   :ref:`ctlike` will then compute the Test Statistic value for that
   source and dump the result in the log file:

   .. code-block:: none

      2019-04-02T13:58:58:  Name ......................: Crab
      2019-04-02T13:58:58:  Instruments ...............: all
      2019-04-02T13:58:58:  Test Statistic ............: 21529.9862965408

   The Test Statistic value will also be added as new attribute
   ``ts`` to the XML result file:

   .. code-block:: xml

      <source name="Crab" type="PointSource" ts="21529.986" tscalc="1">
