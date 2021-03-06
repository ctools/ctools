.. _start_unbinned:

Doing an unbinned analysis
--------------------------

  .. admonition:: What you will learn

     You will learn how to **adjust a parametrised model to the events without
     binning the data**.

     Gamma-ray events are rare, hence the counts cubes generated by
     :ref:`ctbin` may be sparsly populated, having many empty pixels, in
     particular at high energies. In that case it may be worth to
     analyse the events directly with an unbinned maximum likelihood
     analysis.

     An unbinned analysis is generally preferred over a binned analysis for
     short observation times (a few tens of hours) or if you want to assure
     that the analysis results are not biased by the selected binning.

An alternative analysis technique consists of working directly on the event
list without binning the events in a counts cube.
You do this with the :ref:`ctlike` tool as follows:

.. code-block:: bash

   $ ctlike
   Input event list, counts cube or observation definition XML file [cntcube.fits] selected_events.fits
   Calibration database [prod2]
   Instrument response function [South_0.5h]
   Input model definition XML file [models.xml] $CTOOLS/share/models/crab.xml
   Output model definition XML file [crab_results.xml]

You will recognise that :ref:`ctlike` runs much faster in unbinned mode
compared to binned mode.
This is understandable as the selected event list contains
only 20877 events, while the binned counts cube you used before had
200 x 200 x 20 = 800000 bins. As unbinned maximum likelihood fitting loops
over the events (while binned maximum likelihood loops over the bins),
there are much less operations to perform in unbinned than in binned mode
(there is some additional overhead in unbinned mode that comes from
integrating the models over the region of interest, yet this is negligible
compared to the operations needed when looping over all pixels). So as long
as you work with small event lists, unbinned mode is faster (this
typically holds up to few tens of hours of observing time).
Unbinned :ref:`ctlike` should also be more precise as no binning is performed,
hence there is no loss of information due to histogramming.

Below you see the corresponding output from the ``ctlike.log`` file. The fitted
parameters are essentially identical to the ones found in binned mode.
The slight difference with respect to the binned analysis may be explained
by the different event sample that were used for the analysis: while
binned likelihood works on rectangular counts cubes, unbinned likelihood works
on circular event selection regions. It is thus not possible to select exactly
the same events for both analyses.

.. code-block:: none

   2019-04-02T14:15:10: +=================================+
   2019-04-02T14:15:10: | Maximum likelihood optimisation |
   2019-04-02T14:15:10: +=================================+
   2019-04-02T14:15:10:  >Iteration   0: -logL=143782.837, Lambda=1.0e-03
   2019-04-02T14:15:10:  >Iteration   1: -logL=143779.346, Lambda=1.0e-03, delta=3.491, step=1.0e+00, max(|grad|)=5.346881 [Index:7]
   2019-04-02T14:15:10:  >Iteration   2: -logL=143779.343, Lambda=1.0e-04, delta=0.003, step=1.0e+00, max(|grad|)=-0.055672 [Index:3]
   2019-04-02T14:15:10:
   2019-04-02T14:15:10: +=========================================+
   2019-04-02T14:15:10: | Maximum likelihood optimisation results |
   2019-04-02T14:15:10: +=========================================+
   2019-04-02T14:15:10: === GOptimizerLM ===
   2019-04-02T14:15:10:  Optimized function value ..: 143779.343
   2019-04-02T14:15:10:  Absolute precision ........: 0.005
   2019-04-02T14:15:10:  Acceptable value decrease .: 2
   2019-04-02T14:15:10:  Optimization status .......: converged
   2019-04-02T14:15:10:  Number of parameters ......: 10
   2019-04-02T14:15:10:  Number of free parameters .: 4
   2019-04-02T14:15:10:  Number of iterations ......: 2
   2019-04-02T14:15:10:  Lambda ....................: 1e-05
   2019-04-02T14:15:10:  Maximum log likelihood ....: -143779.343
   2019-04-02T14:15:10:  Observed events  (Nobs) ...: 22708.000
   2019-04-02T14:15:10:  Predicted events (Npred) ..: 22707.995 (Nobs - Npred = 0.00519125738355797)
   2019-04-02T14:15:10: === GModels ===
   2019-04-02T14:15:10:  Number of models ..........: 2
   2019-04-02T14:15:10:  Number of parameters ......: 10
   2019-04-02T14:15:10: === GModelSky ===
   2019-04-02T14:15:10:  Name ......................: Crab
   2019-04-02T14:15:10:  Instruments ...............: all
   2019-04-02T14:15:10:  Observation identifiers ...: all
   2019-04-02T14:15:10:  Model type ................: PointSource
   2019-04-02T14:15:10:  Model components ..........: "PointSource" * "PowerLaw" * "Constant"
   2019-04-02T14:15:10:  Number of parameters ......: 6
   2019-04-02T14:15:10:  Number of spatial par's ...: 2
   2019-04-02T14:15:10:   RA .......................: 83.6331 [-360,360] deg (fixed,scale=1)
   2019-04-02T14:15:10:   DEC ......................: 22.0145 [-90,90] deg (fixed,scale=1)
   2019-04-02T14:15:10:  Number of spectral par's ..: 3
   2019-04-02T14:15:10:   Prefactor ................: 5.88338676901458e-16 +/- 1.02452856089807e-17 [1e-23,1e-13] ph/cm2/s/MeV (free,scale=1e-16,gradient)
   2019-04-02T14:15:10:   Index ....................: -2.49375950219757 +/- 0.0149889370322137 [-0,-5]  (free,scale=-1,gradient)
   2019-04-02T14:15:10:   PivotEnergy ..............: 300000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
   2019-04-02T14:15:10:  Number of temporal par's ..: 1
   2019-04-02T14:15:10:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
   2019-04-02T14:15:10:  Number of scale par's .....: 0
   2019-04-02T14:15:10: === GCTAModelIrfBackground ===
   2019-04-02T14:15:10:  Name ......................: CTABackgroundModel
   2019-04-02T14:15:10:  Instruments ...............: CTA
   2019-04-02T14:15:10:  Observation identifiers ...: all
   2019-04-02T14:15:10:  Model type ................: "PowerLaw" * "Constant"
   2019-04-02T14:15:10:  Number of parameters ......: 4
   2019-04-02T14:15:10:  Number of spectral par's ..: 3
   2019-04-02T14:15:10:   Prefactor ................: 1.0018169793538 +/- 0.0133053833141539 [0.001,1000] ph/cm2/s/MeV (free,scale=1,gradient)
   2019-04-02T14:15:10:   Index ....................: -0.00708154249642314 +/- 0.00805278228449961 [-5,5]  (free,scale=1,gradient)
   2019-04-02T14:15:10:   PivotEnergy ..............: 1000000 [10000,1000000000] MeV (fixed,scale=1000000,gradient)
   2019-04-02T14:15:10:  Number of temporal par's ..: 1
   2019-04-02T14:15:10:   Normalization ............: 1 (relative value) (fixed,scale=1,gradient)

.. note::
   Many tools or scripts can also be used in unbinned mode, including
   :ref:`csresmap`, :ref:`ctbutterfly` and :ref:`csspec` that were used
   earlier. It is sufficient to replace the input counts cube by an event
   list to activate unbinned mode for these tools.
