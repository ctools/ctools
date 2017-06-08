.. _1dc_first_unbinned:

Fitting the model components directly to the events
---------------------------------------------------

  .. admonition:: What you will learn

     You will learn how to analyse the event data using an **unbinned maximum
     likelihood analysis**.

     Use an unbinned maximum likelihood analysis if you have doubts about the
     impact of the selected binning on your analysis and if the number of
     observations is not too large (typically a few tens of 30 minute
     observations). Or if you simply to not worry about computation time.

Instead of binning the events into a counts cube and computing the effective
:ref:`instrument response functions <glossary_irf>`
you can directly fit the model to the events using an
unbinned maximum likelihood analysis. The advantage of the unbinned analysis
is that the event time remains an explicit quantity of the analysis and that
the different observations together with their
:ref:`instrument response functions <glossary_irf>`
are kept separate. In addition, no binning is performed and the analysis takes
full advantage of the information carried by each individual event. The
drawback is that the computing time scales linearly with the number of events
and the number of observations, becoming evetually prohibitive for the analysis
of large volumes of data (but for the analysis of small junks of data the
unbinned analysis is in fact faster than the stacked analysis).

To perform an unbinned model fitting you provide the
:ref:`observation definition file <glossary_obsdef>`
referencing the event lists as input to the :ref:`ctlike` tool:

.. code-block:: bash

  $ ctlike
  Input event list, counts cube or observation definition XML file [events.fits] obs_selected.xml
  Input model definition XML file [$CTOOLS/share/models/crab.xml] models_cutoff.xml
  Output model definition XML file [crab_results.xml] unbinned_results_cutoff.xml

The tool will take a few minutes (on Mac OS X) to perform the model fitting,
and will write the results into an updated
:ref:`model definition file <glossary_moddef>`
containing the fitted model parameters and their statistical uncertainties.
You may inspect the log file ``ctlike.log`` to verify that the model fit
converged properly, as illustrated in the example below:

.. code-block:: bash

   2017-06-08T09:56:42: +=================================+
   2017-06-08T09:56:42: | Maximum likelihood optimisation |
   2017-06-08T09:56:42: +=================================+
   2017-06-08T09:57:14:  >Iteration   0: -logL=22938810.009, Lambda=1.0e-03
   2017-06-08T09:57:43:  >Iteration   1: -logL=22914874.739, Lambda=1.0e-03, delta=23935.270, step=1.0e+00, max(|grad|)=52907.253765 [Index:14]
   2017-06-08T09:58:11:  >Iteration   2: -logL=22911216.020, Lambda=1.0e-04, delta=3658.719, step=1.0e+00, max(|grad|)=-13199.549393 [DEC:1]
   2017-06-08T09:58:41:  >Iteration   3: -logL=22909769.414, Lambda=1.0e-05, delta=1446.606, step=1.0e+00, max(|grad|)=-11987.289647 [DEC:1]
   2017-06-08T09:59:10:  >Iteration   4: -logL=22909405.587, Lambda=1.0e-06, delta=363.826, step=1.0e+00, max(|grad|)=-8985.645939 [DEC:1]
   2017-06-08T09:59:39:  >Iteration   5: -logL=22909337.513, Lambda=1.0e-07, delta=68.074, step=1.0e+00, max(|grad|)=-5891.649861 [DEC:1]
   2017-06-08T10:00:07:  >Iteration   6: -logL=22909309.777, Lambda=1.0e-08, delta=27.736, step=1.0e+00, max(|grad|)=-3762.571320 [DEC:1]
   2017-06-08T10:00:36:  >Iteration   7: -logL=22909295.909, Lambda=1.0e-09, delta=13.868, step=1.0e+00, max(|grad|)=-2458.462512 [RA:7]
   2017-06-08T10:01:05:  >Iteration   8: -logL=22909288.709, Lambda=1.0e-10, delta=7.200, step=1.0e+00, max(|grad|)=-1835.957201 [RA:7]
   2017-06-08T10:01:34:  >Iteration   9: -logL=22909284.904, Lambda=1.0e-11, delta=3.805, step=1.0e+00, max(|grad|)=-1368.259068 [RA:7]
   2017-06-08T10:02:04:  >Iteration  10: -logL=22909282.863, Lambda=1.0e-12, delta=2.041, step=1.0e+00, max(|grad|)=-1022.653174 [RA:7]
   2017-06-08T10:02:35:  >Iteration  11: -logL=22909281.751, Lambda=1.0e-13, delta=1.112, step=1.0e+00, max(|grad|)=-765.437910 [RA:7]
   2017-06-08T10:03:05:  >Iteration  12: -logL=22909281.140, Lambda=1.0e-14, delta=0.611, step=1.0e+00, max(|grad|)=-572.313029 [RA:7]
   2017-06-08T10:03:33:  >Iteration  13: -logL=22909280.802, Lambda=1.0e-15, delta=0.338, step=1.0e+00, max(|grad|)=-430.018412 [RA:7]
   2017-06-08T10:04:02:  >Iteration  14: -logL=22909280.613, Lambda=1.0e-16, delta=0.189, step=1.0e+00, max(|grad|)=-321.908506 [RA:7]
   2017-06-08T10:04:30:  >Iteration  15: -logL=22909280.507, Lambda=1.0e-17, delta=0.106, step=1.0e+00, max(|grad|)=-240.266994 [RA:7]
   2017-06-08T10:04:58:  >Iteration  16: -logL=22909280.448, Lambda=1.0e-18, delta=0.059, step=1.0e+00, max(|grad|)=-180.004153 [RA:7]
   2017-06-08T10:05:27:  >Iteration  17: -logL=22909280.414, Lambda=1.0e-19, delta=0.033, step=1.0e+00, max(|grad|)=-135.633545 [RA:7]
   2017-06-08T10:05:55:  >Iteration  18: -logL=22909280.395, Lambda=1.0e-20, delta=0.019, step=1.0e+00, max(|grad|)=-102.119554 [RA:7]
   2017-06-08T10:06:24:  >Iteration  19: -logL=22909280.385, Lambda=1.0e-21, delta=0.011, step=1.0e+00, max(|grad|)=-76.657216 [RA:7]
   2017-06-08T10:06:52:  >Iteration  20: -logL=22909280.379, Lambda=1.0e-22, delta=0.006, step=1.0e+00, max(|grad|)=-57.439799 [RA:7]
   2017-06-08T10:07:20:  >Iteration  21: -logL=22909280.375, Lambda=1.0e-23, delta=0.003, step=1.0e+00, max(|grad|)=-43.064181 [RA:7]
   2017-06-08T10:07:49:
   2017-06-08T10:07:49: +=========================================+
   2017-06-08T10:07:49: | Maximum likelihood optimisation results |
   2017-06-08T10:07:49: +=========================================+
   2017-06-08T10:07:49: === GOptimizerLM ===
   2017-06-08T10:07:49:  Optimized function value ..: 22909280.375
   2017-06-08T10:07:49:  Absolute precision ........: 0.005
   2017-06-08T10:07:49:  Acceptable value decrease .: 2
   2017-06-08T10:07:49:  Optimization status .......: converged
   2017-06-08T10:07:49:  Number of parameters ......: 17
   2017-06-08T10:07:49:  Number of free parameters .: 11
   2017-06-08T10:07:49:  Number of iterations ......: 21
   2017-06-08T10:07:49:  Lambda ....................: 1e-24
   2017-06-08T10:07:49:  Maximum log likelihood ....: -22909280.375
   2017-06-08T10:07:49:  Observed events  (Nobs) ...: 3349183.000
   2017-06-08T10:07:49:  Predicted events (Npred) ..: 3349177.996 (Nobs - Npred = 5.00393125554547)

Also the :ref:`ctbutterfly` tool and the :ref:`csspec` and :ref:`csresmap`
scripts can be run into unbinned mode by providing the
:ref:`observation definition file <glossary_obsdef>`
instead of the counts cube on input.
For illustration, the butterfly diagrams and spectra obtained for ``Src001``
and ``Src002`` using an unbinned maximum likelihood analysis are shown
below.

.. figure:: first_spectrum_cutoff_unbinned.png
   :width: 600px
   :align: center

   *Butterfly diagrams determined with ctbutterfly and spectral points obtained with csspec using an unbinned analysis*
