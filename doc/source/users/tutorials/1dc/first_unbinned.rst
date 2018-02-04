.. _1dc_first_unbinned:

Fitting the model components directly to the events
---------------------------------------------------

  .. admonition:: What you will learn

     You will learn how to analyse the event data using an **unbinned maximum
     likelihood analysis**.

     Use an unbinned maximum likelihood analysis if you have doubts about the
     impact of the selected binning on your analysis and if the number of
     observations is not too large (typically a few tens of 30 minute
     observations). Or if you simply do not worry about computation time.

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
   Input event list, counts cube or observation definition XML file [cntcube.fits] obs_selected.xml
   Input model definition XML file [bkgcube_eplaw_iem.xml] models_eplaw.xml
   Output model definition XML file [results_stacked_eplaw_iem.xml] results_unbinned_eplaw.xml

The tool will take a few minutes (on Mac OS X) to perform the model fitting,
and will write the results into an updated
:ref:`model definition file <glossary_moddef>`
containing the fitted model parameters and their statistical uncertainties.
You may inspect the log file ``ctlike.log`` to verify that the model fit
converged properly, as illustrated in the example below:

.. code-block:: none

   2018-01-26T13:57:31: +=================================+
   2018-01-26T13:57:31: | Maximum likelihood optimisation |
   2018-01-26T13:57:31: +=================================+
   2018-01-26T13:57:59:  >Iteration   0: -logL=20951221.333, Lambda=1.0e-03
   2018-01-26T13:58:24:  >Iteration   1: -logL=20920374.352, Lambda=1.0e-03, delta=30846.981, step=1.0e+00, max(|grad|)=74106.669528 [Index:14]
   2018-01-26T13:58:49:  >Iteration   2: -logL=20914723.637, Lambda=1.0e-04, delta=5650.715, step=1.0e+00, max(|grad|)=6059.304173 [Index:14]
   2018-01-26T13:59:15:  >Iteration   3: -logL=20912552.166, Lambda=1.0e-05, delta=2171.471, step=1.0e+00, max(|grad|)=4361.007623 [RA:0]
   2018-01-26T13:59:41:  >Iteration   4: -logL=20911666.571, Lambda=1.0e-06, delta=885.595, step=1.0e+00, max(|grad|)=2602.512475 [RA:0]
   2018-01-26T14:00:07:  >Iteration   5: -logL=20911523.902, Lambda=1.0e-07, delta=142.669, step=1.0e+00, max(|grad|)=1076.703629 [RA:0]
   2018-01-26T14:00:31:  >Iteration   6: -logL=20911508.855, Lambda=1.0e-08, delta=15.047, step=1.0e+00, max(|grad|)=333.013402 [RA:0]
   2018-01-26T14:00:56:  >Iteration   7: -logL=20911507.544, Lambda=1.0e-09, delta=1.311, step=1.0e+00, max(|grad|)=93.512904 [RA:0]
   2018-01-26T14:01:22:  >Iteration   8: -logL=20911507.529, Lambda=1.0e-10, delta=0.015, step=1.0e+00, max(|grad|)=24.977311 [RA:0]
   2018-01-26T14:01:46:  >Iteration   9: -logL=20911507.529, Lambda=1.0e-11, delta=0.000, step=1.0e+00, max(|grad|)=6.442211 [RA:0]
   2018-01-26T14:02:12:
   2018-01-26T14:02:12: +=========================================+
   2018-01-26T14:02:12: | Maximum likelihood optimisation results |
   2018-01-26T14:02:12: +=========================================+
   2018-01-26T14:02:12: === GOptimizerLM ===
   2018-01-26T14:02:12:  Optimized function value ..: 20911507.529
   2018-01-26T14:02:12:  Absolute precision ........: 0.005
   2018-01-26T14:02:12:  Acceptable value decrease .: 2
   2018-01-26T14:02:12:  Optimization status .......: converged
   2018-01-26T14:02:12:  Number of parameters ......: 17
   2018-01-26T14:02:12:  Number of free parameters .: 11
   2018-01-26T14:02:12:  Number of iterations ......: 9
   2018-01-26T14:02:12:  Lambda ....................: 1e-12
   2018-01-26T14:02:12:  Maximum log likelihood ....: -20911507.529
   2018-01-26T14:02:12:  Observed events  (Nobs) ...: 3084595.000
   2018-01-26T14:02:12:  Predicted events (Npred) ..: 3084591.999 (Nobs - Npred = 3.00117048900574)

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
