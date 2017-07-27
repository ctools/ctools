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

   2017-07-27T16:28:26: +=================================+
   2017-07-27T16:28:26: | Maximum likelihood optimisation |
   2017-07-27T16:28:26: +=================================+
   2017-07-27T16:28:54:  >Iteration   0: -logL=20949416.217, Lambda=1.0e-03
   2017-07-27T16:29:21:  >Iteration   1: -logL=20918674.227, Lambda=1.0e-03, delta=30741.990, step=1.0e+00, max(|grad|)=74249.049334 [Index:14]
   2017-07-27T16:29:49:  >Iteration   2: -logL=20913114.563, Lambda=1.0e-04, delta=5559.664, step=1.0e+00, max(|grad|)=6728.499004 [RA:0]
   2017-07-27T16:30:16:  >Iteration   3: -logL=20910994.152, Lambda=1.0e-05, delta=2120.411, step=1.0e+00, max(|grad|)=4949.846735 [RA:0]
   2017-07-27T16:30:44:  >Iteration   4: -logL=20910115.048, Lambda=1.0e-06, delta=879.104, step=1.0e+00, max(|grad|)=2889.871293 [RA:0]
   2017-07-27T16:31:12:  >Iteration   5: -logL=20909966.036, Lambda=1.0e-07, delta=149.012, step=1.0e+00, max(|grad|)=1190.543059 [RA:0]
   2017-07-27T16:31:39:  >Iteration   6: -logL=20909948.012, Lambda=1.0e-08, delta=18.024, step=1.0e+00, max(|grad|)=550.745751 [DEC:1]
   2017-07-27T16:32:06:  >Iteration   7: -logL=20909946.600, Lambda=1.0e-09, delta=1.412, step=1.0e+00, max(|grad|)=208.667347 [DEC:1]
   2017-07-27T16:32:33:  >Iteration   8: -logL=20909946.575, Lambda=1.0e-10, delta=0.025, step=1.0e+00, max(|grad|)=51.069713 [DEC:1]
   2017-07-27T16:33:00:  >Iteration   9: -logL=20909946.573, Lambda=1.0e-11, delta=0.002, step=1.0e+00, max(|grad|)=16.914402 [DEC:1]
   2017-07-27T16:33:26:
   2017-07-27T16:33:26: +=========================================+
   2017-07-27T16:33:26: | Maximum likelihood optimisation results |
   2017-07-27T16:33:26: +=========================================+
   2017-07-27T16:33:26: === GOptimizerLM ===
   2017-07-27T16:33:26:  Optimized function value ..: 20909946.573
   2017-07-27T16:33:26:  Absolute precision ........: 0.005
   2017-07-27T16:33:26:  Acceptable value decrease .: 2
   2017-07-27T16:33:26:  Optimization status .......: converged
   2017-07-27T16:33:26:  Number of parameters ......: 17
   2017-07-27T16:33:26:  Number of free parameters .: 11
   2017-07-27T16:33:26:  Number of iterations ......: 9
   2017-07-27T16:33:26:  Lambda ....................: 1e-12
   2017-07-27T16:33:26:  Maximum log likelihood ....: -20909946.573
   2017-07-27T16:33:26:  Observed events  (Nobs) ...: 3084639.000
   2017-07-27T16:33:26:  Predicted events (Npred) ..: 3084635.972 (Nobs - Npred = 3.0282306545414)

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
