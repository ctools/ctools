.. _1dc_first_stacked:

Performing a stacked analysis
-----------------------------

In the previous example you have fit the model to the data using an unbinned
maximum likelihood method, which directly operates on the event lists without
any binning of the data. This is the most precise method to analysis the CTA
data, yet as you may have recognised, it is a time consuming method, in
particular if you want to analyse a large number of events.

An alternative analysis method is the stacked analysis, where the data from
all observations are binned into a single
:ref:`counts cube <glossary_countscube>`,
spanned by
Right Ascension (or Galactic longitude),
Declination (or Galactic latitude), and energy.
For observation durations exceeding a few 10 hours, a stacked analysis is
generally faster than an unbinned analysis, but any temporal information
is lost.

To perform a stacked analyses, you start with creating a
:ref:`counts cube <glossary_countscube>`
by typing

.. code-block:: bash

   $ ctbin
   Input event list or observation definition XML file [events.fits] obs.xml
   First coordinate of image center in degrees (RA or galactic l) (0-360) [0] 0.0
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [0] 0.0
   Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [GAL]
   Image scale (in degrees/pixel) [0.02]
   Size of the X axis in pixels [400]
   Size of the Y axis in pixels [400]
   Algorithm for defining energy bins (FILE|LIN|LOG) [LOG]
   Start value for first energy bin in TeV [0.1]
   Stop value for last energy bin in TeV [100]
   Number of energy bins (1-200) [20]
   Output counts cube file [cntcube.fits]

This creates a 3-dimensional counts cube in Galactic coordinates, centred on
the Galactic centre and 8 degrees x 8 degrees wide, with 20 logarithmically
spaced energy bins between 100 GeV and 100 TeV.

As next step you have to compute the effective response function for the
counts cube. The effective response function is composed of an exposure
cube, a point spread function cube, and a background cube. You generate
the cubes by typing

.. code-block:: bash

   $ ctexpcube
   Input event list or observation definition XML file [NONE] obs.xml
   Input counts cube file to extract exposure cube definition [NONE] cntcube.fits
   Output exposure cube file [expcube.fits]

.. code-block:: bash

   $ ctpsfcube
   Input event list or observation definition XML file [NONE] obs.xml
   Input counts cube file to extract PSF cube definition [NONE]
   First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 0.0
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 0.0
   Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] GAL
   Image scale (in degrees/pixel) [1.0]
   Size of the X axis in pixels [10]
   Size of the Y axis in pixels [10]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins [20]
   Output PSF cube file [psfcube.fits]

.. code-block:: bash

   $ ctbkgcube
   Input event list or observation definition XML file [NONE] obs.xml
   Input counts cube file to extract background cube definition [NONE] cntcube.fits
   Input model definition XML file [NONE] models.xml
   Output background cube file [bkgcube.fits]
   Output model definition XML file [NONE] stacked_models.xml

.. note::
   For convenience, the dimensions of the exposure and background cubes were
   taken identical to the dimensions of the counts cube. This is however not
   a requirement, and each of the cubes may have a different dimension and
   size. For the point spread function cube a coarser spatial binning was
   used to keep the size of the cube at a manageable level. The point spread
   function cube varies in fact only slowly over the field of view, and a
   coarse spatial binning is sufficient to capture that variability.

With this prepatory work finished, you can now perform a binned maximum
likelihood using :ref:`ctlike`. Instead of the
:ref:`Observation Definition File <glossary_obsdef>`
specified for the unbinned analysis, you now need to specify the
:ref:`counts cube <glossary_countscube>`
on input, and :ref:`ctlike` will then automatically query for the response
cubes:

.. code-block:: bash

   $ ctlike
   Input event list, counts cube or observation definition XML file [obs.xml] cntcube.fits
   Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits
   Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits
   Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits
   Input model definition XML file [models.xml] stacked_models.xml
   Output model definition XML file [results.xml]

You may recognise that :ref:`ctlike` now runs significantly faster.
An inspection of the log file ``ctlike.log`` demonstrates that the model fit
converged properly:

.. code-block:: bash

   2017-02-02T16:39:58: +=================================+
   2017-02-02T16:39:58: | Maximum likelihood optimisation |
   2017-02-02T16:39:58: +=================================+
   2017-02-02T16:40:14:  >Iteration   0: -logL=-1011176.684, Lambda=1.0e-03
   2017-02-02T16:40:30:  >Iteration   1: -logL=-1018314.686, Lambda=1.0e-03, delta=7138.002, max(|grad|)=-22791.121954 [DEC:7]
   2017-02-02T16:40:46:  >Iteration   2: -logL=-1020732.797, Lambda=1.0e-04, delta=2418.112, max(|grad|)=-16741.951327 [DEC:7]
   2017-02-02T16:41:02:  >Iteration   3: -logL=-1022040.008, Lambda=1.0e-05, delta=1307.210, max(|grad|)=-10516.727948 [DEC:7]
   2017-02-02T16:41:19:  >Iteration   4: -logL=-1022455.353, Lambda=1.0e-06, delta=415.345, max(|grad|)=-4793.902211 [DEC:7]
   2017-02-02T16:41:34:  >Iteration   5: -logL=-1022499.179, Lambda=1.0e-07, delta=43.827, max(|grad|)=-1120.934327 [DEC:7]
   2017-02-02T16:41:50:  >Iteration   6: -logL=-1022499.914, Lambda=1.0e-08, delta=0.735, max(|grad|)=-84.494765 [RA:6]
   2017-02-02T16:42:06:  >Iteration   7: -logL=-1022499.918, Lambda=1.0e-09, delta=0.003, max(|grad|)=-8.207385 [RA:6]
   2017-02-02T16:42:22:
   2017-02-02T16:42:22: +=========================================+
   2017-02-02T16:42:22: | Maximum likelihood optimisation results |
   2017-02-02T16:42:22: +=========================================+
   2017-02-02T16:42:22: === GOptimizerLM ===
   2017-02-02T16:42:22:  Optimized function value ..: -1022499.918
   2017-02-02T16:42:22:  Absolute precision ........: 0.005
   2017-02-02T16:42:22:  Acceptable value decrease .: 2
   2017-02-02T16:42:22:  Optimization status .......: converged
   2017-02-02T16:42:22:  Number of parameters ......: 16
   2017-02-02T16:42:22:  Number of free parameters .: 10
   2017-02-02T16:42:22:  Number of iterations ......: 7
   2017-02-02T16:42:22:  Lambda ....................: 1e-10
   2017-02-02T16:42:22:  Maximum log likelihood ....: 1022499.918
   2017-02-02T16:42:22:  Observed events  (Nobs) ...: 3136493.000
   2017-02-02T16:42:22:  Predicted events (Npred) ..: 3136492.994 (Nobs - Npred = 0.00608432479202747)

Similar to the unbinned analysis you can use
:ref:`csresmap`
to compute a residual map for a counts cube by typing

.. code-block:: bash

   $ csresmap
   Input event list, counts cube, or observation definition XML file [obs.xml] cntcube.fits
   Input model cube file (generated with ctmodel) [NONE]
   Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits
   Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits
   Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits
   Input model definition XML file [results.xml]
   Residual map computation algorithm (SUB|SUBDIV|SUBDIVSQRT) [SUB]
   Output residual map file [resmap.fits]

which generates a residual map that is indistinguishable to the map shown for
the unbinned analysis (see below).

.. figure:: first_skymap_residual_stacked.png
   :width: 400px
   :align: center

   *Residual sky map after subtraction of the fitted model for a stacked analysis*

You can also use
:ref:`ctbutterfly`
to compute a butterfly diagram for a counts cube by typing

.. code-block:: bash

   $ ctbutterfly
   Input event list, counts cube or observation definition XML file [obs.xml] cntcube.fits
   Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits
   Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits
   Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits
   Source of interest [Src002] Src001
   Input model definition XML file [results.xml]
   Start value for first energy bin in TeV [0.1]
   Stop value for last energy bin in TeV [100.0]
   Output ASCII file [butterfly_src002.txt] butterfly_stacked_src001.txt

and you can use
:ref:`csspec`
to derive a source spectrum from a counts cube by typing

.. code-block:: bash

   $ csspec
   Input event list, counts cube, or observation definition XML file [obs.xml] cntcube.fits
   Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits
   Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits
   Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits
   Input model definition XML file [results.xml]
   Source name [Src002] Src001
   Binning algorithm (LIN|LOG|FILE) [LOG]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins (0=unbinned) [10]
   Output spectrum file [spectrum_src002.fits] spectrum_stacked_src001.fits

The results are again indistinguishable to the result for unbinned analysis:

.. figure:: first_spectrum_stacked.png
   :width: 600px
   :align: center

   *Butterfly diagrams determined with ctbutterfly and spectral points determined with csspec for Src001 (top) and Src002 (bottom) using a stacked analysis*

