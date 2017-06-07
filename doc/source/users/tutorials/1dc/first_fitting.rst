.. _1dc_first_fitting:

Fitting the model components to the counts cube
-----------------------------------------------

  .. admonition:: What you will learn

     You will learn how to **fit a parametric model** to the counts cube using
     a maximum likelihood algorithm.

     You will also learn how to display the fit results in form of
     `ds9 <http://ds9.si.edu>`_
     region files, butterfly diagrams and spectral points.

     Please note that display tools are **not** part of ctools, yet some
     scripts for result display that use the ``matplotlib`` Python module can be
     found in the
     ``$CTOOLS/share/examples/python``
     folder (see :ref:`how to display results<1dc_howto_display>`).

Now you are ready to fit the model to the counts cube and to determine its
maximum likelihood parameters.

You do this with the :ref:`ctlike` tool that adjusts all parameters in the
:ref:`model definition file <glossary_moddef>`
that have the attribute ``free`` set to ``"1"``.
In the current example, the free model parameters are the positions and spectral
parameters of the two point sources and the spectral normalisation of the
background component.
You run the :ref:`ctlike` tool as follows:

.. code-block:: bash

   $ ctlike
   Input event list, counts cube or observation definition XML file [events.fits] cntcube.fits
   Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits
   Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits
   Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits
   Input model definition XML file [$CTOOLS/share/models/crab.xml] stacked_models.xml
   Output model definition XML file [crab_results.xml] stacked_results.xml

The tool will take a few minutes (on Mac OS X) to perform the model fitting,
and will write the results into an updated
:ref:`model definition file <glossary_moddef>`
containing the fitted model parameters and their statistical uncertainties.
You may inspect the log file ``ctlike.log`` to verify that the model fit
converged properly, as illustrated in the example below:

.. code-block:: bash

   2017-06-06T23:54:05: +=================================+
   2017-06-06T23:54:05: | Maximum likelihood optimisation |
   2017-06-06T23:54:05: +=================================+
   2017-06-06T23:54:18:  >Iteration   0: -logL=242350.348, Lambda=1.0e-03
   2017-06-06T23:54:31:  >Iteration   1: -logL=222183.040, Lambda=1.0e-03, delta=20167.308, step=1.0e+00, max(|grad|)=50064.145481 [Index:13]
   2017-06-06T23:54:44:  >Iteration   2: -logL=218671.485, Lambda=1.0e-04, delta=3511.555, step=1.0e+00, max(|grad|)=-11302.678938 [RA:6]
   2017-06-06T23:54:57:  >Iteration   3: -logL=217673.527, Lambda=1.0e-05, delta=997.958, step=1.0e+00, max(|grad|)=-9593.644236 [RA:6]
   2017-06-06T23:55:10:  >Iteration   4: -logL=217505.751, Lambda=1.0e-06, delta=167.776, step=1.0e+00, max(|grad|)=-6625.379605 [RA:6]
   2017-06-06T23:55:23:  >Iteration   5: -logL=217465.972, Lambda=1.0e-07, delta=39.779, step=1.0e+00, max(|grad|)=-4315.893743 [RA:6]
   2017-06-06T23:55:37:  >Iteration   6: -logL=217446.648, Lambda=1.0e-08, delta=19.324, step=1.0e+00, max(|grad|)=3158.486974 [RA:0]
   2017-06-06T23:55:50:  >Iteration   7: -logL=217437.175, Lambda=1.0e-09, delta=9.473, step=1.0e+00, max(|grad|)=2234.551298 [RA:0]
   2017-06-06T23:56:03:  >Iteration   8: -logL=217432.730, Lambda=1.0e-10, delta=4.445, step=1.0e+00, max(|grad|)=1578.966812 [RA:0]
   2017-06-06T23:56:17:  >Iteration   9: -logL=217430.542, Lambda=1.0e-11, delta=2.187, step=1.0e+00, max(|grad|)=1164.679821 [RA:0]
   2017-06-06T23:56:31:  >Iteration  10: -logL=217429.374, Lambda=1.0e-12, delta=1.168, step=1.0e+00, max(|grad|)=816.726019 [RA:0]
   2017-06-06T23:56:44:  >Iteration  11: -logL=217428.835, Lambda=1.0e-13, delta=0.539, step=1.0e+00, max(|grad|)=543.038385 [RA:0]
   2017-06-06T23:56:58:  >Iteration  12: -logL=217428.635, Lambda=1.0e-14, delta=0.201, step=1.0e+00, max(|grad|)=290.465021 [RA:0]
   2017-06-06T23:57:12:  >Iteration  13: -logL=217428.574, Lambda=1.0e-15, delta=0.061, step=1.0e+00, max(|grad|)=223.932827 [RA:0]
   2017-06-06T23:57:25:  >Iteration  14: -logL=217428.528, Lambda=1.0e-16, delta=0.046, step=1.0e+00, max(|grad|)=228.672746 [RA:0]
   2017-06-06T23:57:38:  >Iteration  15: -logL=217428.480, Lambda=1.0e-17, delta=0.049, step=1.0e+00, max(|grad|)=179.805198 [RA:0]
   2017-06-06T23:57:52:  >Iteration  16: -logL=217428.459, Lambda=1.0e-18, delta=0.021, step=1.0e+00, max(|grad|)=103.278218 [RA:0]
   2017-06-06T23:58:05:  >Iteration  17: -logL=217428.452, Lambda=1.0e-19, delta=0.006, step=1.0e+00, max(|grad|)=62.896707 [RA:0]
   2017-06-06T23:58:18:  >Iteration  18: -logL=217428.448, Lambda=1.0e-20, delta=0.004, step=1.0e+00, max(|grad|)=41.622300 [RA:0]
   2017-06-06T23:58:31:
   2017-06-06T23:58:31: +=========================================+
   2017-06-06T23:58:31: | Maximum likelihood optimisation results |
   2017-06-06T23:58:31: +=========================================+
   2017-06-06T23:58:31: === GOptimizerLM ===
   2017-06-06T23:58:31:  Optimized function value ..: 217428.448
   2017-06-06T23:58:31:  Absolute precision ........: 0.005
   2017-06-06T23:58:31:  Acceptable value decrease .: 2
   2017-06-06T23:58:31:  Optimization status .......: converged
   2017-06-06T23:58:31:  Number of parameters ......: 16
   2017-06-06T23:58:31:  Number of free parameters .: 10
   2017-06-06T23:58:31:  Number of iterations ......: 18
   2017-06-06T23:58:31:  Lambda ....................: 1e-21
   2017-06-06T23:58:31:  Maximum log likelihood ....: -217428.448
   2017-06-06T23:58:31:  Observed events  (Nobs) ...: 2204717.000
   2017-06-06T23:58:31:  Predicted events (Npred) ..: 2204716.997 (Nobs - Npred = 0.00284357415512204)

You may also convert the fitted model positions into a `ds9 <http://ds9.si.edu>`_
region file using the :ref:`csmodelinfo` script so that you can overlay the
fit results over a sky map:

.. code-block:: bash

   $ csmodelinfo pnt_type=circle free_color=black show_labels=no
   Input model definition XML file [model.xml] stacked_results.xml
   Output DS9 region file [ds9.reg] positions.reg

The command line arguments ``pnt_type``, ``free_color`` and ``show_labels``
enable to fine tune the parameters in the `ds9 <http://ds9.si.edu>`_
region file. In this case, the positions are marked by black circles without
showing the source names.

The following image shows a zoom of the sky map that comprises both point
sources, with the initial source positions determined by :ref:`cssrcdetect`
as green crosses and the positions fitted by :ref:`ctlike` as black circles.
Obviously, the initial positions were already near the fitted positions,
which is required to assure the proper convergence of the fit.

.. figure:: first_skymap_fitted.png
   :width: 600px
   :align: center

   *Background subtracted sky map of the events recorded around the Galactic Centre during the Galactic Plane Survey with the fitted positions of the sources shown as black circles*

You can also convert the spectral parameters of the point sources into a
butterfly diagram for each source using the :ref:`ctbutterfly` tool.
The butterfly diagram shows the envelope of all spectral models that are
statistically compatible with the data.
You create the butterfly diagram for the first source using

.. code-block:: bash

  $ ctbutterfly
  Input event list, counts cube or observation definition XML file [events.fits] cntcube.fits
  Input exposure cube file (only needed for stacked analysis) [ctexpcube.fits] expcube.fits
  Input PSF cube file (only needed for stacked analysis) [psfcube.fits]
  Input background cube file (only needed for stacked analysis) [bkgcube.fits]
  Source of interest [Crab] Src001
  Input model definition XML file [$CTOOLS/share/models/crab.xml] stacked_results.xml
  Start value for first energy bin in TeV [0.1]
  Stop value for last energy bin in TeV [100.0]
  Output ASCII file [butterfly.txt] butterfly_src001.txt

and for the second source using

.. code-block:: bash

   $ ctbutterfly
   Input event list, counts cube or observation definition XML file [cntcube.fits]
   Input exposure cube file (only needed for stacked analysis) [expcube.fits]
   Input PSF cube file (only needed for stacked analysis) [psfcube.fits]
   Input background cube file (only needed for stacked analysis) [bkgcube.fits]
   Source of interest [Src001] Src002
   Input model definition XML file [stacked_results.xml]
   Start value for first energy bin in TeV [0.1]
   Stop value for last energy bin in TeV [100.0]
   Output ASCII file [butterfly_src001.txt] butterfly_src002.txt

The butterfly diagrams for both sources are displayed in the figure below.
The figure also shows spectral points for each source that were determined
using the :ref:`csspec` script.
You create the spectrum for the first source using

.. code-block:: bash

   $ csspec
   Input event list, counts cube, or observation definition XML file [events.fits] cntcube.fits
   Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits
   Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits
   Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits
   Input model definition XML file [$CTOOLS/share/models/crab.xml] stacked_results.xml
   Source name [Crab] Src001
   Binning algorithm (LIN|LOG|FILE) [LOG]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins (0=unbinned) [20] 10
   Output spectrum file [spectrum.fits] spectrum_src001.fits

and for the second source using

.. code-block:: bash

   $ csspec
   Input event list, counts cube, or observation definition XML file [cntcube.fits]
   Input exposure cube file (only needed for stacked analysis) [expcube.fits]
   Input PSF cube file (only needed for stacked analysis) [psfcube.fits]
   Input background cube file (only needed for stacked analysis) [bkgcube.fits]
   Input model definition XML file [stacked_results.xml]
   Source name [Src001] Src002
   Binning algorithm (LIN|LOG|FILE) [LOG]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins (0=unbinned) [10]
   Output spectrum file [spectrum_src001.fits] spectrum_src002.fits

The :ref:`csspec` script divided here the data into ten logarithmically
spaced energy bins and determined the source flux in each of the bins using
a maximum likelihood model fit.

.. figure:: first_spectrum_stacked.png
   :width: 600px
   :align: center

   *Butterfly diagrams determined with ctbutterfly and spectral points determined with csspec for Src001 (red) and Src002 (blue)*

Obviously, ``Src001`` has a spectral cut-off (red flux points) and hence is not
adequately described by a power law model.
Replacing the spectral model of ``Src001`` by an exponentially cut-off
power law improves the fit to the data, as illustrated by the figure below.

.. figure:: first_spectrum_cutoff_stacked.png
   :width: 600px
   :align: center

   *Butterfly diagrams determined with ctbutterfly for an exponentially cut-off power law for Src001 (red)*
