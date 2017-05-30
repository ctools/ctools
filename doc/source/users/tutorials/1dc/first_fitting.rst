.. _1dc_first_fitting:

Fitting the model components to the data
----------------------------------------

Now you are ready to fit the model to the data and to determine its maximum
likelihood parameters.

You do this with the :ref:`ctlike` tool that adjusts all parameters in the
:ref:`model definition file <glossary_moddef>`
that have the attribute ``free`` set to ``"1"``.
In the current example, the free model parameters are the positions and spectral
parameters of the two point sources and the spectral normalisation of the
background component.
You run the :ref:`ctlike` tool as follows:

.. code-block:: bash

   $ ctlike
   Input event list, counts cube or observation definition XML file [events.fits] obs_selected.xml
   Input model definition XML file [$CTOOLS/share/models/crab.xml] models.xml
   Output model definition XML file [crab_results.xml] results.xml

The tool will take a few minutes (on Mac OS X) to perform the model fitting,
and will write the results into an updated
:ref:`model definition file <glossary_moddef>`
containing the fitted model parameters and their statistical uncertainties.
You may inspect the log file ``ctlike.log`` to verify that the model fit
converged properly, as illustrated in the example below:

.. code-block:: bash

   2017-05-30T14:04:49: +=================================+
   2017-05-30T14:04:49: | Maximum likelihood optimisation |
   2017-05-30T14:04:49: +=================================+
   2017-05-30T14:05:14:  >Iteration   0: -logL=22930903.429, Lambda=1.0e-03
   2017-05-30T14:05:38:  >Iteration   1: -logL=22906772.079, Lambda=1.0e-03, delta=24131.350, step=1.0e+00, max(|grad|)=53111.535756 [Index:13]
   2017-05-30T14:06:00:  >Iteration   2: -logL=22903418.347, Lambda=1.0e-04, delta=3353.732, step=1.0e+00, max(|grad|)=-11638.756860 [RA:0]
   2017-05-30T14:06:23:  >Iteration   3: -logL=22902335.786, Lambda=1.0e-05, delta=1082.561, step=1.0e+00, max(|grad|)=-9993.636957 [RA:0]
   2017-05-30T14:06:46:  >Iteration   4: -logL=22902144.463, Lambda=1.0e-06, delta=191.323, step=1.0e+00, max(|grad|)=-7053.128713 [RA:0]
   2017-05-30T14:07:09:  >Iteration   5: -logL=22902098.908, Lambda=1.0e-07, delta=45.555, step=1.0e+00, max(|grad|)=4767.092710 [RA:6]
   2017-05-30T14:07:33:  >Iteration   6: -logL=22902077.092, Lambda=1.0e-08, delta=21.816, step=1.0e+00, max(|grad|)=3512.238312 [RA:6]
   2017-05-30T14:07:56:  >Iteration   7: -logL=22902066.319, Lambda=1.0e-09, delta=10.774, step=1.0e+00, max(|grad|)=2564.648439 [RA:6]
   2017-05-30T14:08:21:  >Iteration   8: -logL=22902060.969, Lambda=1.0e-10, delta=5.350, step=1.0e+00, max(|grad|)=1861.818754 [RA:6]
   2017-05-30T14:08:45:  >Iteration   9: -logL=22902058.299, Lambda=1.0e-11, delta=2.670, step=1.0e+00, max(|grad|)=1346.424825 [RA:6]
   2017-05-30T14:09:08:  >Iteration  10: -logL=22902056.961, Lambda=1.0e-12, delta=1.338, step=1.0e+00, max(|grad|)=971.176097 [RA:6]
   2017-05-30T14:09:31:  >Iteration  11: -logL=22902056.287, Lambda=1.0e-13, delta=0.674, step=1.0e+00, max(|grad|)=699.274062 [RA:6]
   2017-05-30T14:09:54:  >Iteration  12: -logL=22902055.947, Lambda=1.0e-14, delta=0.340, step=1.0e+00, max(|grad|)=502.893171 [RA:6]
   2017-05-30T14:10:18:  >Iteration  13: -logL=22902055.774, Lambda=1.0e-15, delta=0.172, step=1.0e+00, max(|grad|)=361.372646 [RA:6]
   2017-05-30T14:10:41:  >Iteration  14: -logL=22902055.687, Lambda=1.0e-16, delta=0.088, step=1.0e+00, max(|grad|)=259.519694 [RA:6]
   2017-05-30T14:11:04:  >Iteration  15: -logL=22902055.642, Lambda=1.0e-17, delta=0.045, step=1.0e+00, max(|grad|)=186.309625 [RA:6]
   2017-05-30T14:11:27:  >Iteration  16: -logL=22902055.620, Lambda=1.0e-18, delta=0.023, step=1.0e+00, max(|grad|)=133.719821 [RA:6]
   2017-05-30T14:11:50:  >Iteration  17: -logL=22902055.608, Lambda=1.0e-19, delta=0.012, step=1.0e+00, max(|grad|)=95.962554 [RA:6]
   2017-05-30T14:12:13:  >Iteration  18: -logL=22902055.602, Lambda=1.0e-20, delta=0.006, step=1.0e+00, max(|grad|)=68.863372 [RA:6]
   2017-05-30T14:12:35:  >Iteration  19: -logL=22902055.599, Lambda=1.0e-21, delta=0.003, step=1.0e+00, max(|grad|)=49.416455 [RA:6]
   2017-05-30T14:12:58:
   2017-05-30T14:12:58: +=========================================+
   2017-05-30T14:12:58: | Maximum likelihood optimisation results |
   2017-05-30T14:12:58: +=========================================+
   2017-05-30T14:12:58: === GOptimizerLM ===
   2017-05-30T14:12:58:  Optimized function value ..: 22902055.599
   2017-05-30T14:12:58:  Absolute precision ........: 0.005
   2017-05-30T14:12:58:  Acceptable value decrease .: 2
   2017-05-30T14:12:58:  Optimization status .......: converged
   2017-05-30T14:12:58:  Number of parameters ......: 16
   2017-05-30T14:12:58:  Number of free parameters .: 10
   2017-05-30T14:12:58:  Number of iterations ......: 19
   2017-05-30T14:12:58:  Lambda ....................: 1e-22
   2017-05-30T14:12:58:  Maximum log likelihood ....: -22902055.599
   2017-05-30T14:12:58:  Observed events  (Nobs) ...: 3348255.000
   2017-05-30T14:12:58:  Predicted events (Npred) ..: 3348253.996 (Nobs - Npred = 1.00355448853225)

You may also convert the fitted model positions into a `ds9 <http://ds9.si.edu>`_
region file using the :ref:`csmodelinfo` script so that you can overlay the
fit results over a sky map:

.. code-block:: bash

   $ csmodelinfo pnt_type=circle free_color=black show_labels=no
   Input model definition XML file [model.xml] results.xml
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
The butterfly diagram shows the envelope of all power laws that are
statistically compatible with the data.
You create the butterfly diagram for the first source using

.. code-block:: bash

   $ ctbutterfly
   Input event list, counts cube or observation definition XML file [events.fits] obs_selected.xml
   Source of interest [Crab] Src001
   Input model definition XML file [$CTOOLS/share/models/crab.xml] results.xml
   Start value for first energy bin in TeV [0.1]
   Stop value for last energy bin in TeV [100.0]
   Output ASCII file [butterfly.txt] butterfly_src001.txt

and for the second source using

.. code-block:: bash

   $ ctbutterfly
   Input event list, counts cube or observation definition XML file [events.fits] obs_selected.xml
   Source of interest [Crab] Src002
   Input model definition XML file [$CTOOLS/share/models/crab.xml] results.xml
   Start value for first energy bin in TeV [0.1]
   Stop value for last energy bin in TeV [100.0]
   Output ASCII file [butterfly.txt] butterfly_src002.txt

The butterfly diagrams for both sources are displayed in the figure below:

.. figure:: first_spectrum.png
   :width: 600px
   :align: center

   *Butterfly diagrams determined with ctbutterfly and spectral points determined with csspec for Src001 (top) and Src002 (bottom)*

The figure also shows spectral points for each source that were determined
using the :ref:`csspec` script.
You create the spectrum for the first source using

.. code-block:: bash

   $ csspec
   Input event list, counts cube, or observation definition XML file [events.fits] obs_selected.xml
   Input model definition XML file [$CTOOLS/share/models/crab.xml] results.xml
   Source name [Crab] Src001
   Binning algorithm (LIN|LOG|FILE) [LOG]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins (0=unbinned) [20] 10
   Output spectrum file [spectrum.fits] spectrum_src001.fits

and for the second source using

.. code-block:: bash

   $ csspec
   Input event list, counts cube, or observation definition XML file [events.fits] obs_selected.xml
   Input model definition XML file [$CTOOLS/share/models/crab.xml] results.xml
   Source name [Crab] Src002
   Binning algorithm (LIN|LOG|FILE) [LOG]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins (0=unbinned) [20] 10
   Output spectrum file [spectrum.fits] spectrum_src002.fits

The :ref:`csspec` script divided here the data into ten logarithmically
spaced energy bins and determined the source flux in each of the bins using
a maximum likelihood model fit.
