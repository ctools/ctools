.. _1dc_first_fitting:

Fitting the model components to the data
----------------------------------------

Now you are ready to fit the model to the data and to determine its maximum
likelihood parameters.
You do this with the :ref:`ctlike` tool that adjusts all parameters in the
:ref:`model definition XML file <glossary_moddef>`
that have the attribute ``free`` set to ``"1"``.
In the current example, the free model parameters are the positions and spectral
parameters of the two point sources and the spectral normalisation of the
background component.
You run the :ref:`ctlike` tool as follows:

.. code-block:: bash

   $ ctlike
   Input event list, counts cube or observation definition XML file [events.fits] obs.xml
   Input model definition XML file [$CTOOLS/share/models/crab.xml] models.xml
   Output model definition XML file [crab_results.xml] results.xml

The tool will take a few minutes to perform the model fitting, and will write
the results into an updated
:ref:`model definition XML file <glossary_moddef>`
containing the fitted model parameters and their statistical uncertainties.
You may inspect the log file ``ctlike.log`` to verify that the model fit
converged properly, as illustrated in the example below:

.. code-block:: bash

   2016-12-02T12:06:03: +=================================+
   2016-12-02T12:06:03: | Maximum likelihood optimisation |
   2016-12-02T12:06:03: +=================================+
   2016-12-02T12:06:19:  >Iteration   0: -logL=13460292.282, Lambda=1.0e-03
   2016-12-02T12:06:34:  >Iteration   1: -logL=13447430.338, Lambda=1.0e-03, delta=12861.944, max(|grad|)=38643.139247 [Index:13]
   2016-12-02T12:06:49:  >Iteration   2: -logL=13444219.743, Lambda=1.0e-04, delta=3210.595, max(|grad|)=18223.630095 [RA:0]
   2016-12-02T12:07:05:  >Iteration   3: -logL=13442740.130, Lambda=1.0e-05, delta=1479.613, max(|grad|)=12560.528406 [RA:0]
   2016-12-02T12:07:20:  >Iteration   4: -logL=13442178.499, Lambda=1.0e-06, delta=561.631, max(|grad|)=5836.127526 [RA:0]
   2016-12-02T12:07:35:  >Iteration   5: -logL=13442096.088, Lambda=1.0e-07, delta=82.411, max(|grad|)=1183.057138 [RA:0]
   2016-12-02T12:07:50:  >Iteration   6: -logL=13442094.149, Lambda=1.0e-08, delta=1.939, max(|grad|)=55.233391 [RA:0]
   2016-12-02T12:08:06:  >Iteration   7: -logL=13442094.147, Lambda=1.0e-09, delta=0.002, max(|grad|)=-3.893441 [RA:6]
   2016-12-02T12:08:21:
   2016-12-02T12:08:21: +=========================================+
   2016-12-02T12:08:21: | Maximum likelihood optimisation results |
   2016-12-02T12:08:21: +=========================================+
   2016-12-02T12:08:21: === GOptimizerLM ===
   2016-12-02T12:08:21:  Optimized function value ..: 13442094.147
   2016-12-02T12:08:21:  Absolute precision ........: 0.005
   2016-12-02T12:08:21:  Acceptable value decrease .: 2
   2016-12-02T12:08:21:  Optimization status .......: converged
   2016-12-02T12:08:21:  Number of parameters ......: 16
   2016-12-02T12:08:21:  Number of free parameters .: 10
   2016-12-02T12:08:21:  Number of iterations ......: 7
   2016-12-02T12:08:21:  Lambda ....................: 1e-10
   2016-12-02T12:08:21:  Maximum log likelihood ....: -13442094.147
   2016-12-02T12:08:21:  Observed events  (Nobs) ...: 1899265.000
   2016-12-02T12:08:21:  Predicted events (Npred) ..: 1899225.997 (Nobs - Npred = 39.0025871018879)

You may also convert the fitted model positions into a `ds9 <http://ds9.si.edu>`_
region file using the :ref:`csmodelinfo` script so that you can overlay the
fit results over a sky map:

.. code-block:: bash

   $ csmodelinfo
   Input model definition XML file [model.xml] results.xml
   Output DS9 region file [ds9.reg] positions.reg

The following image shows a zoom of the sky map that comprises both point
sources, with the initial source positions determined by :ref:`cssrcdetect`
as green crosses and the positions fitted by :ref:`ctlike` as blue circles.
Obviously, the initial positions were already nearby the fitted position,
which is required to assure the proper convergence of the fit.

.. figure:: first_skymap_fitted.png
   :width: 600px
   :align: center

   *Background subtracted sky map of the events recorded during the Galactic Plane Survey around the Galactic Centre with the fitted positions of the sources shown as blue circles*

You can also convert the spectral parameters of the point sources into a
butterfly diagram for each source using the :ref:`ctbutterfly` tool.
The butterfly diagram shows the envelope of all power laws that are
statistically compatible with the data.
You create the butterfly diagram for the first source using:

.. code-block:: bash

   $ ctbutterfly
   Input event list, counts cube or observation definition XML file [events.fits] obs.xml
   Source of interest [Crab] Src001
   Input model definition XML file [$CTOOLS/share/models/crab.xml] results.xml
   Start value for first energy bin in TeV [0.1]
   Stop value for last energy bin in TeV [100.0]
   Output ASCII file [butterfly.txt] butterfly_src001.txt

You create the butterfly diagram for the second source by selecting ``Src002``
as ``"Source of interest"``.
The butterfly diagrams for both sources are displayed in the figure below:

.. figure:: first_spectrum.png
   :width: 600px
   :align: center

   *Butterfly diagrams determined with ctbutterfly and spectral points determined with csspec for Src001 (top) and Src002 (bottom)*

The figure also shows spectral points for each source that were determined
using the :ref:`csspec` script.
You create the spectrum for the first source using:

.. code-block:: bash

   $ csspec
   Input event list, counts cube, or observation definition XML file [events.fits] obs.xml
   Input model definition XML file [results.xml] results.xml
   Source name [Crab] Src001
   Binning algorithm (LIN|LOG|FILE) [LOG]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins (0=unbinned) [10]
   Output spectrum file [spectrum.fits] spectrum_src001.fits

Alike the :ref:`ctbutterfly` tool, you have to run the :ref:`csspec` script
for each source for which you want to obtain spectral points.
The script will divide the data into a number of logarithmically spaced energy
bins and determine the flux in each of the bins using a maximum likelihood
model fit.
Ten energy bins were used in this example.
