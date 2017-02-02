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

The tool will take several minutes to perform the model fitting, and will write
the results into an updated
:ref:`model definition XML file <glossary_moddef>`
containing the fitted model parameters and their statistical uncertainties.
You may inspect the log file ``ctlike.log`` to verify that the model fit
converged properly, as illustrated in the example below:

.. code-block:: bash

   2017-02-02T15:11:01: +=================================+
   2017-02-02T15:11:01: | Maximum likelihood optimisation |
   2017-02-02T15:11:01: +=================================+
   2017-02-02T15:12:10:  >Iteration   0: -logL=39119340.262, Lambda=1.0e-03
   2017-02-02T15:13:12:  >Iteration   1: -logL=39113674.616, Lambda=1.0e-03, delta=5665.646, max(|grad|)=25710.424562 [Index:13]
   2017-02-02T15:14:14:  >Iteration   2: -logL=39111144.224, Lambda=1.0e-04, delta=2530.392, max(|grad|)=-17928.193858 [DEC:1]
   2017-02-02T15:15:16:  >Iteration   3: -logL=39109792.705, Lambda=1.0e-05, delta=1351.519, max(|grad|)=-10987.858986 [DEC:1]
   2017-02-02T15:16:19:  >Iteration   4: -logL=39109377.992, Lambda=1.0e-06, delta=414.713, max(|grad|)=-4554.316982 [DEC:1]
   2017-02-02T15:17:20:  >Iteration   5: -logL=39109337.418, Lambda=1.0e-07, delta=40.574, max(|grad|)=-910.537530 [DEC:1]
   2017-02-02T15:18:21:  >Iteration   6: -logL=39109336.839, Lambda=1.0e-08, delta=0.579, max(|grad|)=-70.481486 [DEC:1]
   2017-02-02T15:19:26:  >Iteration   7: -logL=39109336.837, Lambda=1.0e-09, delta=0.002, max(|grad|)=-3.284527 [DEC:1]
   2017-02-02T15:20:30:
   2017-02-02T15:20:30: +=========================================+
   2017-02-02T15:20:30: | Maximum likelihood optimisation results |
   2017-02-02T15:20:30: +=========================================+
   2017-02-02T15:20:30: === GOptimizerLM ===
   2017-02-02T15:20:30:  Optimized function value ..: 39109336.837
   2017-02-02T15:20:30:  Absolute precision ........: 0.005
   2017-02-02T15:20:30:  Acceptable value decrease .: 2
   2017-02-02T15:20:30:  Optimization status .......: converged
   2017-02-02T15:20:30:  Number of parameters ......: 16
   2017-02-02T15:20:30:  Number of free parameters .: 10
   2017-02-02T15:20:30:  Number of iterations ......: 7
   2017-02-02T15:20:30:  Lambda ....................: 1e-10
   2017-02-02T15:20:30:  Maximum log likelihood ....: -39109336.837
   2017-02-02T15:20:30:  Observed events  (Nobs) ...: 7736578.000
   2017-02-02T15:20:30:  Predicted events (Npred) ..: 7736555.998 (Nobs - Npred = 22.002333926037)

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
