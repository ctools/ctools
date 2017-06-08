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

   2017-06-08T08:20:32: +=================================+
   2017-06-08T08:20:32: | Maximum likelihood optimisation |
   2017-06-08T08:20:32: +=================================+
   2017-06-08T08:20:46:  >Iteration   0: -logL=242767.783, Lambda=1.0e-03
   2017-06-08T08:21:00:  >Iteration   1: -logL=222518.221, Lambda=1.0e-03, delta=20249.562, step=1.0e+00, max(|grad|)=50533.167265 [Index:13]
   2017-06-08T08:21:14:  >Iteration   2: -logL=218997.783, Lambda=1.0e-04, delta=3520.439, step=1.0e+00, max(|grad|)=-13548.755375 [DEC:1]
   2017-06-08T08:21:29:  >Iteration   3: -logL=217965.881, Lambda=1.0e-05, delta=1031.901, step=1.0e+00, max(|grad|)=-11344.447813 [DEC:1]
   2017-06-08T08:21:43:  >Iteration   4: -logL=217775.661, Lambda=1.0e-06, delta=190.220, step=1.0e+00, max(|grad|)=-7778.127593 [DEC:1]
   2017-06-08T08:21:57:  >Iteration   5: -logL=217725.868, Lambda=1.0e-07, delta=49.793, step=1.0e+00, max(|grad|)=-4796.822100 [DEC:1]
   2017-06-08T08:22:11:  >Iteration   6: -logL=217701.511, Lambda=1.0e-08, delta=24.357, step=1.0e+00, max(|grad|)=-3267.160620 [RA:6]
   2017-06-08T08:22:25:  >Iteration   7: -logL=217688.626, Lambda=1.0e-09, delta=12.885, step=1.0e+00, max(|grad|)=-2438.542445 [RA:6]
   2017-06-08T08:22:40:  >Iteration   8: -logL=217681.375, Lambda=1.0e-10, delta=7.250, step=1.0e+00, max(|grad|)=-1905.466691 [RA:6]
   2017-06-08T08:22:54:  >Iteration   9: -logL=217677.367, Lambda=1.0e-11, delta=4.008, step=1.0e+00, max(|grad|)=-1338.131602 [RA:6]
   2017-06-08T08:23:08:  >Iteration  10: -logL=217675.391, Lambda=1.0e-12, delta=1.977, step=1.0e+00, max(|grad|)=-1117.908664 [RA:6]
   2017-06-08T08:23:22:  >Iteration  11: -logL=217674.202, Lambda=1.0e-13, delta=1.188, step=1.0e+00, max(|grad|)=-683.041271 [RA:6]
   2017-06-08T08:23:36:  >Iteration  12: -logL=217673.712, Lambda=1.0e-14, delta=0.490, step=1.0e+00, max(|grad|)=-481.729534 [RA:6]
   2017-06-08T08:23:50:  >Iteration  13: -logL=217673.495, Lambda=1.0e-15, delta=0.217, step=1.0e+00, max(|grad|)=-385.504388 [RA:6]
   2017-06-08T08:24:04:  >Iteration  14: -logL=217673.327, Lambda=1.0e-16, delta=0.168, step=1.0e+00, max(|grad|)=-314.151335 [RA:6]
   2017-06-08T08:24:18:  >Iteration  15: -logL=217673.231, Lambda=1.0e-17, delta=0.096, step=1.0e+00, max(|grad|)=-227.509872 [RA:6]
   2017-06-08T08:24:32:  >Iteration  16: -logL=217673.170, Lambda=1.0e-18, delta=0.060, step=1.0e+00, max(|grad|)=-176.781601 [RA:6]
   2017-06-08T08:24:47:  >Iteration  17: -logL=217673.138, Lambda=1.0e-19, delta=0.032, step=1.0e+00, max(|grad|)=-124.109305 [RA:6]
   2017-06-08T08:25:01:  >Iteration  18: -logL=217673.124, Lambda=1.0e-20, delta=0.014, step=1.0e+00, max(|grad|)=-88.895131 [RA:6]
   2017-06-08T08:25:15:  >Iteration  19: -logL=217673.117, Lambda=1.0e-21, delta=0.007, step=1.0e+00, max(|grad|)=-63.724112 [RA:6]
   2017-06-08T08:25:29:  >Iteration  20: -logL=217673.112, Lambda=1.0e-22, delta=0.006, step=1.0e+00, max(|grad|)=61.601206 [DEC:7]
   2017-06-08T08:25:43:  >Iteration  21: -logL=217673.107, Lambda=1.0e-23, delta=0.005, step=1.0e+00, max(|grad|)=58.318115 [DEC:7]
   2017-06-08T08:25:57:  >Iteration  22: -logL=217673.101, Lambda=1.0e-24, delta=0.006, step=1.0e+00, max(|grad|)=50.936104 [DEC:7]
   2017-06-08T08:26:12:  >Iteration  23: -logL=217673.097, Lambda=1.0e-25, delta=0.004, step=1.0e+00, max(|grad|)=46.423925 [DEC:7]
   2017-06-08T08:26:26:
   2017-06-08T08:26:26: +=========================================+
   2017-06-08T08:26:26: | Maximum likelihood optimisation results |
   2017-06-08T08:26:26: +=========================================+
   2017-06-08T08:26:26: === GOptimizerLM ===
   2017-06-08T08:26:26:  Optimized function value ..: 217673.097
   2017-06-08T08:26:26:  Absolute precision ........: 0.005
   2017-06-08T08:26:26:  Acceptable value decrease .: 2
   2017-06-08T08:26:26:  Optimization status .......: converged
   2017-06-08T08:26:26:  Number of parameters ......: 16
   2017-06-08T08:26:26:  Number of free parameters .: 10
   2017-06-08T08:26:26:  Number of iterations ......: 23
   2017-06-08T08:26:26:  Lambda ....................: 1e-26
   2017-06-08T08:26:26:  Maximum log likelihood ....: -217673.097
   2017-06-08T08:26:26:  Observed events  (Nobs) ...: 2205594.000
   2017-06-08T08:26:26:  Predicted events (Npred) ..: 2205593.994 (Nobs - Npred = 0.00562562886625528)

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
   Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits
   Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits
   Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits
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
adequately described by a power law model. You should therefore replace the
power law in the
:ref:`model definition file <glossary_moddef>`
by an exponentially cutoff power law, as shown below:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Src001" type="PointSource">
       <spectrum type="ExponentialCutoffPowerLaw">
         <parameter name="Prefactor"    scale="1e-18" value="5.7"  min="1e-07" max="1000.0" free="1"/>
         <parameter name="Index"        scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
         <parameter name="CutoffEnergy" scale="1e7"   value="1.0"  min="0.01"  max="1000.0" free="1"/>
         <parameter name="PivotEnergy"  scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA"  value="266.4045" scale="1" free="1" />
         <parameter name="DEC" value="-28.9945" scale="1" free="1" />
       </spatialModel>
     </source>
     ...
   </source_library>

Fitting this model to the data improves the fit and the resulting butterfly
diagram follows now reasonably well the spectral points:

.. figure:: first_spectrum_cutoff_stacked.png
   :width: 600px
   :align: center

   *Butterfly diagrams determined with ctbutterfly for an exponentially cut-off power law for Src001 (red)*
