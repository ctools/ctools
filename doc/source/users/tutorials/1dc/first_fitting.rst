.. _1dc_first_fitting:

Fitting the model components to the data
----------------------------------------

.. code-block:: bash

   $ ctlike
   Input event list, counts cube or observation definition XML file [events.fits] obs.xml
   Input model XML file [$CTOOLS/share/models/crab.xml] models.xml
   Output model XML file [crab_results.xml] results.xml


.. figure:: first_skymap_fitted.png
   :width: 600px
   :align: center

   *Background subtracted sky map of the events recorded during the Galactic Plane Survey around the Galactic Centre with the fitted positions of the sources shown as blue circles*

.. code-block:: bash

   $ ctbutterfly
   Input event list, counts cube or observation definition XML file [events.fits] obs.xml
   Source of interest [Crab] Src001
   Input model XML file [$CTOOLS/share/models/crab.xml] ctlike_unbinned_bkg_3ptsrc.xml
   Start value for first energy bin in TeV [0.1]
   Stop value for last energy bin in TeV [100.0]
   Output ASCII file [butterfly.txt] butterfly_src001.txt

.. code-block:: bash

   $ csspec
   Input event list, counts cube, or observation definition XML file [events.fits] obs.xml
   Input model XML file [results.xml] ctlike_unbinned_bkg_3ptsrc.xml
   Source name [Crab] Src001
   Binning algorithm (LIN|LOG|FILE) [LOG]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins (0=unbinned) [10]
   Output spectrum file [spectrum.fits] spectrum_src001.fits

.. figure:: first_spectrum.png
   :width: 600px
   :align: center

   *Butterfly diagrams determined with ctbutterfly and spectral points determined with csspec for Src001 (top) and Src002 (bottom)*

