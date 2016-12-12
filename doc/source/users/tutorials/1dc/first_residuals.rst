.. _1dc_first_residuals:

Inspecting the fit residuals
----------------------------

After doing the model fit you should investigate the residuals to verify that
the model components properly describe the observed event distribution.
You do this with the :ref:`csresmap` script by providing the output
:ref:`model definition XML file <glossary_moddef>`
``results.xml`` produced by :ref:`ctlike` as input:

.. code-block:: bash

   $ csresmap
   Input event list, counts cube, or observation definition XML file [events.fits] obs.xml
   Input model definition XML file [$CTOOLS/share/models/crab.xml] results.xml
   First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 0.0
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 0.0
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Coordinate System (CEL|GAL) [CEL] GAL
   Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR]
   Size of the X axis in pixels [200] 500
   Size of the Y axis in pixels [200] 250
   Pixel size (deg/pixel) [0.02] 0.02
   Residual map computation algorithm (SUB|SUBDIV|SUBDIVSQRT) [SUBDIV] SUB
   Output residual map file [resmap.fits]

This produces a residual map that you can display for example with
`ds9 <http://ds9.si.edu>`_.
The figure below shows the result, with the fitted source positions indicated
by white circles.
The residuals near the positions of both point sources are relatively flat,
demonstrating that both excesses are well described by point sources.
Residuals along the Galactic plane indicate that not all sources that are
present in the data are described by the model, which is explained by the
relatively high detection threshold that you used in the source detection
step.

.. figure:: first_skymap_residual.png
   :width: 600px
   :align: center

   *Residual sky map after subtraction of the fitted model*
