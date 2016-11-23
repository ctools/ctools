.. _1dc_sky_map:

Generating a sky map from the events
------------------------------------

Once the observations are selected, the first thing that you want to do is to
visualize the content of the observations by generating a sky map from the
events.
You do this with the :ref:`ctskymap` tool by typing:

.. code-block:: bash

   $ ctskymap
   Input event list or observation definition XML file [events.fits] obs.xml
   First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 0.0
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 0.0
   Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] GAL
   Image scale (in degrees/pixel) [0.02]
   Size of the X axis in pixels [200] 500
   Size of the Y axis in pixels [200] 250
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Output skymap file [skymap.fits]

We generate here a sky map centred on the Galactic Centre in Galactic
coordinates using a cartesian projection.
The sky map is 10° wide and 5° high, with an image scale of 0.02° per pixel.
All events between 0.1 and 100 TeV have been collected in the skymap.
The :ref:`ctskymap` tool writes the result into the FITS file ``skymap.fits``
that will be created in the ``my_first_analysis`` folder.
The FITS file can be displayed using for example ds9 (see below):

.. figure:: first_skymap.png
   :width: 600px
   :align: center

   *Sky map of the events recorded during the Galactic Plane Survey around the Galactic Centre*

Bla
