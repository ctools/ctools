.. _sec_binning_cta:

Binning CTA data
~~~~~~~~~~~~~~~~

As next analysis step you will bin the data in a counts cube using 
:ref:`ctbin`.
A counts cube is a 3 dimensional data cube, spanned by
Right Ascension (or Galactic longitude), Declination (or Galactic latitude),
and energy (typically logarithmically spaced, but this is under user
control).

:ref:`ctbin` is executed by typing:

.. code-block:: bash

  $ ctbin
  Input event list or observation definition XML file [events.fits] 
  First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 
  Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 
  Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR] 
  Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] 
  Image scale (in degrees/pixel) [0.02] 
  Size of the X axis in pixels [200] 
  Size of the Y axis in pixels [200] 
  Algorithm for defining energy bins (FILE|LIN|LOG) [LOG] 
  Start value for first energy bin in TeV [0.1] 
  Stop value for last energy bin in TeV [100.0] 
  Number of energy bins [20] 
  Output counts cube file [cntcube.fits] 

The counts cube will be centred on the location of the Crab (Right Ascension
83.63 degrees, Declination 22.01 degrees) and will be aligned in celestial
coordinates. A cartesian projection has been selected. The counts cube has 
200 x 200 spatial pixels of 0.02 x 0.02 degrees in size, hence it covers a 
total area of 4 x 4 degrees.

The counts cube will contain 20 maps, which are logarithmically spaced
in energy, and which cover the energy range from 0.1 TeV to 100 TeV. In this
example, the counts cube will be saved as ``cntcube.fits`` in the working
directory. In addition to the counts cube, that is stored as the primary
image extension, the FITS file also contains an extension named ``EBOUNDS``
that defines the energy boundaries that were used, and an extension ``GTI``
that defines the Good Time Intervals that have been used. The following
image shows the resulting FITS file. The ``EBOUNDS`` table has 20 rows, one
for each energy bin, while the ``GTI`` table has just a single row, indicating
the start and stop time of the simulated data.

.. figure:: cntmap-fits.jpg
   :width: 600px
   :align: center

   *Counts cube FITS file*


An image of the first bin, covering the energy range 100 - 141 GeV, is 
shown below:

.. figure:: cntmap-map.jpg
   :height: 400px
   :align: center

   *Counts cube for first energy bin*

For illustration, the last few lines of the log file ``ctbin.log`` are 
reproduced below:

.. code-block:: xml

  2015-12-07T20:52:41: +=================+
  2015-12-07T20:52:41: | Bin observation |
  2015-12-07T20:52:41: +=================+
  2015-12-07T20:52:41: === CTA observation ===
  2015-12-07T20:52:41:  Events in list ............: 22685
  2015-12-07T20:52:41:  Events in cube ............: 18195
  2015-12-07T20:52:41:  Event bins outside RoI ....: 0
  2015-12-07T20:52:41:  Events outside cube area ..: 4490
  2015-12-07T20:52:41:  Events outside energy bins : 0

From the 22685 events that have been simulated and stored in the 
``events.fits`` file, 18195 lie within the cube boundaries and are thus put
into the resulting counts cube. The counts cube is stored in a cartesian
projection in a World Coordinate System (WCS) compliant format.
