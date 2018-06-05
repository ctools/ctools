.. _1dc_sky_map:

Generating a sky map from the events
------------------------------------

  .. admonition:: What you will learn

     You will learn how to **generate a sky map** from the selected events.

     A background model based on the background template included in the
     :ref:`instrument response functions <glossary_irf>`
     may be subtracted from the sky map to display only celestial excess events,
     but be aware that for the moment the background subtraction is
     computationally intensive.

After having selected the observations and the events you can begin your
analysis.

The first thing you want to do is to visualise the content of the
observations by generating a sky map from the events. You do this with the
:ref:`ctskymap` tool by typing:

.. code-block:: bash

   $ ctskymap
   Input event list or observation definition XML file [events.fits] obs_selected.xml
   First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 0.0
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 0.0
   Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] GAL
   Image scale (in degrees/pixel) [0.02]
   Size of the X axis in pixels [200] 400
   Size of the Y axis in pixels [200] 400
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Background subtraction method (NONE|IRF) [NONE]
   Output skymap file [skymap.fits]

This generates a sky map centred on the Galactic Centre in Galactic
coordinates using a cartesian projection.
The sky map is 8 degrees wide and 8 degrees high, with an image scale of
0.02 degrees per pixel.
All events between 100 GeV and 100 TeV are collected in the sky map.
The sky map is written into the FITS file ``skymap.fits`` that is created in
the working directory.
The sky map, displayed using
`ds9 <http://ds9.si.edu>`_ with a square-root color scaling and small Gaussian
kernel smoothing applied, is shown below:

.. figure:: first_skymap.png
   :width: 400px
   :align: center

   *Sky map of the events recorded around the Galactic Centre during the Galactic Plane Survey*

The sky map shows a wide-spread distribution of events with several sources
superimposed.
Many of the events originate from an irreducable background that hampers the
recognition of the gamma-ray sources.
To describe the irreducable background in the CTA data, templates of the
background event distribution are included in the
:ref:`instrument response functions <glossary_irf>`.
These templates can be used by :ref:`ctskymap` to subtract the irreducable
background contribution from the sky map.
The background subtraction is enabled by running the :ref:`ctskymap` with
the background subtraction method set to ``IRF``, as shown in the following
example:

.. code-block:: bash

   $ ctskymap
   Input event list or observation definition XML file [obs_selected.xml]
   First coordinate of image center in degrees (RA or galactic l) (0-360) [0.0]
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [0.0]
   Projection method (AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN) [CAR]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [GAL]
   Image scale (in degrees/pixel) [0.02]
   Size of the X axis in pixels [400]
   Size of the Y axis in pixels [400]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Background subtraction method (NONE|IRF|RING) [NONE] IRF
   Output skymap file [skymap.fits] skymap_irf.fits

The figure below shows the resulting sky map.
The irreducable background has been subtracted from the sky map and the sources
of gamma-ray emission are now clearly discernable.

.. figure:: first_skymap_irf.png
   :width: 400px
   :align: center

   *IRF background subtracted sky map of the events recorded around the Galactic Centre during the Galactic Plane Survey*

In reality, however, the distribution of the irreducable background may only be
purely known, and it may be necessary to estimate the background from the data
themselves.
This can be achieved by the so called ring-background method, that estimates the
background rate from a ring around a given position.
The ring-background can be used by specifying ``RING`` as the background
subtraction method. :ref:`ctskymap` will then query for the source region and
ring radii.

.. warning::
   The estimation of the ring-background by :ref:`ctskymap` is time consuming.
   Be prepared to wait for about half an hour for the results on Max OS X.

To create a ring-background subtracted sky map, type the following:

.. code-block:: bash

   $ ctskymap
   Input event list or observation definition XML file [obs_selected.xml]
   First coordinate of image center in degrees (RA or galactic l) (0-360) [0.0]
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [0.0]
   Projection method (AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN) [CAR]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [GAL]
   Image scale (in degrees/pixel) [0.02]
   Size of the X axis in pixels [400]
   Size of the Y axis in pixels [400]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Background subtraction method (NONE|IRF|RING) [IRF] RING
   Source region radius for estimating on-counts (degrees) [0.2]
   Inner ring radius (degrees) [0.6]
   Outer ring radius (degrees) [0.8]
   Output skymap file [skymap_irf.fits] skymap_ring.fits

The figure below shows the resulting sky map.

.. figure:: first_skymap_ring.png
   :width: 400px
   :align: center

   *Ring background subtracted sky map of the events recorded around the Galactic Centre during the Galactic Plane Survey*

Note the dark zones above and below the Galactic centre that due to pixels
that have bright sources in the rings that are used for their background
estimation.
This dark zones can be avoided by specifying so-called exclusion regions for
the bright sources in the sky map (see the section :ref:`Generating a sky map <start_skymap>`).
