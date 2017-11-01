.. _classical_analysis:

Performing a classical analysis
--------------------------------

  .. admonition:: What you will learn

     You will learn how to perform a *classical* atmospheric Cherenkov analysis
     of a source. In thys type of analysis you do not need to have a 3D (spatial
     and spectral) model of the residual hadronic background. First you create a
     skymap to identify a source region. Then, you extract the measured counts
     from the source region and from one or more background regions, and perform
     a 1D maximum likelihood analysis to derive the source's spectrum.

We will assume that you have selected the observations available around your
source of interest. You can do this using the :ref:`csobsselect` script.

We will derive the spectrum of a source with known position, i.e., the supernova
remnants Cassiopeia A. First, you want to create a skymap to identify the source
region. You can do this using the :ref:`ctskymap` tool.

.. code-block:: bash

    $ ctskymap
    Input event list or observation definition XML file [events.fits] obs.xml
    First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 350.85
    Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 58.815
    Projection method (AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN) [CAR] TAN
    Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL]
    Image scale (in degrees/pixel) [0.02]
    Size of the X axis in pixels [200] 250
    Size of the Y axis in pixels [200] 250
    Lower energy limit (TeV) [0.1]
    Upper energy limit (TeV) [100.0] 50.
    Background subtraction method (NONE|IRF|RING) [NONE] RING
    Source region radius for estimating on-counts (degrees) [0.0] 0.2
    Inner ring radius (degrees) [0.6]
    Outer ring radius (degrees) [0.8]
    Output skymap file [skymap.fits]

In the spirit of the classical analysis we have chosen to derive the background
using the RING method (e.g., Berge et al. 2007 A&A 466 1219B). The RING method
is fast and robust against linear gradients of the background rates in the
the camera, but requires a model (from IRFs) of the background acceptance as a
function of position in the sky. You need to choose inner/outer radii such that
you avoid emission from a source when deriving the background. This means that
the inner radius must be much larger than the source's size, or the instrument
PSF for a pointlike source. This has produced a FITS file ``skymap.fits`` that
contains three images of the region around the source. The primary image shows
the excess counts, i.e., the total number of counts minus the estimated
background counts. The BACKGROUND image shows the number of estimated background
counts. Finally, the SIGNIFICANCE image shows the significance of the excess,
calculated according to Li&Ma 1983 ApJ 272 317, Eq. 17.

.. figure:: classic_analysis_skymap.png
   :width: 400px
   :align: center

   *Sky map of the significance of a gamma-ray excess around Cas A. The green circle shows a circular region with 0.2 deg radius centered at the source's position.*

You can note several features in the skymap. There is significant emission
around Cas A extending beyond the 0.2 deg region that we have computed
the source counts from. Also, there is a ring with negative significance (i.e.,
a count deficit) at offsets between 0.6 deg and 0.8 deg from the source. This is
an artefact due to the fact that when computing the background in this region
the region around Cas A was falling into the RING used for the background
estimation.

We will recompute the skymap addressing these issues. We will enlarge the source
region to a radius of 0.25 deg. We will also provide an input exclusion region
to avoid Cas A when calculating the background. To this end you can use two
exclusion formats: a ds9 region file, or a FITS WCS map. This is what a ds9
region file to excluce Cas A looks like.

.. code-block:: bash

    # Region file format: DS9 version 4.1
    global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
    fk5
    circle(350.85,58.815,900.000")

In fact we could have done this from the beginning since Cas A is a know source.
In general you will need to iterate until you have found all the significant
gamma-ray emission regions and added them to the exlusion regions/map, which is
then necessary for spectral extraction.

We rerun :ref:`ctskymap` with the new parameters.

.. code-block:: bash

    $ ctskymap inexclusion=CasA-exclusion.reg
    Input event list or observation definition XML file [obs.xml]
    First coordinate of image center in degrees (RA or galactic l) (0-360) [350.85]
    Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [58.815]
    Projection method (AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN) [TAN]
    Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL]
    Image scale (in degrees/pixel) [0.02]
    Size of the X axis in pixels [250]
    Size of the Y axis in pixels [250]
    Lower energy limit (TeV) [0.1]
    Upper energy limit (TeV) [50.]
    Background subtraction method (NONE|IRF|RING) [RING]
    Source region radius for estimating on-counts (degrees) [0.2] 0.25
    Inner ring radius (degrees) [0.6]
    Outer ring radius (degrees) [0.8] 0.85
    Output skymap file [skymap.fits] skymap-exclusion.fits

Below you can see the new significance map with the source/exclusion region.

.. figure:: classic_analysis_skymap_exclusion.png
   :width: 400px
   :align: center

   *Sky map of the significance of a gamma-ray excess around Cas A. The green circle shows a circular region with 0.25 deg radius centered at the source's position, that is excluded from the background estimation.*

For a classical spectral analysis we need to derive count spectra for the source
region and for background regions. This is accomplished by the :ref:`csphagen`
script. This script saves the source (On) and background (Off) count spectra
in `OGIP format <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node5.html>`_,
along with the instrument response refashioned according to this format
conventions.