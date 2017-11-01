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

In the spirit of the classical analysis we have chosen to derive the background
using the RING method. You need to choose inner/outer radii such that you avoid
emission from a source when deriving the background. This means that the inner
radius must be larger than the source's size, or the instrument PSF for a
pointlike source. This has produced a FITS file ``skymap.fits`` that contains
three images of the region around the source. The primary image shows the excess
counts, i.e., the total number of counts minus the estimated background counts.
The BACKGROUND image shows the number of estimated background counts. Finally,
the SIGNIFICANCE image shows the significance of the excess, calculated
according to Li&Ma
