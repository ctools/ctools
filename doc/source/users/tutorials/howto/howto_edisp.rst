.. _1dc_howto_edisp:

How to take the energy dispersion into account?
-----------------------------------------------

  .. admonition:: What you will learn

     You will learn how to **take the energy dispersion into account** when
     fitting a model to the data using a maximum likelihood analysis.

     Although the effect of the energy dispersion is often neglegible there
     may be cases where you want to consider energy dispersion in an analysis,
     for example if you are analysing the data down to very low energies.

  .. warning::
     Energy dispersion is fully implemented in ctools but is computationally
     intensive. So be aware that the **tools and scripts will take a substantial
     amount of computing time**.

If you are doing a stacked analysis the first thing you need is an energy
dispersion cube. You generate this cube using the :ref:`ctedispcube` tool:

.. code-block:: bash

   $ ctedispcube
   Input event list or observation definition XML file [NONE] obs_selected.xml
   Input counts cube file to extract energy dispersion cube definition [NONE]
   First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 0.0
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 0.0
   Projection method (AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN) [CAR]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] GAL
   Image scale (in degrees/pixel) [1.0]
   Size of the X axis in pixels [10]
   Size of the Y axis in pixels [10]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins [20] 30
   Output energy dispersion cube file [edispcube.fits]

This produces an
:ref:`energy dispersion cube <glossary_edispcube>`
FITS file that contains the weighted energy dispersion as function of
sky position and energy.
Since the energy dispersion varies only slowly over the field of view it is
sufficient to sample that variation at a large spatial scale of typically one
degree.

Now you can run the :ref:`ctlike` tool by adding the attribute ``edisp=yes``
to the command line. As you see below, this instructs :ref:`ctlike` to query
for the energy dispersion cube and to take the energy dispersion into account
during the fit:

.. code-block:: bash

   $ ctlike edisp=yes
   Input event list, counts cube or observation definition XML file [events.fits] cntcube.fits
   Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits
   Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits
   Input energy dispersion cube file (only needed for stacked analysis) [NONE] edispcube.fits
   Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits
   Input model definition XML file [$CTOOLS/share/models/crab.xml] stacked_models.xml
   Output model definition XML file [crab_results.xml] stacked_results_edisp.xml

If you are doing an unbinned analysis you do not need to generate an energy
dispersion cube and you can directly run :ref:`ctlike` with the ``edisp=yes``
attribute:

.. code-block:: bash

   $ ctlike edisp=yes
   Input event list, counts cube or observation definition XML file [events.fits] obs_selected.xml
   Input model definition XML file [$CTOOLS/share/models/crab.xml] models.xml
   Output model definition XML file [crab_results.xml] unbinned_models_edisp.xml

The following tools and scripts do accept the optional ``edisp=yes`` attribute:

* :ref:`ctbutterfly`
* :ref:`cterror`
* :ref:`ctlike`
* :ref:`ctmodel`
* :ref:`ctprob`
* :ref:`cttsmap`
* :ref:`ctulimit`
* :ref:`cslightcrv`
* :ref:`csphasecrv`
* :ref:`cspull`
* :ref:`csresmap`
* :ref:`cssens`
* :ref:`csspec`
* :ref:`cstsdist`
* :ref:`cstsmapsplit`

