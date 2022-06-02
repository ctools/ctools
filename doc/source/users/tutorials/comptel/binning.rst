.. _comptel_binning:

Binning the data
----------------

  .. admonition:: What you will learn

     You will learn how to bin COMPTEL data for an analysis with ctools.


Data need to be binned for the analysis with ctools. The data are binned
in three-dimensional data cubes, called ``DRI`` (Data Required for Imaging).
``DRI`` datasets are spanned by the scatter directions :math:`\chi` and
:math:`\psi` and the Compton scatter angle :math:`\bar{\varphi}`.
Binned data products for an observation include event data (``DRE``),
geometry factors (``DRG``) and a background model (``DRB``). In addition,
a two-dimensional exposure map (``DRX``) is needed. Finally, the binned
data products also include the Instrument Response Function in the so-called
``IAQ`` format, which is a two-dimensional map that provides the interaction
probabilities as function of :math:`\varphi_{\rm geo}` and :math:`\bar{\varphi}`

You bin the data with the :ref:`comobsbin` tool. Binned data products will be
written in a so-called data store so that they can be reused in different analyses.
:ref:`comobsbin` will only create binned data files that do not yet exist
in the datastore. The example below shows how you create binned data products
for viewing period ``0001``.

.. code-block:: bash

   $ comobsbin
   Input observation definition file [obs.xml]
   Response type (MODEL|SIM2|SIM3) [MODEL]
   Algorithm for defining energy bins (FILE|LIN|LOG) [FILE] LOG
   Minimum energy (MeV) (0.05-50.0) [0.75]
   Maximum energy (MeV) (0.05-50.0) [30.0]
   Number of energy bins (1-20) [1] 4
   Phase expression in the format phasemin0-phasemax0;phasemin1-phasemax1;... [NONE]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL]
   Projection method (AIT|AZP|CAR|GLS|MOL|SFL|SIN|STG|TAN) [TAN]
   Output folder for files [$COMDATASTORE] binned
   Output observation definition file [obs_binned.xml]

.. warning::

   Note that the creation of binned data products is relatively time
   consuming, in particular for geometry (``DRG``) and exposure
   files (``DRX``). :ref:`comobsbin` may typically take between
   30 min and 2 hours to process the data for one viewing period,
   depending on the operating system and processor speed.

:ref:`comobsbin` produces on output an
:ref:`observation definition file <glossary_obsdef>`
which in the example is named ``obs_binned.xml``. The content of this file
is shown below. According to the user parameter, the energy interval
0.75-30 MeV was split into four logarithmically spaced energy bins,
and a model Instrument Response Function was computed for each of the
energy bins. :ref:`comobsbin` also generates a Phibar-normalised geometry
function as a first order approximation of the instrumental background model.

.. code-block:: bash

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <observation_list title="observation list">
     <observation name="CRAB" id="vp0001_0_000750-001886keV" instrument="COM">
       <parameter name="DRE" file="binned/vp0001_0_cel_dre_000750-001886keV.fits" />
       <parameter name="DRB" file="binned/vp0001_0_cel_drb_000750-001886keV.fits" />
       <parameter name="DRG" file="binned/vp0001_0_cel_drg.fits" />
       <parameter name="DRX" file="binned/vp0001_0_drx.fits" />
       <parameter name="IAQ" value="binned/iaq_000750-001886keV.fits" />
     </observation>
     <observation name="CRAB" id="vp0001_0_001886-004743keV" instrument="COM">
       <parameter name="DRE" file="binned/vp0001_0_cel_dre_001886-004743keV.fits" />
       <parameter name="DRB" file="binned/vp0001_0_cel_drb_001886-004743keV.fits" />
       <parameter name="DRG" file="binned/vp0001_0_cel_drg.fits" />
       <parameter name="DRX" file="binned/vp0001_0_drx.fits" />
       <parameter name="IAQ" value="binned/iaq_001886-004743keV.fits" />
     </observation>
     <observation name="CRAB" id="vp0001_0_004743-011929keV" instrument="COM">
       <parameter name="DRE" file="binned/vp0001_0_cel_dre_004743-011929keV.fits" />
       <parameter name="DRB" file="binned/vp0001_0_cel_drb_004743-011929keV.fits" />
       <parameter name="DRG" file="binned/vp0001_0_cel_drg.fits" />
       <parameter name="DRX" file="binned/vp0001_0_drx.fits" />
       <parameter name="IAQ" value="binned/iaq_004743-011929keV.fits" />
     </observation>
     <observation name="CRAB" id="vp0001_0_011929-029999keV" instrument="COM">
       <parameter name="DRE" file="binned/vp0001_0_cel_dre_011929-029999keV.fits" />
       <parameter name="DRB" file="binned/vp0001_0_cel_drb_011929-029999keV.fits" />
       <parameter name="DRG" file="binned/vp0001_0_cel_drg.fits" />
       <parameter name="DRX" file="binned/vp0001_0_drx.fits" />
       <parameter name="IAQ" value="binned/iaq_011929-029999keV.fits" />
     </observation>
   </observation_list>
