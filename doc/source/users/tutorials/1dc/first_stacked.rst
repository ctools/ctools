.. _1dc_first_stacked:

Stacking the data
-----------------

  .. admonition:: What you will learn

     You will learn how to **prepare the data for analysis** by stacking all
     events into a counts cube and by computing the effective
     :ref:`instrument response function <glossary_irf>`
     for this counts cube.

     Please note that for the moment we do not consider the energy dispersion
     since for many cases the effect of the energy dispersion is negligible,
     and taking energy dispersion into account is computationally intensive. If
     you want to learn how to take the energy dispersion into account please
     read :ref:`this section<1dc_howto_edisp>`.

To analyse the selected observations, we recommend to stack the events into
a counts cube. You do this using the :ref:`ctbin` tool:

.. code-block:: bash

   $ ctbin
   Input event list or observation definition XML file [events.fits] obs_selected.xml
   First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 0.0
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 0.0
   Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] GAL
   Image scale (in degrees/pixel) [0.02]
   Size of the X axis in pixels [200] 300
   Size of the Y axis in pixels [200] 300
   Algorithm for defining energy bins (FILE|LIN|LOG) [LOG]
   Start value for first energy bin in TeV [0.1]
   Stop value for last energy bin in TeV [100.0]
   Number of energy bins (1-200) [20] 30
   Output counts cube file [cntcube.fits]

The :ref:`ctbin` tool creates a 3-dimensional counts cube in Galactic
coordinates, centred on the Galactic centre and 6 degrees x 6 degrees wide,
with 30 logarithmically spaced energy bins between 100 GeV and 100 TeV.

Since the :ref:`ctbin` tool combines observations that potentially have
different
:ref:`instrument response functions <glossary_irf>`
and exposure times into a single counts cube, you have to compute the
effective
:ref:`instrument response function <glossary_irf>`
for this counts cube before you can analyse the data.
For each component of the
:ref:`instrument response function <glossary_irf>`
there is a specific tool to perform this computation.

First, :ref:`ctexpcube` computes the exposure of the stacked counts cube
which is the effective area multiplied by the livetime for each observation.
You run :ref:`ctexpcube` as follows:

.. code-block:: bash

   $ ctexpcube
   Input event list or observation definition XML file [NONE] obs_selected.xml
   Input counts cube file to extract exposure cube definition [NONE] cntcube.fits
   Output exposure cube file [expcube.fits]

This produces an
:ref:`exposure cube <glossary_expcube>`
FITS file that contains the exposure as function of sky position and energy.
In the example above you extracted the binning of the exposure cube from
the counts cube, which considerably reduces the number of User parameters
that are queried by the tool.

.. note::

   The binning of the exposure cube does not need to correspond to the binning
   of the counts cube. In any case, exposure values will be determined by
   interpolation from the values stored in the exposure cube file. The same
   is true for the point spread function and background cubes that are
   described below, or the energy dispersion cube that is described
   :ref:`here<1dc_howto_edisp>`.

Next, :ref:`ctpsfcube` computes the weighted Point Spread Function of the
stacked counts cube.
You run :ref:`ctpsfcube` as follows:

.. code-block:: bash

   $ ctpsfcube
   Input event list or observation definition XML file [NONE] obs_selected.xml
   Input counts cube file to extract PSF cube definition [NONE]
   First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63] 0.0
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01] 0.0
   Projection method (AIT|AZP|CAR|MER|MOL|STG|TAN) [CAR]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL] GAL
   Image scale (in degrees/pixel) [1.0]
   Size of the X axis in pixels [10]
   Size of the Y axis in pixels [10]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins [20] 30
   Output PSF cube file [psfcube.fits]

This produces a
:ref:`point spread function cube <glossary_psfcube>`
FITS file that contains the weighted point spread function as function of
sky position and energy.
You may have noted in the example that the definiton of the
:ref:`point spread function cube <glossary_psfcube>`
has not been extracted from the counts cube, since this would lead to a
large FITS file on output.
The point spread function varies in fact only slowly over the field of view
of the camera, and consequently it is sufficient to sample that variation
at a large spatial scale of typically one degree.

Finally, :ref:`ctbkgcube` computes a
:ref:`background cube <glossary_bkgcube>`
that predicts the number of background events in the counts cube.
You run :ref:`ctbkgcube` as follows:

.. code-block:: bash

   $ ctbkgcube
   Input event list or observation definition XML file [NONE] obs_selected.xml
   Input counts cube file to extract background cube definition [NONE] cntcube.fits
   Input model definition XML file [NONE] models.xml
   Output background cube file [bkgcube.fits]
   Output model definition XML file [NONE] stacked_models.xml

This produces a
:ref:`background cube <glossary_bkgcube>`
FITS file that contains the predicted background rate as function of sky
position and energy.
The tool also produces a
:ref:`model definition file <glossary_moddef>`
on output that can serve as input for a maximum likelihood analysis. The file
is a copy of the input
:ref:`model definition file <glossary_moddef>`
where the input background model will be replaced by a background model of
type ``CTACubeBackground``. Below is the content of the
:ref:`model definition file <glossary_moddef>`
that was generated by :ref:`ctbkgcube`.

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Src001" type="PointSource">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor" value="1" error="0" scale="5.7e-18" min="0" free="1" />
         <parameter name="Index" value="1" error="-0" scale="-2.48" min="-4.03225806451613" max="4.03225806451613" free="1" />
         <parameter name="PivotEnergy" value="1" scale="300000" free="0" />
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA" value="266.424004498437" error="0" scale="1" free="1" />
         <parameter name="DEC" value="-29.0049010253548" error="0" scale="1" free="1" />
       </spatialModel>
     </source>
     <source name="Src002" type="PointSource">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor" value="1" error="0" scale="5.7e-18" min="0" free="1" />
         <parameter name="Index" value="1" error="-0" scale="-2.48" min="-4.03225806451613" max="4.03225806451613" free="1" />
         <parameter name="PivotEnergy" value="1" scale="300000" free="0" />
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA" value="266.831945177213" error="0" scale="1" free="1" />
         <parameter name="DEC" value="-28.1460284439951" error="0" scale="1" free="1" />
       </spatialModel>
     </source>
     <source name="BackgroundModel" type="CTACubeBackground">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor" value="1" error="0" scale="1" min="0.01" max="100" free="1" />
         <parameter name="Index" value="0" error="0" scale="1" min="-5" max="5" free="1" />
         <parameter name="PivotEnergy" value="1" scale="1000000" free="0" />
       </spectrum>
     </source>
   </source_library>
