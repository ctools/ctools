.. _sec_stacked_reponse:

Pre-computing the response
--------------------------

  .. admonition:: What you will learn

     You will learn how to **compute the instrument response for a counts
     cube**. Pre-computing the instrument response will speed-up
     the model fitting later.

To speed-up the model fitting it is recommended to pre-compute the instrument
response for the counts cube. Specifically, you should compute an exposure
cube, a point spread function cube and a background cube.

The
:ref:`exposure cube <glossary_expcube>`
is computed using the :ref:`ctexpcube` tool and provides the effective area
multiplied by the livetime of the observation.
You run :ref:`ctexpcube` as follows:

.. code-block:: bash

   $ ctexpcube
   Input event list or observation definition XML file [NONE] selected_events.fits
   Calibration database [prod2]
   Instrument response function [South_0.5h]
   Input counts cube file to extract exposure cube definition [NONE] cntcube.fits
   Output exposure cube file [expcube.fits]

This produces the FITS file ``expcube.fits`` that contains the exposure as
function of sky position and energy. In the example above you extracted the
binning of the exposure cube from the counts cube, which considerably reduces
the number of parameters that are queried by the tool.

  .. note::

     The binning of the exposure cube does not need to correspond to the
     binning of the counts cube. In any case, exposure values will be
     determined by interpolation from the values stored in the exposure cube
     file. The same is true for the point spread function and background cubes
     that are described below.

Next, we use the :ref:`ctpsfcube` tool to compute the
:ref:`point spread function (PSF) cube <glossary_psfcube>` for the counts
cube.
You run :ref:`ctpsfcube` as follows:

.. code-block:: bash

   $ ctpsfcube
   Input event list or observation definition XML file [NONE] selected_events.fits
   Calibration database [prod2]
   Instrument response function [South_0.5h]
   Input counts cube file to extract PSF cube definition [NONE]
   First coordinate of image center in degrees (RA or galactic l) (0-360) [83.63]
   Second coordinate of image center in degrees (DEC or galactic b) (-90-90) [22.01]
   Projection method (AIT|AZP|CAR|GLS|MER|MOL|SFL|SIN|STG|TAN) [CAR]
   Coordinate system (CEL - celestial, GAL - galactic) (CEL|GAL) [CEL]
   Image scale (in degrees/pixel) [1.0]
   Size of the X axis in pixels [10]
   Size of the Y axis in pixels [10]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Number of energy bins [20]
   Output PSF cube file [psfcube.fits]

This produces the FITS file ``psfcube.fits`` that contains the point spread
function as function of sky position and energy. You may have noted in the
example that the definiton of the
:ref:`point spread function cube <glossary_psfcube>`
has not been extracted from the counts cube, since this would lead to a
large FITS file on output.
The point spread function varies in fact only slowly over the field of view
of the camera, and consequently it is sufficient to sample that variation
at a large spatial scale of typically one degree.

Finally, we use the :ref:`ctbkgcube` tool to compute the
:ref:`background cube <glossary_bkgcube>`.
You run :ref:`ctbkgcube` as follows:

.. code-block:: bash

   $ ctbkgcube
   Input event list or observation definition XML file [NONE] selected_events.fits
   Calibration database [prod2]
   Instrument response function [South_0.5h]
   Input counts cube file to extract background cube definition [NONE] cntcube.fits
   Input model definition XML file [NONE] $CTOOLS/share/models/crab.xml
   Output background cube file [bkgcube.fits]
   Output model definition XML file [NONE] models.xml

This produces the FITS file ``bkgcube.fits`` that contains the predicted
background rate as function of sky position and energy.
The tool also produces the
:ref:`model definition file <glossary_moddef>`
``models.xml``
on output that will serve as input for the maximum likelihood analysis that
will follow.
The file is a copy of the input
:ref:`model definition file <glossary_moddef>`
``$CTOOLS/share/models/crab.xml``
where the input background model has been replaced by a background model of
type ``CTACubeBackground``. Below is the content of the ``models.xml`` file:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Crab" type="PointSource">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor" value="5.7" error="0" scale="1e-16" min="1e-07" max="1000" free="1" />
         <parameter name="Index" value="2.48" error="0" scale="-1" min="0" max="5" free="1" />
         <parameter name="PivotEnergy" value="0.3" scale="1000000" min="0.01" max="1000" free="0" />
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA" value="83.6331" scale="1" min="-360" max="360" free="0" />
         <parameter name="DEC" value="22.0145" scale="1" min="-90" max="90" free="0" />
       </spatialModel>
     </source>
     <source name="BackgroundModel" type="CTACubeBackground" instrument="CTA,HESS,MAGIC,VERITAS">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor" value="1" error="0" scale="1" min="0.01" max="100" free="1" />
         <parameter name="Index" value="0" error="0" scale="1" min="-5" max="5" free="1" />
         <parameter name="PivotEnergy" value="1" scale="1000000" free="0" />
       </spectrum>
     </source>
   </source_library>
