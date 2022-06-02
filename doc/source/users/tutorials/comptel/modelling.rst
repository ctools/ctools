.. _comptel_modelling:

Modeling the data
-----------------

  .. admonition:: What you will learn

     You will learn how to generate a
     :ref:`model definition file <glossary_moddef>`
     for the analysis of COMPTEL data.


As next step you need to create a model to describe the COMPTEL data.
This step is accomplished using the :ref:`comobsmodel` script. The
example below shows how you can generate a model for a point source at
the position of the Crab.

.. code-block:: bash

   $ comobsmodel
   Input observation definition file [obs.xml] obs_binned.xml
   Right Ascension of point source (deg) [NONE] 83.6331
   Declination of point source (deg) [NONE] 22.0145
   Name of point source [NONE] Crab
   Bremsstrahlung component (NONE|MAP|CUBE) [NONE]
   Inverse Compton component (NONE|MAP|CUBE) [NONE]
   Isotropic component (NONE|CONST|CONSTFIX) [NONE]
   Output model definition file [models.xml]

:ref:`comobsmodel` produces on output a
:ref:`model definition file <glossary_moddef>`
which in the example is named ``models.xml``. The content of this file
is shown below. The model comprises a point source with power law spectrum at the
position of the Crab and a background model component for each of the four energy
bins. By default the background model is of type ``DRBPhibarBins`` which defines
a free scaling factor for each of the :math:`\bar{\varphi}` layers of the COMPTEL data space.

.. code-block:: bash

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Crab" type="PointSource" tscalc="1">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor" value="1" error="0" scale="0.002" min="5e-23" free="1" />
         <parameter name="Index" value="1" error="-0" scale="-2" min="-5" max="5" free="1" />
         <parameter name="PivotEnergy" value="1" scale="1" free="0" />
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA" value="83.6331" error="0" scale="1" free="1" />
         <parameter name="DEC" value="22.0145" error="0" scale="1" free="1" />
       </spatialModel>
     </source>
     <source name="Background_vp0001_0_000750-001886keV" type="DRBPhibarBins" instrument="COM" id="vp0001_0_000750-001886keV">
       <parameter name="Normalization" value="1" scale="1" min="0" max="1000" free="0" />
       <parameter name="Normalization" value="1" scale="1" min="0" max="1000" free="0" />
       <parameter name="Normalization" value="1" scale="1" min="0" max="1000" free="0" />
       <parameter name="Normalization" value="1" scale="1" min="0" max="1000" free="0" />
       <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
       ...
       <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     </source>
     <source name="Background_vp0001_0_001886-004743keV" type="DRBPhibarBins" instrument="COM" id="vp0001_0_001886-004743keV">
       <parameter name="Normalization" value="1" scale="1" min="0" max="1000" free="0" />
       <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
       ...
       <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     </source>
     <source name="Background_vp0001_0_004743-011929keV" type="DRBPhibarBins" instrument="COM" id="vp0001_0_004743-011929keV">
       <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
       ...
       <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     </source>
     <source name="Background_vp0001_0_011929-029999keV" type="DRBPhibarBins" instrument="COM" id="vp0001_0_011929-029999keV">
       <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
       ...
       <parameter name="Normalization" value="1" error="0" scale="1" min="0" max="1000" free="1" />
     </source>
   </source_library>

.. note::

   To shorten the display, identical lines were replaced by three dots.

