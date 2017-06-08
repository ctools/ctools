.. _1dc_howto_ulimit:

How to compute upper limits?
----------------------------

  .. admonition:: What you will learn

     You will learn how to **determine the upper flux limit for a source**
     using the :ref:`ctulimit` tool.

You determine upper limits for gamma-ray sources using the :ref:`ctulimit`
tool.

Suppose you want to know the upper limit on the gamma-ray flux for
a putative source at position :math:`l=0.1^\circ` and :math:`b=0.4^\circ`
(corresponding to :math:`\alpha=266.0753^\circ` and
:math:`\delta=-28.6420^\circ`).
For this you have to add a test source with the respective coordinates to
the
:ref:`model definition file <glossary_moddef>`,
for example

.. code-block:: xml

   <source name="Test" type="PointSource">
     <spectrum type="PowerLaw">
       <parameter name="Prefactor"   value="1" scale="5.7e-18" min="0"  free="1" />
       <parameter name="Index"       value="1" scale="-2.48"   min="-5" max="5" free="0" />
       <parameter name="PivotEnergy" value="1" scale="300000" free="0" />
     </spectrum>
     <spatialModel type="PointSource">
       <parameter name="RA"  value="266.0753" scale="1" free="0" />
       <parameter name="DEC" value="-28.6420" scale="1" free="0" />
     </spatialModel>
   </source>

and then run the :ref:`ctulimit` tool as follows:

.. code-block:: bash

   $ ctulimit
   Input event list, counts cube or observation definition XML file [events.fits] cntcube.fits
   Input exposure cube file (only needed for stacked analysis) [NONE] expcube.fits
   Input PSF cube file (only needed for stacked analysis) [NONE] psfcube.fits
   Input background cube file (only needed for stacked analysis) [NONE] bkgcube.fits
   Source of interest [Crab] Test
   Input model definition XML file [$CTOOLS/share/models/crab.xml] stacked_models_cutoff_iem_up.xml

After the run, the upper gamma-ray flux limit for this test source can be
extracted from the :ref:`ctulimit` log file:

.. code-block:: bash

   2017-06-08T15:39:54: +=====================+
   2017-06-08T15:39:54: | Upper limit results |
   2017-06-08T15:39:54: +=====================+
   2017-06-08T15:39:54:  Differential flux limit ...: 2.4438738640373e-20 ph/cm2/s/MeV at 1 TeV
   2017-06-08T15:39:54:  Integral flux limit .......: 1.64945554704815e-14 ph/cm2/s within [1-100] TeV
   2017-06-08T15:39:54:  Energy flux limit .........: 7.26289446286815e-14 erg/cm2/s within [1-100] TeV

.. warning::
   It is important to fix the spectral index of the test source in the
   upper limit computation, since for a source that is not detected, the
   index of the spectral law cannot be constrained. In other words, you have
   to make an assumption about the spectral energy density distribution when
   computing the upper flux limit of a source. The spectral index is fixed
   by setting ``free="0"`` for the ``Index`` parameter.
