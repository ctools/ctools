.. _howto_ulimit:

How to compute upper limits?
----------------------------

  .. admonition:: What you will learn

     You will learn how to **determine the upper flux limit for a source**
     using the :ref:`ctulimit` tool.

You determine upper limits for gamma-ray sources using the :ref:`ctulimit`
tool.

Suppose you want to know the upper limit on the gamma-ray flux for
a putative source at position :math:`l=0.1^\circ` and :math:`b=0.4^\circ`.
For this you have to add a test source with the respective coordinates to
the
:ref:`model definition file <glossary_moddef>`.
To do this, copy the 1DC :ref:`model definition <glossary_moddef>` file

.. code-block:: bash

   $ cp $CTOOLS/share/models/1dc_howto.xml ulimit_model.xml

to :download:`ulimit_model.xml <ulimit_model.xml>` and add the following
source to it:

.. literalinclude:: ulimit_model.xml
   :language: xml
   :lines: 3-13

Now run the :ref:`ctulimit` tool as follows:

.. code-block:: bash

   $ ctulimit
   Input event list, counts cube or observation definition XML file [events.fits] cntcube.fits
   Input exposure cube file [NONE] expcube.fits
   Input PSF cube file [NONE] psfcube.fits
   Input background cube file [NONE] bkgcube.fits
   Source of interest [Crab] Test
   Input model definition XML file [$CTOOLS/share/models/crab.xml] ulimit_model.xml

After the run, the upper gamma-ray flux limit for this test source can be
extracted from the :ref:`ctulimit` log file:

.. code-block:: bash

   2019-04-04T10:31:20: +=====================+
   2019-04-04T10:31:20: | Upper limit results |
   2019-04-04T10:31:20: +=====================+
   2019-04-04T10:31:20:  Differential flux limit ...: 1.3961626479978e-20 ph/cm2/s/MeV at 1 TeV
   2019-04-04T10:31:20:  Integral flux limit .......: 9.42318774389192e-15 ph/cm2/s within [1-100] TeV
   2019-04-04T10:31:20:  Energy flux limit .........: 4.14922476753972e-14 erg/cm2/s within [1-100] TeV

.. warning::
   It is important to fix the spectral index of the test source in the
   upper limit computation, since for a source that is not detected, the
   index of the spectral law cannot be constrained. In other words, you have
   to make an assumption about the spectral energy density distribution when
   computing the upper flux limit of a source. The spectral index is fixed
   by setting ``free="0"`` for the ``Index`` parameter.
