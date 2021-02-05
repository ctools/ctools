.. _howto_ts:

How to compute the signifiance of a source?
-------------------------------------------

  .. admonition:: What you will learn

     You will learn how to **compute the significance of a source detection**
     using the :ref:`ctlike` tool.

To compute the detection significance of a given source you need to add the
attribute ``tscalc="1"`` to the ``<source>`` tag in the
:ref:`model definition file <glossary_moddef>`
and execute :ref:`ctlike` using the modified file.
Copy the 1DC :ref:`model definition <glossary_moddef>` file

.. code-block:: bash

   $ cp $CTOOLS/share/models/1dc_howto.xml ts_model.xml

and modify the file by adding ``tscalc="1"`` to the ``<source>`` tag for
source ``Src1``. The first lines of the :download:`ts_model.xml <ts_model.xml>`
:ref:`model definition XML file <glossary_moddef>`
should then look like this:

.. literalinclude:: ts_model.xml
   :language: xml
   :lines: 1-14

Then run :ref:`ctlike` to compute the so-called Test Statistic, which is twice
the difference between the log-likelihood value obtained for a fit that includes
``Src001`` and a fit that excludes ``Src001``. The Test Statistic value will
follow a :math:`\chi^2_p` distribution with :math:`p` degrees of
freedom, where :math:`p` is the number of free parameters of ``Src001``.
For :math:`p=1`, the source significance is the square root of the Test
Statistic.

.. code-block:: bash

   $ ctlike
   Input event list, counts cube or observation definition XML file [events.fits] cntcube.fits
   Input exposure cube file [NONE] expcube.fits
   Input PSF cube file [NONE] psfcube.fits
   Input background cube file [NONE] bkgcube.fits
   Input model definition XML file [$CTOOLS/share/models/crab.xml] ts_model.xml
   Output model definition XML file [crab_results.xml] result.xml

An excerpt of the :ref:`ctlike` log file is shown below

.. code-block:: none

   2018-01-25T12:46:12: === GModelSky ===
   2018-01-25T12:46:12:  Name ......................: Src001
   2018-01-25T12:46:12:  Instruments ...............: all
   2018-01-25T12:46:12:  Test Statistic ............: 7659.72844228207

The Test Statistic values are also written into the
:ref:`model definition file <glossary_moddef>`
created by :ref:`ctlike`:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Src001" type="PointSource" ts="7659.728" tscalc="1">
       ...
     </source>
     ...
   </source_library>
