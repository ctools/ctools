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
source ``Src1``. Your
:ref:`model definition XML file <glossary_moddef>`
should then look like this:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Src001" type="PointSource" tscalc="1">
       <spectrum type="ExponentialCutoffPowerLaw">
         <parameter name="Prefactor"    scale="1e-18" value="5.7"  min="1e-07" max="1000.0" free="1"/>
         <parameter name="Index"        scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
         <parameter name="CutoffEnergy" scale="1e7"   value="1.0"  min="0.01"  max="1000.0" free="1"/>
         <parameter name="PivotEnergy"  scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA"  scale="1" value="266.424" free="1" />
         <parameter name="DEC" scale="1" value="-29.005" free="1" />
       </spatialModel>
     </source>
     ...
   </source_library>

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
   2018-01-25T12:46:12:  Test Statistic ............: 7659.64342272468

The Test Statistic values are also written into the
:ref:`model definition file <glossary_moddef>`
created by :ref:`ctlike`:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Src001" type="PointSource" ts="7659.643" tscalc="1">
       ...
     </source>
     ...
   </source_library>
