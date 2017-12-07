.. _1dc_howto_ts:

How to compute the signifiance of a source?
-------------------------------------------

  .. admonition:: What you will learn

     You will learn how to **compute the significance of a source detection**
     using the :ref:`ctlike` tool.

To compute the detection significance of a given source you need to add the
attribute ``tscalc="1"`` to the ``<source>`` tag in the
:ref:`model definition file <glossary_moddef>`
and execute :ref:`ctlike` using the modified file.
The structure of the XML file with significance computation requested for
``Src001`` is shown below.

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
         <parameter name="RA" value="266.424004498437" error="0" scale="1" free="1" />
         <parameter name="DEC" value="-29.0049010253548" error="0" scale="1" free="1" />
       </spatialModel>
     </source>
     ...
   </source_library>

:ref:`ctlike` will compute the so-called Test Statistic, which is twice the
difference between the log-likelihood value obtained for a fit that includes
``Src001`` and a fit that excludes ``Src001``. The Test Statistic value will
follow a :math:`\chi^2_p` distribution with :math:`p` degrees of
freedom, where :math:`p` is the number of free parameters of ``Src001``.
For :math:`p=1`, the source significance is the square root of the Test
Statistic.

An excerpt of the :ref:`ctlike` log file is shown below where significance
computation was requested for ``Src001`` and ``Src002``. The computation was
done using a stacked analysis.

.. code-block:: bash

   2017-07-28T01:26:56: === GModelSky ===
   2017-07-28T01:26:56:  Name ......................: Src001
   2017-07-28T01:26:56:  Instruments ...............: all
   2017-07-28T01:26:56:  Test Statistic ............: 15404.1562257314
   ...
   2017-07-28T01:26:56: === GModelSky ===
   2017-07-28T01:26:56:  Name ......................: Src002
   2017-07-28T01:26:56:  Instruments ...............: all
   2017-07-28T01:26:56:  Test Statistic ............: 3419.3020299729

The Test Statistic values are also written into the
:ref:`model definition file <glossary_moddef>`
created by :ref:`ctlike`:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Src001" type="PointSource" ts="15404.156" tscalc="1">
       ...
     </source>
     <source name="Src002" type="PointSource" ts="3419.302" tscalc="1">
       ...
     </source>
     ...
   </source_library>

