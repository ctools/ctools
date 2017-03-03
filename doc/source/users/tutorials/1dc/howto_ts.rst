.. _1dc_howto_ts:

How to compute the signifiance of a source?
-------------------------------------------

To compute the significance of a given source you need to add the attribute
``tscalc="1"`` to the ``<source>`` tag in the
:ref:`model definition file <glossary_moddef>`
and execute :ref:`ctlike` using this file.
The structure of the XML file with significance computation requested for
``Src001`` is shown below.

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Src001" type="PointSource" tscalc="1">
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
     ...
   </source_library>

:ref:`ctlike` will compute the so-called Test Statistic, which is twice the
difference between the log-likelihood value obtained for a fit that includes
``Src001`` and a fit that excludes ``Src001``. The Test Statistic value will
follow a :math:`\chi^2_p` distribution with :math:`p` degrees of
freedom, where :math:`p` is the number of free parameters of ``Src001``.
For :math:`p=1`, the source significance is the square root of the Test
Statistic.

An excerpt of the :ref:`ctlike` log file is shown below where signifance
computation was requested for ``Src001`` and ``Src002`` and the diffuse
Galactic emission component ``IEM``. The computation has been done using
a stacked analysis.

.. code-block:: bash

   2017-03-03T23:36:07: === GModelSky ===
   2017-03-03T23:36:07:  Name ......................: Src001
   2017-03-03T23:36:07:  Instruments ...............: all
   2017-03-03T23:36:07:  Test Statistic ............: 12650.600597267
   ...
   2017-03-03T23:36:07: === GModelSky ===
   2017-03-03T23:36:07:  Name ......................: Src002
   2017-03-03T23:36:07:  Instruments ...............: all
   2017-03-03T23:36:07:  Test Statistic ............: 2154.31897658715
   ...
   2017-03-03T23:36:07: === GModelSky ===
   2017-03-03T23:36:07:  Name ......................: IEM
   2017-03-03T23:36:07:  Instruments ...............: all
   2017-03-03T23:36:07:  Test Statistic ............: 43070.3063685261

The Test Statistic values are also written into the
:ref:`model definition file <glossary_moddef>`
created by :ref:`ctlike`:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <source_library title="source library">
     <source name="Src001" type="PointSource" ts="12650.601" tscalc="1">
       ...
     </source>
     <source name="Src002" type="PointSource" ts="2154.319" tscalc="1">
       ...
     </source>
     <source name="IEM" type="DiffuseSource" ts="43070.306" tscalc="1">
       ...
     </source>
     ...
   </source_library>

