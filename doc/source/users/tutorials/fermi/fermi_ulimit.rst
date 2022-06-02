.. _fermi_ulimit:

Derive upper flux limits
------------------------

  .. admonition:: What you will learn

     You will learn how to use :ref:`ctulimit` to determine upper flux limits
     from Fermi-LAT data using a likelihood profile method.

Suppose there is another potential source in the field of view but you cannot
see this source in a counts map, hence you want to derive an upper flux limit
for the source. To do this you first have to add an additional source to the
:ref:`model definition file <glossary_moddef>`. An example for a source at
:math:`\alpha=122^\circ` and :math:`\delta=-42^\circ` is shown below:

.. code-block:: xml

   <?xml version="1.0" standalone="no"?>
   <source_library title="source library">
     <source type="PointSource" name="Test">
       <spectrum type="PowerLaw">
         <parameter name="PhotonFlux" scale="1e-12" value="1.0"      min="1e-07" max="1e+07"     free="1"/>
         <parameter name="Index"      scale="1.0"   value="-2.0"     min="-5.0"  max="+5.0"      free="0"/>
         <parameter name="LowerLimit" scale="1.0"   value="100.0"    min="10.0"  max="1000000.0" free="0"/>
         <parameter name="UpperLimit" scale="1.0"   value="100000.0" min="10.0"  max="1000000.0" free="0"/>
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA"  scale="1.0" value="122.00" min="-360" max="360" free="0"/>
         <parameter name="DEC" scale="1.0" value="-42.00" min="-90"  max="90"  free="0"/>
       </spatialModel>
     </source>
     <source type="PointSource" name="Vela">
     ...
   </source_library>

Now you can run the :ref:`ctulimit` tool. Since the tool is normally tuned
for Cherenkov Telescope data analysis, energies are by default given in TeV.
To specify a reference energy of the differential flux limit of 1 GeV you
need to provide the hidden parameter ``eref=0.001``. You can also set the
flux interval for the intergal flux computation to 100 MeV - 100 GeV by
specifying ``emin=0.0001 emax=0.1`` (all values are given in TeV):

.. code-block:: bash

   $ ctulimit eref=0.001 emin=0.0001 emax=0.1
   Input event list, counts cube or observation definition XML file [events.fits] obs.xml
   Source of interest [Crab] Test
   Input model definition XML file [$CTOOLS/share/models/crab.xml] models_ulimit.xml

The results of the upper limit computation can then be extracted from the
``ctulimit.log`` file that is created by :ref:`ctulimit`:

.. code-block:: bash

   2019-04-04T15:03:22: +=====================+
   2019-04-04T15:03:22: | Upper limit results |
   2019-04-04T15:03:22: +=====================+
   2019-04-04T15:03:22:  Differential flux limit ...: 1.77270316013682e-13 ph/cm2/s/MeV at 0.001 TeV
   2019-04-04T15:03:22:  Integral flux limit .......: 1.77093045697669e-09 ph/cm2/s within [0.0001-0.1] TeV
   2019-04-04T15:03:22:  Energy flux limit .........: 1.96192914922622e-12 erg/cm2/s within [0.0001-0.1] TeV
