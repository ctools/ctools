.. _howto_comptel_ulimit:

Derive upper flux limits for a source
-------------------------------------

  .. admonition:: What you will learn

     You will learn how to use :ref:`ctulimit` to determine upper flux limits
     from COMPTEL data using a likelihood profile method.

Suppose we have tried to get the flux from a source, for example the Vela
pulsar, but we have not managed to detect the source and now want to determine
an upper flux limit for that object. You can do this with the :ref:`ctulimit`
tool.

First you have to retreive the COMPTEL data for the Vela pulsar. You can
do this by following
:ref:`the tutorial that explains how to download COMPTEL data <howto_comptel_download>`.
Suppose you have selected the data from viewing period 8. You then have to
create the following
:ref:`observation definition file <glossary_obsdef>`:

.. code-block:: xml

   <?xml version="1.0" standalone="no"?>
   <observation_list title="observation library">
     <observation name="Vela" id="100001" instrument="COM">
       <parameter name="DRE" file="m50866_dre.fits"/>
       <parameter name="DRB" file="m35330_drg.fits"/>
       <parameter name="DRG" file="m35330_drg.fits"/>
       <parameter name="DRX" file="m32504_drx.fits"/>
       <parameter name="IAQ" value="ENERG(0.75-1.0)MeV"/>
     </observation>
     <observation name="Vela" id="200001" instrument="COM">
       <parameter name="DRE" file="m50867_dre.fits"/>
       <parameter name="DRB" file="m35330_drg.fits"/>
       <parameter name="DRG" file="m35330_drg.fits"/>
       <parameter name="DRX" file="m32504_drx.fits"/>
       <parameter name="IAQ" value="ENERG(1.0-3.0)MeV"/>
     </observation>
     <observation name="Vela" id="300001" instrument="COM">
       <parameter name="DRE" file="m50868_dre.fits"/>
       <parameter name="DRB" file="m35330_drg.fits"/>
       <parameter name="DRG" file="m35330_drg.fits"/>
       <parameter name="DRX" file="m32504_drx.fits"/>
       <parameter name="IAQ" value="ENERG(3.0-10.0)MeV"/>
     </observation>
     <observation name="Vela" id="400001" instrument="COM">
       <parameter name="DRE" file="m50869_dre.fits"/>
       <parameter name="DRB" file="m35330_drg.fits"/>
       <parameter name="DRG" file="m35330_drg.fits"/>
       <parameter name="DRX" file="m32504_drx.fits"/>
       <parameter name="IAQ" value="ENERG(10.0-30.0)MeV"/>
     </observation>
   </observation_list>

As next step you have to create a
:ref:`model definition file <glossary_moddef>`.
You can copy the file used for the Crab and adjust the point source
component as shown below. Note that you need to fix the spectral index, and
you should also constrain the ``Prefactor`` to positive values.

.. code-block:: xml

   <?xml version="1.0" standalone="no"?>
   <source_library title="source library">
     <source name="Vela" type="PointSource">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor"   scale="1e-3" value="1.0"  min="0.0"  max="1000.0" free="1"/>
         <parameter name="Index"       scale="-1"   value="1.5"  min="0.0"  max="+5.0"   free="0"/>
         <parameter name="PivotEnergy" scale="1.0"  value="1.0"  min="0.01" max="1000.0" free="0"/>
       </spectrum>
       <spatialModel type="PointSource">
         <parameter name="RA"  scale="1.0" value="128.84" min="-360" max="360" free="0"/>
         <parameter name="DEC" scale="1.0" value="-45.18" min="-90"  max="90"  free="0"/>
       </spatialModel>
     </source>
     <source name="Background(0.75-1.0)MeV" type="DRBFitting" instrument="COM" id="100001">
     ...
   </source_library>

Now you can run the :ref:`ctulimit` tool as follows. Since the tool is normally
tuned for Cherenkov Telescope data analysis, energies are by default given in
TeV. To specify a reference energy of the differential flux limit of 1 MeV you
need to provide the hidden parameter ``eref=0.000001``. You can also set the
flux interval for the intergal flux computation to 0.75 - 30 MeV by
specifying ``emin=0.000000075 emax=0.00003`` (all values are given in TeV):

.. code-block:: bash

   $ ctulimit eref=0.000001 emin=0.000000075 emax=0.00003
   Input event list, counts cube or observation definition XML file [events.fits] obs.xml
   Source of interest [Crab] Vela
   Input model definition XML file [$CTOOLS/share/models/crab.xml] models.xml

The results of the upper limit computation can then be extracted from the
``ctulimit.log`` file that is created by :ref:`ctulimit`:

.. code-block:: bash

   2017-08-27T08:13:44: +=====================+
   2017-08-27T08:13:44: | Upper limit results |
   2017-08-27T08:13:44: +=====================+
   2017-08-27T08:13:44:  Differential flux limit ...: 9.26955448510129e-06 ph/cm2/s/MeV at 1e-06 TeV
   2017-08-27T08:13:44:  Integral flux limit .......: 6.43104918005e-05 ph/cm2/s within [7.5e-08-3e-05] TeV
   2017-08-27T08:13:44:  Energy flux limit .........: 1.54555137999306e-10 erg/cm2/s within [7.5e-08-3e-05] TeV
