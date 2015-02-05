.. _sec_connecting_irf:

Connecting observations to dedicated instrument response functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   To execute the following examples on your machine you need to
   install a specific calibration database. Download the file
   :download:`caldb.tar.gz <caldb.tar.gz>` into a local working
   directory, decompress and untar the folder, and set the calibration
   database to this folder by typing:

   .. code-block:: bash

      $ export CALDB=$PWD/caldb

   Note that the response functions in the calibration database
   are all identical to the ``cta_dummy_irf`` response, yet there
   a multiple database and response entries so that the usage of
   specific response functions can be emulated.

So far the same instrument response function has been used for all
observations in our analysis examples.
You may have recognised that we always specified the calibration database
``dummy`` and the instrument response function ``cta_dummy_irf`` that
are shipped with the ctools release.
When combining multiple observations in a joint analysis we may however
need to specify dedicated instrument response functions for the
different observations.
This can be done by adding the calibration information to the observation
definition XML file:

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <observation_list title="observation library">
    <observation name="Crab" id="00001" instrument="CTA">
      <parameter name="EventList" file="events1.fits"/>
      <parameter name="Calibration" database="north" response="cta_north_irf"/>
    </observation>
    <observation name="Crab" id="00002" instrument="CTA">
      <parameter name="EventList" file="events2.fits"/>
      <parameter name="Calibration" database="south" response="cta_south_irf"/>
    </observation>
  </observation_list>

Each observation now has a ``Calibration`` parameter with the attributes
``database`` and ``response``.
To see which calibration databases and response functions are available on
your system you can run the :ref:`cscaldb` script (use the ``debug=yes``
option to get the output of the script on the console):

.. code-block:: bash

   $ cscaldb debug=yes
   Parfile cscaldb.par not found. Create default parfile.
   2015-02-05T14:26:50: +============+
   2015-02-05T14:26:50: | Parameters |
   2015-02-05T14:26:50: +============+
   2015-02-05T14:26:50:  chatter ...................: 2
   2015-02-05T14:26:50:  clobber ...................: yes
   2015-02-05T14:26:50:  debug .....................: yes
   2015-02-05T14:26:50:  mode ......................: ql
   2015-02-05T14:26:50: 
   2015-02-05T14:26:50: +==============+
   2015-02-05T14:26:50: | Mission: cta |
   2015-02-05T14:26:50: +==============+
   2015-02-05T14:26:50: === Response functions in database "dummy" ===
   2015-02-05T14:26:50: cta_dummy_irf
   2015-02-05T14:26:50: cta_prod1_irf
   2015-02-05T14:26:50: cta_prod2_irf
   2015-02-05T14:26:50: 
   2015-02-05T14:26:50: === Response functions in database "north" ===
   2015-02-05T14:26:50: cta_north_irf
   2015-02-05T14:26:50: 
   2015-02-05T14:26:50: === Response functions in database "south" ===
   2015-02-05T14:26:50: cta_south_irf

You will see that for the CTA observatory three databases exist.
The ``dummy`` database contains three response functions
(``cta_dummy_irf``, ``cta_prod1_irf``, ``cta_prod2_irf``) while the
``north`` and ``south`` databases contain a single response function.
You may now do a :ref:`ctlike` fit using the XML file with the specific
response definitions:

.. code-block:: bash

   $ ctlike
   Event list, counts cube or observation definition file [events.fits] obs_irf.xml
   Source model [$CTOOLS/share/models/crab.xml] 
   Source model output file [crab_results.xml]

To see the usage of the individual response functions you may inspect the
log file. You will notice that each observation now has a specific response
function and that the filename of the response information differs for
both observations.

.. code-block:: xml

   2015-02-05T14:48:46: +==============+
   2015-02-05T14:48:46: | Observations |
   2015-02-05T14:48:46: +==============+
   2015-02-05T14:48:46: === GObservations ===
   2015-02-05T14:48:46:  Number of observations ....: 2
   2015-02-05T14:48:46:  Number of predicted events : 0
   2015-02-05T14:48:46: === GCTAObservation ===
   2015-02-05T14:48:46:  Name ......................: Crab
   2015-02-05T14:48:46:  Identifier ................: 00001
   2015-02-05T14:48:46:  Instrument ................: CTA
   2015-02-05T14:48:46:  Event file ................: events1.fits
   2015-02-05T14:48:46:  Event type ................: EventList
   2015-02-05T14:48:46:  Statistics ................: Poisson
   2015-02-05T14:48:46:  Ontime ....................: 1800 s
   2015-02-05T14:48:46:  Livetime ..................: 1710 s
   2015-02-05T14:48:46:  Deadtime correction .......: 0.95
   2015-02-05T14:48:46:  User energy range .........: undefined
   2015-02-05T14:48:46: === GCTAPointing ===
   2015-02-05T14:48:46:  Pointing direction ........: (RA,Dec)=(83.63,21.51)
   2015-02-05T14:48:46: === GCTAResponseIrf ===
   2015-02-05T14:48:46:  Response name .............: cta_north_irf
   2015-02-05T14:48:46:  Energy dispersion .........: Not used
   2015-02-05T14:48:46:  Save energy range .........: undefined
   2015-02-05T14:48:46: === GCaldb ===
   2015-02-05T14:48:46:  Database root .............: /Users/jurgen/Desktop/tmp/beyond/caldb
   2015-02-05T14:48:46:  Selected Mission ..........: CTA
   2015-02-05T14:48:46:  Selected Instrument .......: NORTH
   2015-02-05T14:48:46:  Calibration Index File ....: /Users/jurgen/Desktop/tmp/beyond/caldb/data/cta/north/caldb.indx
   2015-02-05T14:48:46:  Number of entries .........: 4
   2015-02-05T14:48:46: === GCTAAeffPerfTable ===
   2015-02-05T14:48:46:  Filename ..................: /Users/jurgen/Desktop/tmp/beyond/caldb/data/cta/north/bcf/cta_north_irf.dat
   2015-02-05T14:48:46:  Number of energy bins .....: 20
   2015-02-05T14:48:46:  Log10(Energy) range .......: 0.0199526 - 125.893 TeV
   2015-02-05T14:48:46:  Offset angle dependence ...: Fixed sigma=3
   ...
   2015-02-05T14:48:46: === GCTAObservation ===
   2015-02-05T14:48:46:  Name ......................: Crab
   2015-02-05T14:48:46:  Identifier ................: 00002
   2015-02-05T14:48:46:  Instrument ................: CTA
   2015-02-05T14:48:46:  Event file ................: events2.fits
   2015-02-05T14:48:46:  Event type ................: EventList
   2015-02-05T14:48:46:  Statistics ................: Poisson
   2015-02-05T14:48:46:  Ontime ....................: 1800 s
   2015-02-05T14:48:46:  Livetime ..................: 1710 s
   2015-02-05T14:48:46:  Deadtime correction .......: 0.95
   2015-02-05T14:48:46:  User energy range .........: undefined
   2015-02-05T14:48:46: === GCTAPointing ===
   2015-02-05T14:48:46:  Pointing direction ........: (RA,Dec)=(83.63,22.51)
   2015-02-05T14:48:46: === GCTAResponseIrf ===
   2015-02-05T14:48:46:  Response name .............: cta_south_irf
   2015-02-05T14:48:46:  Energy dispersion .........: Not used
   2015-02-05T14:48:46:  Save energy range .........: undefined
   2015-02-05T14:48:46: === GCaldb ===
   2015-02-05T14:48:46:  Database root .............: /Users/jurgen/Desktop/tmp/beyond/caldb
   2015-02-05T14:48:46:  Selected Mission ..........: CTA
   2015-02-05T14:48:46:  Selected Instrument .......: SOUTH
   2015-02-05T14:48:46:  Calibration Index File ....: /Users/jurgen/Desktop/tmp/beyond/caldb/data/cta/south/caldb.indx
   2015-02-05T14:48:46:  Number of entries .........: 4
   2015-02-05T14:48:46: === GCTAAeffPerfTable ===
   2015-02-05T14:48:46:  Filename ..................: /Users/jurgen/Desktop/tmp/beyond/caldb/data/cta/south/bcf/cta_south_irf.dat
   2015-02-05T14:48:46:  Number of energy bins .....: 20
   2015-02-05T14:48:46:  Log10(Energy) range .......: 0.0199526 - 125.893 TeV
   2015-02-05T14:48:46:  Offset angle dependence ...: Fixed sigma=3

A more fine grained control over the response function can be achieved by 
specifying individual filenames for the various response components.
An example for an observation definition XML file is shown below.
This is definitely expert mode, to be used with utmost care.

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <observation_list title="observation library">
    <observation name="Crab" id="00001" instrument="CTA">
      <parameter name="EventList"           file="events1.fits"/>
      <parameter name="EffectiveArea"       file="$(CALDB)/data/cta/dummy/bcf/cta_dummy_irf.dat"/>
      <parameter name="PointSpreadFunction" file="$(CALDB)/data/cta/dummy/bcf/cta_prod1_irf.dat"/>
      <parameter name="EnergyDispersion"    file="$(CALDB)/data/cta/dummy/bcf/cta_prod2_irf.dat"/>
      <parameter name="Background"          file="$(CALDB)/data/cta/north/bcf/cta_north_irf.dat"/>
    </observation>
    <observation name="Crab" id="00002" instrument="CTA">
      <parameter name="EventList"           file="events2.fits"/>
      <parameter name="EffectiveArea"       file="$(CALDB)/data/cta/dummy/bcf/cta_dummy_irf.dat"/>
      <parameter name="PointSpreadFunction" file="$(CALDB)/data/cta/dummy/bcf/cta_prod1_irf.dat"/>
      <parameter name="EnergyDispersion"    file="$(CALDB)/data/cta/dummy/bcf/cta_prod2_irf.dat"/>
      <parameter name="Background"          file="$(CALDB)/data/cta/south/bcf/cta_south_irf.dat"/>
    </observation>
  </observation_list>

Finally, response information may also be provided to combine stacked
observations. An example for the syntax of the observation definition XML 
file is given below:

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <observation_list title="observation library">
    <observation name="Crab" id="00001" instrument="CTA">
      <parameter name="CountsCube"   file="cntcube1.fits"/>
      <parameter name="ExposureCube" file="expcube1.fits"/>
      <parameter name="PsfCube"      file="psfcube1.fits"/>
    </observation>
    <observation name="Crab" id="00002" instrument="CTA">
      <parameter name="CountsCube"   file="cntcube2.fits"/>
      <parameter name="ExposureCube" file="expcube2.fits"/>
      <parameter name="PsfCube"      file="psfcube2.fits"/>
    </observation>
  </observation_list>


