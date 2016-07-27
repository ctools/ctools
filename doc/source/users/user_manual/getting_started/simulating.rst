.. _sec_simulating_cta:

Simulating CTA data
~~~~~~~~~~~~~~~~~~~

You can simulate a CTA observation with the :ref:`ctobssim` tool.

To invoke the tool, type :ref:`ctobssim` at the console prompt
(which is denoted by ``$``).
:ref:`ctobssim` will query for a number of parameters:

.. code-block:: bash

  $ ctobssim
  RA of pointing (degrees) (0-360) [83.63] 
  Dec of pointing (degrees) (-90-90) [22.01] 
  Radius of FOV (degrees) (0-180) [5.0] 
  Start time (MET in s) [0.0] 
  End time (MET in s) [1800.0] 
  Lower energy limit (TeV) [0.1] 
  Upper energy limit (TeV) [100.0] 
  Calibration database [prod2] 
  Instrument response function [South_0.5h] 
  Input model XML file [$CTOOLS/share/models/crab.xml] 
  Output event data file or observation definition XML file [events.fits]

Each line represents a query for one parameter value.
The line starts with a short description of the parameter, followed by 
the default value proposed by ctobssim in squared brackets ``[ ]``.

**If no parameter is entered the default value will be used**
(which is the case for all parameters shown here).
Otherwise, the specified value will overwrite the default value.
The round brackets ``( )`` indicate the range of possible parameter
values (if applicable).


You may have recognised that the environment variable ``$CTOOLS`` has 
been used in the path name of the model. ctools will automatically expand
the environment variables in parameter inputs.

The CTA instrument properties (effective area, PSF width) are taken from
the prod2 database. The response for the southern array using the cuts
optimised for 50 hours of observing time are used.

.. note::

   ctools comes bundled with CTA response functions for the northern and
   the southern array. The response functions are based on a Prod2
   analysis with cuts optimised for 0.5 hours, 5 hours and 50 hours of
   observing time. The following Instrument Response Functions
   (IRFs) are available: ``North_0.5h``, ``North_5h``, ``North_50h``,
   ``South_0.5h``, ``South_5h``, and ``South_50h``.

Events are simulated based on the instrument properties and based on a
source and background model. Only events that fall within the specified
region of interest (ROI), defined as a circle around a sky position in
Right Ascension and Declination (in degrees), will be stored in the output
event data file. The duration of the simulation is taken here to 30 minutes
(or 1800 seconds). Events are simulated for energies between 0.1 and 100 TeV.

The source and background model is defined by the XML file
``$CTOOLS/share/models/crab.xml``:

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <source_library title="source library">
    <source name="Crab" type="PointSource">
      <spectrum type="PowerLaw">
         <parameter name="Prefactor"   scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
         <parameter name="Index"       scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
         <parameter name="PivotEnergy" scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
      </spectrum>
      <spatialModel type="SkyDirFunction">
        <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="0"/>
        <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="0"/>
      </spatialModel>
    </source>
    <source name="CTABackgroundModel" type="CTAIrfBackground" instrument="CTA">
      <spectrum type="PowerLaw">
        <parameter name="Prefactor"   scale="1.0"  value="1.0"  min="1e-3" max="1e+3"   free="1"/>
        <parameter name="Index"       scale="1.0"  value="0.0"  min="-5.0" max="+5.0"   free="1"/>
        <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
      </spectrum>
    </source>
  </source_library>

The model consists of a source library that contains two components:
the Crab nebula and an instrumental background model.

The Crab nebula is modelled by a factorized sky model that has a spectral
and a spatial component (tags ``<spectrum>`` and ``<spatialModel>``,
respectively). The spectrum is modelled by a power law, which is defined by 
three parameters: the ``Prefactor``, the ``Index`` and the ``Scale``.
The spatial model has two parameters: Right Ascension in degrees (``RA``), and 
Declination in degrees (``DEC``). Each parameter has a value and a scale factor, 
the real value of the parameter being the product ``value * scale``. Typically,
``scale`` is chosen so that ``value`` is of the order of 1 (this is relevant for 
model fitting). In addition, ``value`` is bound by a minimum (``min``) and 
maximum (``max``) value, and a parameter may be free (``free="1"``) or fixed
(``free="0"``). The ``min``, ``max``, and ``free`` attributes are not
relevant here for the simulations, but they will be important for the model 
fitting later.

The spectral intensity I(E) (in units of photons/cm2/s/MeV) of the power law is
given by

.. math::
    \frac{dN}{dE} = N_0 \left( \frac{E}{E_0} \right)^{\gamma}

where the parameters in the XML definition have the following mappings:

* :math:`N_0` = ``Prefactor``
* :math:`\gamma` = ``Index``
* :math:`E_0` = ``PivotEnergy``

.. warning::

   Energies are given in the XML file in MeV units. This is a GammaLib
   convention that can not be modified. **So make sure you always use 
   MeV as energy unit in an XML file.**

The instrumental background of CTA is modelled using the background
information provided in the IRF (``CTAIrfBackground``) multipled
by a power law. As it is defined here, the power law represents a
constant of 1, hence the background IRF will be used without any
modification. The power law will become active when fitting the data
later and allows a spectral adjustment of the background model that
may account for uncertainties in the background information provided
in the IRF.

:ref:`ctobssim` has a couple of hidden parameters, the most important one being
certainly ``seed``. ``seed`` is an integer that specifies the seed value
for the random number generator, and changing this parameter will allow to
generate statistically independent Monte Carlo samples of CTA event data.
To use for example a seed value of 41 you should type:

.. code-block:: bash

  $ ctobssim seed=41

:ref:`ctobssim` will write 2 files in the working directory: ``events.fits``
and ``ctobssim.log``. The first file contains the simulated events in FITS 
format and can be inspected using ``fv`` or ``ds9``. The FITS file will 
contain three extensions: an empty primary image, a binary table named 
``EVENTS`` that holds the events (one row per event), and a binary table
named ``GTI`` holding the Good Time Intervals (for the moment a single row
with 2 columns providing the start and the stop time of the simulated time
interval).

The second file produced by :ref:`ctobssim` is a human readable log file that
contains information about the job execution. As example, the last lines
from this file are shown here:

.. code-block:: none

  2016-06-29T10:21:57: +======================+
  2016-06-29T10:21:57: | Simulate observation |
  2016-06-29T10:21:57: +======================+
  2016-06-29T10:21:57: === CTA observation ===
  2016-06-29T10:21:57:  Simulation cone ...........: RA=83.63 deg, Dec=22.01 deg, radius=5.5 deg
  2016-06-29T10:21:57:  Time interval .............: 0 - 1800 s
  2016-06-29T10:21:57:  Photon energy range .......: 100 GeV - 199.526 GeV
  2016-06-29T10:21:57:  Event energy range ........: 100 GeV - 199.526 GeV
  2016-06-29T10:21:57:   Simulation area ..........: 5.75561e+09 cm2
  2016-06-29T10:21:57:   Use model ................: Crab
  2016-06-29T10:21:57:   Normalization ............: 1 [Crab]
  2016-06-29T10:21:57:   Flux .....................: 3.76031e-10 [Crab] photons/cm2/s
  2016-06-29T10:21:57:   Normalized flux ..........: 3.76031e-10 [Crab] photons/cm2/s
  2016-06-29T10:21:57:   Photon rate ..............: 2.16429 photons/s [Crab]
  2016-06-29T10:21:57:   MC source photons ........: 3889 [Crab]
  2016-06-29T10:21:57:   MC source events .........: 1226 [Crab]
  2016-06-29T10:21:57:   MC source events .........: 1226 (all source models)
  2016-06-29T10:21:57:  Photon energy range .......: 199.526 GeV - 398.107 GeV
  ...
  2016-06-29T10:21:57:  MC source photons .........: 10759 [Crab]
  2016-06-29T10:21:57:  MC source events ..........: 3686 [Crab]
  2016-06-29T10:21:57:  MC events outside ROI .....: 0
  2016-06-29T10:21:57:  MC background events ......: 19413
  2016-06-29T10:21:57:  MC events .................: 23099 (all models)

Each line starts with the UTC time at which the line has been written. In
this run, 10759 Crab photons have been thrown during a time interval of 1800
seconds. 3686 of these photons have been registered by CTA as events. In the
same time interval, 19413 background events have been registred by CTA.

.. note::

   :ref:`ctobssim` will split the simulated energy range into a number of
   slices, controlled via the hidden ``eslices`` parameter (ten energy slices
   are used by default). For each energy slice, the simulation area
   will be adapted to the effective area of the array in that energy slice,
   which helps to keep the computing time low. The log file will provide
   information about the simulation in each slice. In the example above, the
   simulation results for the first energy slice are shown, followed by a
   summary of the results for all slices.

You may change the name of the log file using the hidden parameter 
``logfile``:

.. code-block:: bash

  $ ctobssim logfile=my-private-log-file

Furthermore, you may decide on the amount of information provided in the 
log file (the chattiness of the executable) using the hidden parameter 
``chatter``:

.. code-block:: bash

  $ ctobssim chatter=4

``chatter`` can vary between 0 and 4, 0 providing no information while 4 
provides the most detailed information.

**By default, all ctools have a chatter level of 2.**

You may also duplicate the log file information into the console by setting
the hidden ``debug`` parameter to yes:

.. code-block:: bash

  $ ctobssim debug=yes
