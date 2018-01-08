.. _start_simulating:

Simulating event data
---------------------

  .. admonition:: What you will learn

     You will learn how to use the :ref:`ctobssim` tool from the command
     line to **simulate event data**.

     :ref:`ctobssim` draws random events using a model of the celestial
     gamma-ray intensity distribution that was convolved with the
     :ref:`instrument response function <glossary_irf>`.
     The tool adds random events from the expected distribution of the
     residual background to the data.

Let's assume that you want to generate a simulated CTA event list. You
do this by using the :ref:`ctobssim` tool from the command line.

To invoke the tool, type :ref:`ctobssim` at the console prompt
(which is denoted by ``$``).
:ref:`ctobssim` will then query for several parameters that define the
characteristics of the simulation:

.. code-block:: bash

   $ ctobssim
   RA of pointing (degrees) (0-360) [83.63] 83.5
   Dec of pointing (degrees) (-90-90) [22.01] 22.8
   Radius of FOV (degrees) (0-180) [5.0]
   Start time (UTC string, JD, MJD or MET in seconds) [2020-01-01T00:00:00]
   Stop time (UTC string, JD, MJD or MET in seconds) [2020-01-01T01:00:00] 2020-01-01T01:00:00
   Lower energy limit (TeV) [0.03]
   Upper energy limit (TeV) [200.]
   Calibration database [prod2]
   Instrument response function [South_0.5h]
   Input model definition XML file [$CTOOLS/share/models/crab.xml]
   Output event data file or observation definition XML file [events.fits]

Each line represents a query for one parameter value.
The line starts with a short description of the parameter, followed by 
the default value proposed by :ref:`ctobssim` in squared brackets ``[ ]``.

**If no parameter value is entered the default value will be used**.
Otherwise, the specified value will overwrite the default value.
The round brackets ``( )`` indicate the range of possible parameter
values (if applicable).

  .. note::

     Times can be entered in various formats in ctools. Times can be provided
     as

     * UTC date strings (e.g. ``"2020-01-01T00:00:00"``)
     * Modified Julian days (e.g. ``"MJD 58849.0"``)
     * Julian days (e.g. ``"JD 2458849.5"``)
     * Mission elapsed time in seconds, counted from ``2000-01-01T12:00:00``
       (e.g. ``"631108869.18"``)

     Usually times are given in Terrestial Time (``TT``) but may also be
     specified as ``TAI`` or ``UTC`` (e.g. ``"MJD 58849.0 (TAI)"`` or
     ``"MJD 58849.0 (UTC)"``).

You may have recognised that the environment variable ``$CTOOLS`` has 
been used in the path name of the model. ctools will automatically expand
the environment variables in parameter values.

The CTA
:ref:`instrument response function <glossary_irf>` (IRF)
is taken from the ``prod2`` database. The response for the southern array
using the cuts optimised for 0.5 hours of observing time is used.

  .. note::

     ctools comes bundled with CTA
     :ref:`instrument response functions <glossary_irf>` for the northern and
     the southern array. The IRFs are based on a ``prod2``
     analysis with cuts optimised for 0.5 hours, 5 hours and 50 hours of
     observing time. **Be aware that these times do not need to correspond
     to the actual observing time that is simulated.** IRFs optimised for
     short observation times correspond to ``loose`` cuts that keep more
     photons but also more background events, while IRFs optimised for
     long observation times correspond to ``hard`` cuts that limit the
     number of background events at the expense of loosing some photons.
     The following IRFs are available:
     ``North_0.5h``, ``North_5h``, ``North_50h``,
     ``South_0.5h``, ``South_5h``, and ``South_50h``.

Events are simulated based on the
:ref:`instrument response function <glossary_irf>`
and based on a source and background model. Only events that fall within the
specified region of interest (ROI), defined as a circle around a sky position in
Right Ascension and Declination (in degrees), will be stored in the output
event data file. The duration of the simulation is taken here to one hour.
Events are simulated for energies between 30 GeV and 200 TeV.

The source and background model is defined by the
:ref:`model definition XML file <glossary_moddef>`
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
      <spatialModel type="PointSource">
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

The spectral intensity I(E) (in units of
:math:`{\rm photons} \, {\rm cm}^{-2} \, {\rm s}^{-1} \, {\rm MeV}^{-1}`)
of the power law is given by


.. math::
    \frac{dN}{dE} = N_0 \left( \frac{E}{E_0} \right)^{\gamma}

where the parameters in the XML definition have the following mappings:

* :math:`N_0` = ``Prefactor``
* :math:`\gamma` = ``Index``
* :math:`E_0` = ``PivotEnergy``

..

  .. warning::
     **Energies are given in the XML file in MeV units.** This is a GammaLib
     convention that can not be modified. So make sure you always use
     MeV as energy unit in an XML file.

 .. note::
    As customary for IACT observations, the pointing direction
    was slightly offset from the source of interest, i.e.,
    the Crab. This makes it possible to better handle systematics due
    to the instrumental background.

The instrumental background of CTA is modelled using the background
information provided in the
:ref:`instrument response function <glossary_irf>`
(``CTAIrfBackground``), where the energy dependence of the background
model is multipled by a power law. As it is defined here, the power law
represents a constant of 1, hence the background IRF will be used without any
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

..

  .. note::

     Hidden parameters are parameters that are not queried by a tool since
     in general their values is not expected to change frequently. To change
     hidden parameters they have to be given as arguments on the command line.
     Multiple hidden parameters need to be separated by a white space.

:ref:`ctobssim` will write two files in the working directory: ``events.fits``
and ``ctobssim.log``. The first file contains the simulated events in FITS 
format and can be inspected using ``fv`` or ``ds9``. The FITS file will 
contain three extensions: an empty primary image, a binary table named 
``EVENTS`` that holds the events (one row per event), and a binary table
named ``GTI`` holding the Good Time Intervals (for the moment a single row
with two columns providing the start and the stop time of the simulated time
interval).

The second file produced by :ref:`ctobssim` is a human readable log file that
contains information about the job execution. As example, the last lines
from this file are shown here:

.. code-block:: none

   2017-11-28T14:13:40: === CTA observation ===
   2017-11-28T14:13:40:  Simulation cone ...........: RA=83.5 deg, Dec=22.8 deg, radius=5.5 deg
   2017-11-28T14:13:40:  Time interval .............: 6.31109e+08 - 6.31112e+08 s
   2017-11-28T14:13:40:  Photon energy range .......: 30 GeV - 72.3622611060088 GeV
   2017-11-28T14:13:40:  Event energy range ........: 30 GeV - 72.3622611060088 GeV
   2017-11-28T14:13:40:   Simulation area ..........: 1.97769e+09 cm2
   2017-11-28T14:13:40:   Use model ................: Crab
   2017-11-28T14:13:40:   Normalization ............: 1 [Crab]
   2017-11-28T14:13:40:   Flux .....................: 2.5413e-09 [Crab] photons/cm2/s
   2017-11-28T14:13:40:   Normalized flux ..........: 2.5413e-09 [Crab] photons/cm2/s
   2017-11-28T14:13:40:   Photon rate ..............: 5.0259 photons/s [Crab]
   2017-11-28T14:13:40:   MC source photons ........: 18186 [Crab]
   2017-11-28T14:13:40:   MC source events .........: 3544 [Crab]
   2017-11-28T14:13:40:   MC source events .........: 3544 (all source models)
   2017-11-28T14:13:40:  Photon energy range .......: 72.3622611060088 GeV - 174.543227745807 GeV
   ...
   2017-11-28T14:13:40:  MC source photons .........: 47219 [Crab]
   2017-11-28T14:13:40:  MC source events ..........: 11356 [Crab]
   2017-11-28T14:13:49:  MC events outside ROI .....: 0
   2017-11-28T14:13:49:  MC background events ......: 189477
   2017-11-28T14:13:49:  MC identifier 1 ...........: Crab
   2017-11-28T14:13:49:  MC identifier 2 ...........: CTABackgroundModel
   2017-11-28T14:13:49:  MC events .................: 200833 (all models)

Each line starts with the UTC time at which the line has been written. In
this run, 47219 Crab photons have been thrown. 11356 of these photons have been
registered by CTA as events. In the same time interval, 189477 background
events have been registred by CTA.

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

..

  .. note::

     All tools have the hidden parameters ``logfile``, ``chatter``, and
     ``debug`` and you can use these parameters to control the log file
     output. In addition, all tools have the hidden parameter ``clobber``
     that allows to overwrite existing files (set to ``yes`` by default)
     and ``mode`` that defines the mode of automatic parameters (set to
     ``ql`` for *query and learn* by default).

