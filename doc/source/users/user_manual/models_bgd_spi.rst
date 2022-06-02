.. _um_models_bgd_spi:

SPI background models
---------------------

The instrumental background in INTEGRAL/SPI data is modelled by so-called
templates that describe the variation of the background rates with time:

.. math::
   M(p',E',t') = \sum_{i=0}^{N-1} n_i(d',E',t') \times {\rm TPL}_i(d',E',t')

where :math:`M(p',E',t')` is given in units of
:math:`{\rm events} \,\, {\rm s}^{-1} {\rm MeV}^{-1} {\rm sr}^{-1}`.
:math:`d'` is the SPI detector identifier,
:math:`E'` is the measured energy, and
:math:`t'` is the time that is linked via the pointing of the SPI
telescope to an area of the sky.

:math:`{\rm TPL}_i(d',E',t')` is the :math:`i` th out of :math:`N`
templates, and :math:`n_i(d',E',t')` is a normalisation factor for this
template that may depend on the SPI detector identifier, measured energy and
time.


Creating a background model
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since the number of background normalisation factors for INTEGRAL/SPI data
analysis may be large, it is not recommenced to create a
:ref:`model definition XML file <glossary_moddef>`
directly, but to use the ``gammalib.GSPIModelDataSpace()`` constructor in
a Python script to set up the background model.
The ``gammalib.GSPIModelDataSpace()`` constructor takes four arguments:

* a SPI observation group
* the name you want to give to the model component
* a fit method that should be used for the model component
* an index that refers to the template that should be used for the model
  component

In the example below it is assumed that three background templates
were generated using the
`Off-line Scientific Analysis (OSA) <https://www.isdc.unige.ch/integral/analysis#Software>`_
task ``spi_obs_back``, which appends the templates to the SPI observation
group ``og_spi.fits``.
The first template corresponds to the saturated event rates in the Germanium
detectors, the second is a constant rate template, and the third corresponds
to a linearly increasing rate with time.

  .. code-block:: python

     # Load SPI observation
     og = gammalib.GSPIObservation('og_spi.fits')

     # Set background model components
     bgm_model_1 = gammalib.GSPIModelDataSpace(og, 'GEDSAT', 'orbit,dete', 0)
     bgm_model_2 = gammalib.GSPIModelDataSpace(og, 'CONST',  'dete', 1)
     bgm_model_3 = gammalib.GSPIModelDataSpace(og, 'TIME1',  'dete', 2)

     # Append models
     models = gammalib.GModels()
     models.append(bgm_model_1)
     models.append(bgm_model_2)
     models.append(bgm_model_3)

     # Save models
     models.save('model_spi.xml')

For illustration, the typical content of the resulting XML file is as follows:

  .. code-block:: xml

    <source name="GEDSAT" type="DataSpace" instrument="SPI" method="orbit,dete" index="0">
      <parameter name="GEDSAT D000 O0019" value="1" scale="1" free="1"/>
      <parameter name="GEDSAT D001 O0019" value="1" scale="1" free="1"/>
      <parameter name="GEDSAT D002 O0019" value="1" scale="1" free="1"/>
      <parameter name="GEDSAT D003 O0019" value="1" scale="1" free="1"/>
      ...
      <parameter name="GEDSAT D016 O0139" value="1" scale="1" free="1"/>
      <parameter name="GEDSAT D017 O0139" value="1" scale="1" free="1"/>
      <parameter name="GEDSAT D018 O0139" value="1" scale="1" free="1"/>
    </source>
    <source name="CONST" type="DataSpace" instrument="SPI" method="dete" index="1">
      <parameter name="CONST D000" value="1" scale="1" free="1"/>
      <parameter name="CONST D001" value="1" scale="1" free="1"/>
      <parameter name="CONST D002" value="1" scale="1" free="1"/>
      <parameter name="CONST D003" value="1" scale="1" free="1"/>
      <parameter name="CONST D004" value="1" scale="1" free="1"/>
      <parameter name="CONST D005" value="1" scale="1" free="1"/>
      <parameter name="CONST D006" value="1" scale="1" free="1"/>
      <parameter name="CONST D007" value="1" scale="1" free="1"/>
      <parameter name="CONST D008" value="1" scale="1" free="1"/>
      <parameter name="CONST D009" value="1" scale="1" free="1"/>
      <parameter name="CONST D010" value="1" scale="1" free="1"/>
      <parameter name="CONST D011" value="1" scale="1" free="1"/>
      <parameter name="CONST D012" value="1" scale="1" free="1"/>
      <parameter name="CONST D013" value="1" scale="1" free="1"/>
      <parameter name="CONST D014" value="1" scale="1" free="1"/>
      <parameter name="CONST D015" value="1" scale="1" free="1"/>
      <parameter name="CONST D016" value="1" scale="1" free="1"/>
      <parameter name="CONST D017" value="1" scale="1" free="1"/>
      <parameter name="CONST D018" value="1" scale="1" free="1"/>
    </source>
    <source name="TIME1" type="DataSpace" instrument="SPI" method="dete" index="2">
      <parameter name="TIME1 D000" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D001" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D002" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D003" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D004" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D005" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D006" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D007" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D008" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D009" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D010" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D011" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D012" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D013" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D014" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D015" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D016" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D017" value="1" scale="1" free="1"/>
      <parameter name="TIME1 D018" value="1" scale="1" free="1"/>
    </source>


Fit methods
^^^^^^^^^^^

The fit method defines how the indices of the SPI data space are
grouped into template normalisation factors. Grouping is allowed in all three
data space dimensions: for pointings, detectors and energies.

In the example above, individual background normalisation factors for the
``GEDSAT`` template were assigned to each orbit and detector, while for the
``CONST`` and ``TIME1`` template there is one normalisation factor per
detector.  Here some further examples:

* ``orbit,dete,gedfail``: one normalisation factor per orbit and detector and
  Germanium detector failure period
* ``orbit,evtclass``: one normalisation factor per orbit and event class
* ``dete,ebin``: one normalisation factor per detector and energy bin
* ``gedanneal``: one normalisation factor per Germanium detector annealing
  period


Pointing grouping
~~~~~~~~~~~~~~~~~

point
.....

Assigns a template normalisation factor for each pointing.


orbit
.....

Assigns a template normalisation factor for each orbit.

date
....

Assigns a template normalisation factor for each date interval. For example,
``date 4 hours`` assigns a normalisation factor for intervals of 4 hours,
or ``date 1 month`` assigns a normalisation factor for intervals of one month.
The following time units are supported: ``min``, ``hour``, ``day``, ``week``,
``month`` or ``year``.

gedfail
.......

Makes sure that there are different template normalisation factors before and
after the failure of a Germanium detector. This method can either be used
as a fit method, or can be used in combination with ``point``, ``orbit``,
``date`` and ``gedanneal``.

gedanneal
.........

Makes sure that there are different template normalisation factors before and
after the Germanium detector annealing. This method can either be used
as a fit method, or can be used in combination with ``point``, ``orbit``,
``date`` and ``gedfail``.


Detector grouping
~~~~~~~~~~~~~~~~~

dete
....

Assigns a template normalisation factor for each detector.

evtclass
........

Assigns a template normalisation factor for event class. Event classes are
either single detector events (SE), double detector events (ME2) or tripple
detector events (ME3).


Energy grouping
~~~~~~~~~~~~~~~

ebin
....

Assigns a template normalisation factor for each energy bin.
