.. _sec_glossary:

Glossary
========

.. _glossary_bkgcube:

.. topic:: Background cube

   A background cube is a 3-dimensional cube spanned by reconstructed Right
   Ascension or Galactic longitude, reconstructed Declination or Galactic
   latitude, and reconstructed photon energy.
   The background cube :math:`B(p', E')` is computed
   using

   .. math::
      B(p', E') = \frac{\sum_i B_i(p', E') \times \tau_i}
                       {\sum_i \tau_i}

   where
   :math:`B_i(p', E')` is the background rate estimate for observation :math:`i`,
   :math:`p'` is the reconstruction photon arrival direction,
   :math:`E'` is the reconstruction photon energy,
   and :math:`\tau_i` is the livetime of observation :math:`i`.
   The sum is taken over all observations.


.. _glossary_countscube:

.. topic:: Counts cube

   A counts cube is a 3-dimensional data cube spanned by reconstructed Right
   Ascension or Galactic longitude, reconstructed Declination or Galactic
   latitude, and reconstructed energy.


.. _glossary_edispcube:

.. topic:: Energy Dispersion cube

   An energy dispersion cube is a 4-dimensional cube spanned by true Right
   Ascension or Galactic longitude, true Declination or Galactic latitude, true
   photon energy, and migration which is the ratio between reconstructed and true
   photon energy.
   The energy dispersion cube :math:`E_{\rm disp}(E'| p, E)` is computed
   using

   .. math::
      E_{\rm disp}(E'| p, E) =
      \frac{\sum_i E_{\rm disp,i}(E'| p, E) \times A_{\rm eff,i}(p, E) \times \tau_i}
           {\sum_i A_{\rm eff,i}(p, E) \times \tau_i}

   where
   :math:`E_{\rm disp,i}(E'| p, E)` is the energy dispersion for observation :math:`i`,
   :math:`E'` is the reconstruction photon energy,
   :math:`p` is the true photon arrival direction,
   :math:`E` is the true photon energy,
   and :math:`\tau_i` is the livetime of observation :math:`i`.
   The sum is taken over all observations.


.. _glossary_eventlist:

.. topic:: Event list

   An event list is a table comprising the reconstructed properties of the
   observed events. Each row of the table corresponds to an event. The columns
   of the table provide event characteristics, such as trigger time,
   reconstructed arrival direction (in Right Ascension and Declination), and
   the reconstructed event energy.


.. _glossary_expcube:

.. topic:: Exposure cube

   An exposure cube is a 3-dimensional cube spanned by true Right Ascension or
   Galactic longitude, true Declination or Galactic latitude, and true photon
   energy. The exposure
   :math:`X(p, E)` is computed using

   .. math::
      X(p, E) = \sum_i A_{\rm eff,i}(p, E) \times \tau_i

   where
   :math:`A_{\rm eff,i}(p, E)` is the effective area for observation :math:`i`,
   :math:`p` is the true photon arrival direction,
   :math:`E` is the true photon energy,
   and :math:`\tau_i` is the livetime of observation :math:`i`.
   The sum is taken over all observations.


.. _glossary_1dc:

.. topic:: First CTA Data Challenge

   The goal of the first CTA Data Challenge is to enable the CTA Consortium
   Science Working Groups to derive science benchmarks for the CTA Key Science
   objectives.

   The first CTA Data Challenge should provide quantitative estimates of CTA's
   science capabilities that will enable the CTA Consortium Science Working
   Groups to evaluate science trade-offs in the future. The first CTA Data
   Challenge should also result in numerous show cases, including for example
   images, spectra and light curves, that can be used to illustrate CTA's
   science case, and that should enrich the CTA outreach material. And the
   first CTA Data Challenge should also stimulate the enrichment of the CTA
   Science Case.


.. _glossary_gti:

.. topic:: Good Time Interval

   A Good Time Interval is a contiguous time period, defined by a start and
   a stop time, during which the events can be used for a scientific analysis.


.. _glossary_irf:

.. topic:: Instrument Response Functions

   The instrument response functions provide a mathematical description that
   links the measured quantities of an event to the physical quantities of
   the incident photon. The instrument response functions for CTA are factorised
   into the effective area :math:`A_{\rm eff}(p, E, t)` (units :math:`cm^2`),
   the point spread function :math:`PSF(p' | p, E, t)`,
   and the energy dispersion :math:`E_{\rm disp}(E' | p, E, t)`
   following:

   .. math::
      R(p', E', t' | p, E, t) =
      A_{\rm eff}(p, E, t) \times
      PSF(p' | p, E, t) \times
      E_{\rm disp}(E' | p, E, t)


.. _glossary_moddef:

.. topic:: Model definition file

   Source and background model components are defined in ctools by a model
   definition file. The model definition file is an ASCII file in XML format.
   The format of the model definition file is inspired from, and is
   compatible with, the format used by the Fermi-LAT Science Tools.
   The general structure of a model definition file is

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
        <source name="Background" type="CTAIrfBackground" instrument="CTA">
          <spectrum type="PowerLaw">
            <parameter name="Prefactor"   scale="1.0"  value="1.0"  min="1e-3" max="1e+3"   free="1"/>
            <parameter name="Index"       scale="1.0"  value="0.0"  min="-5.0" max="+5.0"   free="1"/>
            <parameter name="PivotEnergy" scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
          </spectrum>
        </source>
      </source_library>

   Each model component is described by the ``<source>`` tag. Each source has
   an mandatory spectral and a spatial component (tags ``<spectrum>`` and
   ``<spatialModel>``) and an optional temporal component (tag ``<temporal>``).
   Parameters that should be adjusted in a maximum likelihood fit should be
   set to ``free="1"``; otherwise they are hold fixed.


.. _glossary_obs:

.. topic:: Observation

   The data are split into observations. Each observation is characterised by
   a stable instrument configuration that can be described by an instrument
   response function. Observations are also known as runs.


.. _glossary_obsdef:

.. topic:: Observation definition file

   Observations are combined in ctools using an observation definition file.
   The observation definition file is an ASCII file in XML format.
   The format of the observation definition file is illustrated below:

   .. code-block:: xml

      <?xml version="1.0" standalone="no"?>
      <observation_list title="observation library">
        <observation name="Crab" id="00001" instrument="CTA">
          <parameter name="EventList" file="events1.fits"/>
        </observation>
        <observation name="Crab" id="00002" instrument="CTA">
          <parameter name="EventList" file="events2.fits"/>
        </observation>
      </observation_list>


.. _glossary_psfcube:

.. topic:: Point Spread Function cube

   A point spread function cube is a 4-dimensional cube spanned by true Right
   Ascension or Galactic longitude, true Declination or Galactic latitude, true
   photon energy, and offset angle between true and reconstructed arrival
   direction of a photon.
   The point spread function cube :math:`PSF(\delta | p, E)` is computed
   using

   .. math::
      PSF(\delta | p, E) = \frac{\sum_i PSF_i(p' | p, E) \times A_{\rm eff,i}(p, E) \times \tau_i}
                                {\sum_i A_{\rm eff,i}(p, E) \times \tau_i}

   where
   :math:`\delta` is the angular separation between the true and measured
   photon directions :math:`p` and :math:`p'`, respectively,
   :math:`A_{\rm eff,i}(p, E)` is the effective area for observation :math:`i`,
   :math:`E` is the true photon energy,
   and :math:`\tau_i` is the livetime of observation :math:`i`.
   The sum is taken over all observations.
