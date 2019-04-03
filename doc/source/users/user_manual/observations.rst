.. _um_observations:

Combining observations
----------------------

Data taking by IACTs is typically split into so-called *runs*, which we refer
here to *observations*. You typically want to combine these observations so
that all data available for a source are used for the analysis. There are two
basic methods how the combination can be done: either by listing all
observations to perform a joint analysis, or by combining all observations
for a co-called *stacked* analysis.


Joint data analysis
^^^^^^^^^^^^^^^^^^^

To perform a joint analysis of the observations you have to list them in a
so-called observation definition XML file. ctools will then compute the joint
log-likelihood

.. math::
   \ln L(M) = \sum_k \ln L_k(M)

for the maximum likelihood estimations, where the sum is taken over the
individual observations :math:`k`. Below you will see the first lines of an
:ref:`observation definition XML file <glossary_obsdef>` that describe the
first three observations of the Galactic Plane survey for the first CTA Data
Challenge:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <observation_list title="observation list">
     <observation name="GPS" id="110000" instrument="CTA">
       <parameter name="EventList" file="$CTADATA/data/baseline/gps/gps_baseline_110000.fits" />
       <parameter name="Calibration" database="1dc" response="South_z40_50h" />
     </observation>
     <observation name="GPS" id="110001" instrument="CTA">
       <parameter name="EventList" file="$CTADATA/data/baseline/gps/gps_baseline_110001.fits" />
       <parameter name="Calibration" database="1dc" response="South_z40_50h" />
     </observation>
     <observation name="GPS" id="110002" instrument="CTA">
       <parameter name="EventList" file="$CTADATA/data/baseline/gps/gps_baseline_110002.fits" />
       <parameter name="Calibration" database="1dc" response="South_z40_50h" />
     </observation>
     ...
   </observation_list>

Each observation has an ``instrument`` and an ``id`` attribute that allow to
uniquely identify the observation. For a given instrument label, the ``id``
attribute must be unique.

Also observations that were binned into counts cubes using :ref:`ctbin` can
be analysed jointly by specifying the counts cubes for all observations in the
observation definiton XML file:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <observation_list title="observation list">
     <observation name="GPS" id="110000" instrument="CTA">
       <parameter name="CountsCube" file="cntcube_cta_110000.fits" />
       <parameter name="Calibration" database="1dc" response="South_z40_50h" />
     </observation>
     <observation name="GPS" id="110001" instrument="CTA">
       <parameter name="CountsCube" file="cntcube_cta_110001.fits" />
       <parameter name="Calibration" database="1dc" response="South_z40_50h" />
     </observation>
     <observation name="GPS" id="110002" instrument="CTA">
       <parameter name="CountsCube" file="cntcube_cta_110002.fits" />
       <parameter name="Calibration" database="1dc" response="South_z40_50h" />
     </observation>
     ...
   </observation_list>

Equally, observations prepared for On/Off analysis using :ref:`csphagen` can
also be analysed jointly:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <observation_list title="observation list">
     <observation name="GPS" id="110000" instrument="CTAOnOff" statistic="wstat">
       <parameter name="Pha_on"  file="onoff_110000_pha_on.fits" />
       <parameter name="Pha_off" file="onoff_110000_pha_off.fits" />
       <parameter name="Arf"     file="onoff_110000_arf.fits" />
       <parameter name="Rmf"     file="onoff_110000_rmf.fits" />
     </observation>
     <observation name="GPS" id="110001" instrument="CTAOnOff" statistic="wstat">
       <parameter name="Pha_on"  file="onoff_110001_pha_on.fits" />
       <parameter name="Pha_off" file="onoff_110001_pha_off.fits" />
       <parameter name="Arf"     file="onoff_110001_arf.fits" />
       <parameter name="Rmf"     file="onoff_110001_rmf.fits" />
     </observation>
     <observation name="GPS" id="110002" instrument="CTAOnOff" statistic="wstat">
       <parameter name="Pha_on"  file="onoff_110002_pha_on.fits" />
       <parameter name="Pha_off" file="onoff_110002_pha_off.fits" />
       <parameter name="Arf"     file="onoff_110002_arf.fits" />
       <parameter name="Rmf"     file="onoff_110002_rmf.fits" />
     </observation>
     ...
   </observation_list>

.. note::
   For On/Off analysis the instrument label has to be suffixed by ``OnOff``,
   i.e. ``CTA`` becomes ``CTAOnOff``, ``HESS`` becomes ``HESSOnOff`` and
   so on.

.. note::
   The optional ``statistic`` attribute allows to specify for each
   observation which statistic should be used. Possible values for unbinned
   or binned analysis are ``poisson``, ``cstat`` (equivalent to ``poisson``),
   ``gaussian`` or ``chi2`` (equivalent to ``gaussian``). Possible values for
   On/Off analysis are ``poisson``, ``cstat`` (equivalent to ``poisson``) or
   ``wstat``.

Observations of different instruments can also be combined for a joint
analysis, as is illustrated in the example below:

.. code-block:: xml

   <?xml version="1.0" standalone="no"?>
   <observation_list title="observation library">
     <observation name="Crab" id="000001" instrument="COM">
       <parameter name="DRE" file="m50439_dre.fits" />
       <parameter name="DRB" file="m34997_drg.fits" />
       <parameter name="DRG" file="m34997_drg.fits" />
       <parameter name="DRX" file="m32171_drx.fits" />
       <parameter name="IAQ" file="ENERG(1.0-3.0)MeV" />
     </observation>
     <observation name="Crab" id="000001" instrument="LAT">
       <parameter name="CountsMap"    file="srcmap.fits" />
       <parameter name="ExposureMap"  file="expmap.fits" />
       <parameter name="LiveTimeCube" file="ltcube.fits" />
       <parameter name="IRF"          value="P8R3_SOURCE_V2" />
     </observation>
     <observation name="Crab" id="000001" instrument="CTA">
       <parameter name="EventList"   file="cta_events.fits" />
       <parameter name="Calibration" database="prod2" response="South_0.5h" />
     </observation>
   </observation_list>

You may also combined unbinned, binned, stacked (see below) or On/Off data in
a joint analysis.


Stacked data analysis
^^^^^^^^^^^^^^^^^^^^^

To speed-up the computations, binned observations can be stacked together,
requiring the computation of average response functions that properly weight
the individual instrument response functions for each observation. Stacking
of observations is possible for binned counts cubes and On/Off spectra.


Stacked counts cubes
~~~~~~~~~~~~~~~~~~~~

By default, :ref:`ctbin` stacks the events from multiple observations into
a single counts cube (to produce counts cubes for each individual observation,
:ref:`ctbin` must be executed with ``stack=no`` option). The response for a
stacked analysis is composed of an exposure cube, a point spread function
cube, an energy dispersion cube and a background cube.

The exposure cube is computed by :ref:`ctexpcube` using

.. math::
   X_\mathrm{cube}(p,E) = \sum_k A_{\mathrm{eff},k}(p,E,t) \times \tau_k

where :math:`A_{\mathrm{eff},k}(p,E,t)` is the effective area and
:math:`\tau_k` the live time of observation :math:`k`. The point spread
function cube is computed by :ref:`ctpsfcube` using

.. math::
   \mathrm{\it PSF}_\mathrm{cube}(p,E,\delta) =
   \frac{\sum_k \mathrm{\it PSF}_k(p'|p,E,t) \times
         A_{\mathrm{eff},k}(p,E,t) \times \tau_k}
        {\sum_k A_{\mathrm{eff},k}(p,E,t) \times \tau_k}

where :math:`\mathrm{\it PSF}_k(p'|p,E,t)` is the point spread function of
observation :math:`k`. The energy dispersion cube is computed by
:ref:`ctedispcube` using

.. math::
   \mathrm{\it D}_\mathrm{cube}(E'|p,E) =
   \frac{\sum_k E_{\mathrm{disp},k}(E'|p,E,t) \times
         A_{\mathrm{eff},k}(p,E,t) \times \tau_k}
        {\sum_k A_{\mathrm{eff},k}(p,E,t) \times \tau_k},

where :math:`E_{\mathrm{disp},k}(E'|p,E,t)` is the energy dispersion of
observation :math:`k`. The background cube is computed by :ref:`ctbkgcube`
using

.. math::
   B_\mathrm{cube}(p', E') = \frac{\sum_k B_k(p',E',t') \times \tau_k}
                                  {\sum_k \tau_k}

where :math:`B_k(p',E',t')` is the background model of observation :math:`k`.
The sum is taken over all observations :math:`k`.

The files for a stacked binned observation are specified as follows in an
:ref:`observation definition XML file <glossary_obsdef>`:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <observation_list title="observation list">
     <observation name="Crab" id="000001" instrument="CTA">
       <parameter name="CountsCube"   file="cntcube.fits" />
       <parameter name="ExposureCube" file="expcube.fits" />
       <parameter name="PsfCube"      file="psfcube.fits" />
       <parameter name="EdispCube"    file="edispcube.fits" />
       <parameter name="BkgCube"      file="bkgcube.fits" />
     </observation>
   </observation_list>

(the ``EdispCube`` parameter is optional).


Stacked On/Off spectra
~~~~~~~~~~~~~~~~~~~~~~

By default, :ref:`csphagen` creates On/Off spectra and response information
for each individual observation, but executing :ref:`csphagen` with the
``stack=yes`` option will lead to the generation of stacked On/Off spectra.
This leads to the combination of events from all observations :math:`k` into
a single On and Off ``PHA`` spectrum

.. math::
   n^\mathrm{on}_i = \sum_k n^\mathrm{on}_{k,i}

   n^\mathrm{off}_i = \sum_k n^\mathrm{off}_{k,i}

where :math:`i` is the energy bin and :math:`k` is the observation.
The effective ``ARF`` is computed using

.. math::
   ARF(E) = \frac{\sum_k ARF_k(E) \times \tau_k}{\sum_k \tau_k}

the effective ``RMF`` using

.. math::
   RMF_{i}(E) = \frac{\sum_k RMF_{k,i}(E) \times ARF_k(E) \times \tau_k}
                     {\sum_k ARF_k(E) \times \tau_k}

the effective background scaling factors using

.. math::
   \alpha_i = \frac{\sum_k \alpha_{k,i} \times b_{k,i} \times \tau_k}
                   {\sum_k b_{k,i} \times \tau_k}

and the effective background response vectors using

.. math::
   b_i = \frac{\sum_k b_{k,i} \times \tau_k}{\sum_k \tau_k}

where :math:`\tau_k` is the live time (or exposure) of observation :math:`k`.

The :ref:`observation definition XML file <glossary_obsdef>` for a stacked
On/Off observation has the same format as for a single On/Off observation:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <observation_list title="observation list">
     <observation name="GPS" id="000001" instrument="CTAOnOff" statistic="wstat">
       <parameter name="Pha_on"  file="onoff_stacked_pha_on.fits" />
       <parameter name="Pha_off" file="onoff_stacked_pha_off.fits" />
       <parameter name="Arf"     file="onoff_stacked_arf.fits" />
       <parameter name="Rmf"     file="onoff_stacked_rmf.fits" />
     </observation>
   </observation_list>

