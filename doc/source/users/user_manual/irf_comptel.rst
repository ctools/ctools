.. _um_irf_comptel:

COMPTEL response functions
--------------------------

Formulation
~~~~~~~~~~~

A COMPTEL event is characterised by an instrument direction spanned by the
angles :math:`(\chi, \psi, \bar{\varphi})`. :math:`(\chi, \psi)` is the
direction of the photon after scattering in the upper detector plane, which is
determined from the photon interaction locations in both detector planes, and

.. math::
   \bar{\varphi} = \arccos \left( 1 - \frac{m_\mathrm{e}c^2}{E_2} + \frac{m_\mathrm{e}c^2}{E_1+E_2} \right)

is the Compton scattering angle as inferred from the energy deposits :math:`E_1`
and :math:`E_2` in the upper and lower detector planes, respectively.
The measured energy of the photon is estimated from the sum

.. math::
   E' = E_1 + E_2

of the energy deposits in both detector planes. The probability that a photon
which interacted in the upper detector plane will encounter a detector of the
lower plane is described by :math:`DRG(\chi, \psi, \bar{\varphi})`, which also
includes any cuts related to the removal of events coming from the Earth limb.

The COMPTEL response is factorised using

.. math::
   R(p',E',t'|p,E,t) = \frac{\tau}{T} \frac{DRX(p)}{T} \times DRG(\chi, \psi, \bar{\varphi})
                       \times IAQ(\chi, \psi, \bar{\varphi} | p, E)

where
:math:`\tau` is the lifetime in units of :math:`s`,
:math:`T` is the ontime in units of :math:`s`,
:math:`DRX(p)` is the exposure in units of :math:`cm^2 \, s`, and
:math:`IAQ(\chi, \psi, \bar{\varphi} | p, E)` quantifies the interaction
probability for a Compton scattering in the upper detector plane followed by
an interaction in the lower detector plane.

We note that :math:`IAQ(\chi, \psi, \bar{\varphi} | p, E)` is azimuthally
symmetric about the source direction, and the IAQ file is stored as a 2D FITS
image providing the interaction probabilities as function of
:math:`\bar{\varphi}` and :math:`\varphi_\mathrm{geo}` for a given energy
range, where :math:`\varphi_\mathrm{geo}` is the angular separation between
:math:`(\chi, \psi)` and :math:`p`.


Analytical Instrument Response Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

COMPTEL instrument response functions may be computed analytically, taken into
account ground calibration information. Analytical computation is implemented
in GammaLib, and is automatically performed when :ref:`comobsbin` is executed.


Simulated Instrument Response Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively, simulated instrument response functions can be used for the data
analysis. Simulated instrument response functions can be specified in the
:ref:`observation definition file <glossary_obsdef>` through the ``IAQ``
parameter as illustrated below.

.. code-block:: xml

   <?xml version="1.0" standalone="no"?>
   <observation_list title="observation library">
     <observation name="Crab" id="100001" instrument="COM">
       <parameter name="DRE" file="m50438_dre.fits"/>
       <parameter name="DRB" file="bgdlix_drb.fits"/>
       <parameter name="DRG" file="m34997_drg.fits"/>
       <parameter name="DRX" file="m32171_drx.fits"/>
       <parameter name="IAQ" value="SIM2(0.75-1.00)MeV(2)deg"/>
     </observation>
     ...
   </observation_list>

A number of simulated Instrument Response Functions are shipped together with
GammaLib. The table below specifies the possible options for the ``value``
field of COMPTEL :ref:`observation definition files <glossary_obsdef>`.

============================== ============= ============== =========
``value``                      Energy range  Phibar binning Type
============================== ============= ============== =========
``SIM2(0.75-1.00)MeV(2)deg``   0.75-1 MeV    2 deg          Continuum
``SIM2(1.00-3.00)MeV(2)deg``   1-3 MeV       2 deg          Continuum
``SIM2(3.00-10.00)MeV(2)deg``  3-10 MeV      2 deg          Continuum
``SIM2(10.00-30.00)MeV(2)deg`` 10-30 MeV     2 deg          Continuum
``SIM2(0.75-30.00)MeV(2)deg``  0.75-30 MeV   2 deg          Continuum
``SIM2(1.00-30.00)MeV(2)deg``  1-30 MeV      2 deg          Continuum
``SIM2(1.809)MeV(2)deg``       1.809 MeV     2 deg          Line
``SIM3(0.75-1.00)MeV(2)deg``   0.75-1 MeV    2 deg          Continuum
``SIM3(1.00-3.00)MeV(2)deg``   1-3 MeV       2 deg          Continuum
``SIM3(3.00-10.00)MeV(2)deg``  3-10 MeV      2 deg          Continuum
``SIM3(10.00-30.00)MeV(2)deg`` 10-30 MeV     2 deg          Continuum
``SIM3(0.75-1.00)MeV(1)deg``   0.75-1 MeV    1 deg          Continuum
``SIM3(1.00-3.00)MeV(1)deg``   1-3 MeV       1 deg          Continuum
``SIM3(3.00-10.00)MeV(1)deg``  3-10 MeV      1 deg          Continuum
``SIM3(10.00-30.00)MeV(1)deg`` 10-30 MeV     1 deg          Continuum
``SIM3(1.00-1.25)MeV(1)deg``   1-1.25 MeV    1 deg          Continuum
``SIM3(1.25-1.50)MeV(1)deg``   1.25-1.5 MeV  1 deg          Continuum
``SIM3(1.50-2.00)MeV(1)deg``   1.5-2 MeV     1 deg          Continuum
``SIM3(2.00-2.50)MeV(1)deg``   2-2.5 MeV     1 deg          Continuum
``SIM3(2.50-3.00)MeV(1)deg``   2.5-3 MeV     1 deg          Continuum
``SIM3(3.00-4.00)MeV(1)deg``   3-4 MeV       1 deg          Continuum
``SIM3(4.00-6.00)MeV(1)deg``   4-6 MeV       1 deg          Continuum
``SIM3(6.00-8.00)MeV(1)deg``   6-8 MeV       1 deg          Continuum
``SIM3(8.00-10.00)MeV(1)deg``  8-10 MeV      1 deg          Continuum
``SIM3(10.00-15.00)MeV(1)deg`` 10-15 MeV     1 deg          Continuum
``SIM3(15.00-30.00)MeV(1)deg`` 15-30 MeV     1 deg          Continuum
``SIM3(0.75-0.90)MeV(1)deg``   0.75-0.9 MeV  1 deg          Continuum
``SIM3(0.90-1.06)MeV(1)deg``   0.9-1.06 MeV  1 deg          Continuum
``SIM3(1.06-1.28)MeV(1)deg``   1.06-1.28 MeV 1 deg          Continuum
``SIM3(1.28-1.50)MeV(1)deg``   1.28-1.50 MeV 1 deg          Continuum
``SIM3(1.50-1.70)MeV(1)deg``   1.50-1.70 MeV 1 deg          Continuum
``SIM3(1.70-1.90)MeV(1)deg``   1.70-1.90 MeV 1 deg          Continuum
``SIM3(1.90-2.10)MeV(1)deg``   1.90-2.10 MeV 1 deg          Continuum
``SIM3(2.10-2.30)MeV(1)deg``   2.10-2.30 MeV 1 deg          Continuum
============================== ============= ============== =========
