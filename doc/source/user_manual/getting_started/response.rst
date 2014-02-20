CTA Instrument Response Functions
---------------------------------

.. note ::

   Instrument response functions are formally not part of ctools as they
   should be provided by the instrument teams. The GammaLib framework on
   which ctools are built comes with a set of basic response functions, but
   for getting the latest instrument response functions you should contact
   the relevant instrument teams. CTA consortium members should 
   `download 
   <https://portal.cta-observatory.org/WG/DM/DM_wiki/DATA_Access/Pages/Science%20Tools.aspx>`_
   the latest calibration database containing Prod1 and Prod2 instrument
   response functions from the consortium Sharepoint site (requires 
   password).


What are instrument response functions?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The instrument response functions provide a mathematical description that
links the measured quantities :math:`\vec{d}` of an event to the physical
quantities :math:`\vec{p}` of the incident photon. The following figure 
illustrates this relationship:

.. figure:: irfs.jpg
   :width: 70%

:math:`I(\vec{p})` is the gamma-ray intensity arriving at Earth as a
function of photon properties :math:`\vec{p}` 
(which usually are true photon energy, true photon incident direction, 
and true photon arrival time),
while :math:`e(\vec{d})` is the expected event rate as function of event 
properties :math:`\vec{d}` (which usually are the measured photon energy,
measured or reconstructued photon incident direction, and measured photon 
arrival time). The expected event rate is obtained by integrating the
product of 
the instrumental response function :math:`R(\vec{d}|\vec{p},\vec{a}`)
and the emitted intensity :math:`I(\vec{p})` over the photon properties
:math:`\vec{p}`.
The argument :math:`\vec{a}` in the response function comprises any 
auxiliary parameter on which the response function may depend on (e.g. 
pointing direction, triggered telescopes, optical efficiencies, 
atmospheric conditions, etc.). All these quantities and hence the 
instrument response function may depend on time.


Data formats
~~~~~~~~~~~~

So far the data format of the CTA instrument response functions has not
yet been defined, and information about instrument response functions is
provided in various ways. Three types of data formats are so far supported:

-  :ref:`sec_cta_rsp_perftable`

-  :ref:`sec_cta_rsp_xspec`

-  :ref:`sec_cta_rsp_rsptable`

All formats provide the instrument response as function of true (and 
sometimes also measured) photon energy, typically from about 20 GeV to
about 125 TeV.
Performance tables provide the instrument response for on-axis sources.
ARF, RMF and PSF files provide the instrument response for a specific
source position (and region) within the field of view.
Response tables provide the instrument response as function of position in 
the field of view. Response tables are thus the most universal form of 
instrument response functions available, and we recommend using this form 
for any analysis.

.. note ::

   Instrument response functions are so far only available for a fixed
   zenith angle of 20 deg. Therefore, no zenith or azimuth angle dependence has 
   been implented in the CTA response functions, but the dependence can be 
   handled by cutting the data into short time segments (typically of 30 
   minutes in length), and by specifying specific response functions for 
   each segment.


Installing the CTA calibration database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After `downloading 
<https://portal.cta-observatory.org/WG/DM/DM_wiki/DATA_Access/Pages/Science%20Tools.aspx>`_
the latest calibration database (only possible for CTA consortium members),
the database is installed using

.. code-block:: bash

  $ [sudo] tar -C $CALDB -zxvf cta-caldb-20140216.tar.gz

or

.. code-block:: bash

  $ [sudo] tar -C $CALDB -xvf cta-caldb-20140216.tar

depending on whether the database has been retrieved as a gzipped file or 
not. Depending on your access rights to the ``$CALDB`` directory,
installation of the calibration database may require root privileges
(type ``sudo`` in this case, otherwise omit this part of the command).

The calibration database must be installed into the directory to which you 
``CALDB`` environment variable points. By default, ctools sets this 
enviroment variable to

.. code-block:: bash

  $CTOOLS/share/caldb

but it may be that some other software installation or configuration 
setting on your system overwrites this location. If in doubt, type

.. code-block:: bash

  echo $CALDB

to find out to which directory your ``CALDB`` environment variable is set. 
You may always overwrite this setting using

.. code-block:: bash

  export CALDB=/my/preferred/caldb/directory

which you can add to your ``.bashrc`` file to set the directory 
permanently.
ctools will use the ``CALDB`` environment variable to find out where the 
calibration database is located. If needed, this mechanism can be 
circumvented (see :ref:`sec_cta_rsp_abspath`).


Specifying the CTA Instrument Response Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The specification of the CTA Instrument Response Functions depends on the 
way how ctools are used. Common to all methods is that the IRFs are 
defined by a response name and a calibration database name. The latter 
may in some cases be the path to a directory on your filesystem.

There are different means to specify the CTA Instrument Response Functions 
when using ctools, and the following section describe the

Using individual event files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(to be written)


Using observation definition files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(to be written)


From within a Python script
^^^^^^^^^^^^^^^^^^^^^^^^^^^

(to be written)


.. _sec_cta_rsp_abspath:

Using absolute path names to instrument response files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(to be written)


Data format details
~~~~~~~~~~~~~~~~~~~

.. _sec_cta_rsp_perftable:

Performance tables
^^^^^^^^^^^^^^^^^^

In the early days, the instrument performances derived from Monte-Carlo
simulations have been summarised in what we call here Performance 
Tables, which are ASCII files that contain as function of energy the
on-axis performance parameters of CTA, such as effective area, point spread
function containment radius, energy resolution, background count rate and
differential sensitivity.

Below an example of a CTA performance table::

  log(E)     Area     r68     r80  ERes. BG Rate    Diff Sens
  -1.7      261.6  0.3621  0.4908 0.5134 0.0189924  6.88237e-11
  -1.5     5458.2  0.2712  0.3685 0.4129 0.1009715  1.72717e-11
  -1.3    15590.0  0.1662  0.2103 0.2721 0.0575623  6.16963e-12
  -1.1    26554.1  0.1253  0.1567 0.2611 0.0213008  2.89932e-12
  -0.9    52100.5  0.1048  0.1305 0.1987 0.0088729  1.39764e-12
  -0.7    66132.1  0.0827  0.1024 0.1698 0.0010976  6.03531e-13
  -0.5   108656.8  0.0703  0.0867 0.1506 0.0004843  3.98147e-13
  -0.3   129833.0  0.0585  0.0722 0.1338 0.0001575  3.23090e-13
  -0.1   284604.3  0.0531  0.0656 0.1008 0.0001367  2.20178e-13
   0.1   263175.3  0.0410  0.0506 0.0831 0.0000210  1.87452e-13
   0.3   778048.6  0.0470  0.0591 0.0842 0.0000692  1.53976e-13
   0.5   929818.8  0.0391  0.0492 0.0650 0.0000146  1.18947e-13
   0.7  1078450.0  0.0335  0.0415 0.0541 0.0000116  1.51927e-13
   0.9  1448579.1  0.0317  0.0397 0.0516 0.0000047  1.42439e-13
   1.1  1899905.0  0.0290  0.0372 0.0501 0.0000081  1.96670e-13
   1.3  2476403.8  0.0285  0.0367 0.0538 0.0000059  2.20695e-13
   1.5  2832570.6  0.0284  0.0372 0.0636 0.0000073  3.22523e-13
   1.7  3534065.3  0.0290  0.0386 0.0731 0.0000135  4.84153e-13
   1.9  3250103.4  0.0238  0.0308 0.0729 0.0000044  6.26265e-13
   2.1  3916071.6  0.0260  0.0354 0.0908 0.0000023  7.69921e-13
   ---------------------------------------------
   1) log(E) = log10(E/TeV) - bin centre
   2) Eff Area - in square metres after background cut (no theta cut)
   3) Ang. Res - 68% containment radius of gamma-ray PSF post cuts - in degrees
   4) Ang. Res - 80% containment radius of gamma-ray PSF post cuts - in degrees
   5) Fractional Energy Resolution (rms)
   6) BG Rate  - inside point-source selection region - post call cuts - in Hz
   7) Diff Sens - differential sensitivity for this bin expressed as E^2 dN/dE
      - in erg cm^-2 s^-1 - for a 50 hours exposure - 5 sigma significance including
      systematics and statistics and at least 10 photons.


.. _sec_cta_rsp_xspec:

ARF, RMF and PSF files
^^^^^^^^^^^^^^^^^^^^^^

Instrument response information for the first CTA Data Challenge (1DC) has been
provided in a format that was heavily inspired by the 
`RMF and ARF file formats 
<http://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html>`_
that have been introduced by
`HEASARC <http://heasarc.gsfc.nasa.gov/>`_
for spectral analysis of X-ray data.

The Redistribution Matrix File (RMF) describes how an incoming photon with 
a given true energy is redistributed in measured energy. In other words, it
describes the energy dispersion of the instrument. The RMF is organised as
a two-dimensional matrix, with the first axis being in true energies while
the second axis represents the measured energies.

The Ancillary Response File (ARF) describes the sensitivity of the 
instrument to photons of a given true energy. The ARF gives the effective
area of the detector system after applying an event selection cut
(theta cut). For computing the ARF, some knowledge of the Point Spread
Function (PSF) is needed, so that the fraction of photons that falls into
the event selection region for a given source can be properly estimated.
Once this is done, no further PSF information is needed for the analysis.

ctools however relies on the full modelling of the instrument, including 
the Point Spread Function, and hence, PSF information has also been 
provided for 1DC. Two data formats have been used for this: a first that 
is based on a simple one-dimensional vector, providing the width of a
2-dimensional Gaussian as function of energy; and a second that is based
on a three-component 2D Gaussian as function of energy and offset angle
in the camera system. While the former was implement by a simply FITS 
table with two columns, the latter was implemented by
:ref:`sec_cta_rsp_rsptable`.


.. _sec_cta_rsp_rsptable:

Response tables
^^^^^^^^^^^^^^^

The CTA response table class ``GCTAResponseTable`` provides a generic 
handle for multi-dimensional response information. It is based on the 
response format used for storing response information for the
*Fermi*/LAT telescope. In this format, all information is stored in
a single row of a FITS binary table. Each element of the row contains
a vector column, that describes the axes of the  multi-dimensional response
cube and the response information. Note that this class may in the future
be promoted to the GammaLib core, as a similar class has been implemented
in the *Fermi*/LAT interface. 
