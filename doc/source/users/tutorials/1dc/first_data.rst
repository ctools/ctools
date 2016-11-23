.. _1dc_first_data:

Data organisation
-----------------

After downloading the
:ref:`First CTA Data Challenge <glossary_1dc>`
data from **TBD** and uncompressing the tarball using

.. code-block:: bash

   $ tar xfvz 1dc.tar.gz

you will have a directory ``1dc`` on your disk with the following content:

.. code-block:: bash

   caldb/
   caldb/data
   caldb/data/cta
   caldb/data/cta/prod2
   ...
   caldb/data/cta/prod3
   caldb/data/cta/prod3/caldb.indx
   caldb/data/cta/prod3/bcf
   ...
   data/
   data/gps_baseline_000001.fits
   data/gps_baseline_000002.fits
   data/gps_baseline_000003.fits
   ...
   data/gps_threshold_000001.fits
   data/gps_threshold_000002.fits
   data/gps_threshold_000003.fits
   ...
   data/egal_baseline_000001.fits
   data/egal_baseline_000002.fits
   data/egal_baseline_000003.fits
   ...
   data/egal_threshold_000001.fits
   data/egal_threshold_000002.fits
   data/egal_threshold_000003.fits
   ...
   data/obs_gps_baseline.xml
   data/obs_gps_threshold.xml
   data/obs_egal_baseline.xml
   data/obs_egal_threshold.xml
   data/obs_all_baseline.xml
   data/obs_all_threshold.xml
   models/
   models/models_gps.xml
   models/models_egal.xml
   models/map_RXJ1713.fits.gz
   models/map_VelaJunior.fits.gz
   models/map_ics.fits.gz
   models/map_pi0.fits.gz

The ``caldb`` folder contains the
:ref:`Instrument Response Functions <glossary_irf>`
that are necessary for the analysis of the simulated CTA data.
The folder contains the ``prod3`` response that should be used for the
:ref:`First CTA Data Challenge <glossary_1dc>`.
For legacy reasons, also the ``prod2`` response functions are shipped with
the tarball.

The ``data`` folder contains the calibrated, reconstructed and background
reduced event data that were procuded for the
:ref:`First CTA Data Challenge <glossary_1dc>`.
Each event file corresponds to the
:ref:`event lists <glossary_eventlist>`
for an
:ref:`observation (or run) <glossary_obs>`, and contains for each event a
time stamp, the reconstructued photon arrival direction, and the reconstructed
photon energy.
The pointing direction during an observation is fixed.
Event files exist for the Galactic Plane Survey (starting with ``gps_``) and
for the extragalactic survey (starting with ``egal_``) for both
the ``baseline`` configuration and the ``threshold`` implementation.

The files ``obs_XXX_YYYYY.xml`` are so called
:ref:`Observation Definition Files <glossary_obsdef>`
that contain the information (or metadata) of a list of observations.
The files are plain ASCII files in XML format that be inspected and
manipulated by any text editor.
The following files are included in the data distribution:

* ``obs_gps_baseline.xml``: Galactic Plane Survey performed with baseline arrays
* ``obs_gps_threshold.xml``: Galactic Plane Survey performed with threshold implementation
* ``obs_egal_baseline.xml``: Extragalactic Survey performed with baseline arrays
* ``obs_egal_threshold.xml``: Extragalactic Survey performed with threshold implementation
* ``obs_all_baseline.xml``: Both surveys performed with baseline arrays
* ``obs_all_threshold.xml``: Both surveys performed with threshold implementation

The ``models`` folder contains the definitions of all source and background
models that have been used for simulating the data.
The file ``models_gps.xml`` is a so called
:ref:`Model Definition File <glossary_moddef>`
that collects the definition of all model components used for the Galactic
Plane Survey simulation.
The file ``models_egal.xml`` is the equivalent file for the Extragalactic
Survey.
The other files in the folder are FITS files containing sky maps that were
used as spatial templates for modelling extended or diffuse emission components.
There files are reference in the
:ref:`Model Definition File <glossary_moddef>`
files.



