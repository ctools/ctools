.. _hess_dr1_select_events:

Selecting the relevant events
-----------------------------

  .. admonition:: What you will learn

     You will learn **how to select the relevant events for a target**
     and **how to determine the useful energy range** for a set of
     observations.

Now you are ready for your first H.E.S.S. data analysis.
As first step you should create a folder that will contain your analysis
results and any intermediate data products.
Let's call this folder ``crab`` since we explain in this example how you
can analyse the data from the Crab nebula.
Create and step into the folder using

.. code-block:: bash

   $ mkdir crab
   $ cd crab

The first step of your analysis consists in selecting the relevant events
from the observations.
In this step you can select a specific energy range, time range or region
of interest.
In the example below you will select events according to the save energy
thresholds that are defined in the effective area component of the
:ref:`instrument response functions <glossary_irf>`.
You do this by setting the hidden parameter ``usethres=DEFAULT``.
In addition, you will select all events within 2 degrees of the pointing
direction:

.. code-block:: bash

   $ ctselect usethres=DEFAULT
   Input event list or observation definition XML file [events.fits] $HESSDATA/obs/obs_crab.xml
   Lower energy limit (TeV) [0.1] NONE
   Radius of ROI around pointing or specified RA/DEC (degrees) (0-180) [3.0] 2.0
   Start time (UTC string, JD, MJD or MET in seconds) [NONE]
   Output event list or observation definition XML file [selected_events.fits] obs_crab_selected.xml

Below is an excerpt of the log file that was created by :ref:`ctselect`.
You will find that the following energy ranges were applied for the
observations:

* ``id=023523``: 0.871 - 100 TeV
* ``id=023526``: 0.708 - 100 TeV
* ``id=023559``: 0.661 - 100 TeV
* ``id=023592``: 0.871 - 100 TeV

.. code-block:: none

   2018-09-11T14:05:15: +=================+
   2018-09-11T14:05:15: | Event selection |
   2018-09-11T14:05:15: +=================+
   2018-09-11T14:05:15: === HESS observation "Crab" (id=023523) ===
   2018-09-11T14:05:15:  Input filename ............: /project-data/hess/hess_dl3_dr1/data/hess_dl3_dr1_obs_id_023523.fits.gz
   2018-09-11T14:05:15:  Event extension name ......: EVENTS
   2018-09-11T14:05:15:  GTI extension name ........: GTI
   2018-09-11T14:05:16: === Event selection ===
   2018-09-11T14:05:16:  Selected energy range .....: 0.870964 - 100 TeV
   2018-09-11T14:05:16:  Requested RoI .............: Centre(RA,DEC)=(83.633333333333, 21.514444444444) deg, Radius=2 deg
   2018-09-11T14:05:16:  RoI of data ...............: Centre(RA,DEC)=(83.633333333333, 21.514444444444) deg, Radius=2 deg
   2018-09-11T14:05:16:  Selected RoI ..............: Centre(RA,DEC)=(83.633333333333, 21.514444444444) deg, Radius=2 deg
   2018-09-11T14:05:16:  cfitsio selection .........: ENERGY >= 0.87096359 && ENERGY <= 100.00000000 && ANGSEP(83.633333,21.514444,RA,DEC) <= 2.000000
   2018-09-11T14:05:16:  FITS filename .............: /var/tmp/tmp.0.QX6f0L[EVENTS][ENERGY >= 0.87096359 && ENERGY <= 100.00000000 && ANGSEP(83.633333,21.514444,RA,DEC) <= 2.000000]
   2018-09-11T14:05:16: === HESS observation "Crab" (id=023526) ===
   2018-09-11T14:05:16:  Input filename ............: /project-data/hess/hess_dl3_dr1/data/hess_dl3_dr1_obs_id_023526.fits.gz
   2018-09-11T14:05:16:  Event extension name ......: EVENTS
   2018-09-11T14:05:16:  GTI extension name ........: GTI
   2018-09-11T14:05:16: === Event selection ===
   2018-09-11T14:05:16:  Selected energy range .....: 0.707946 - 100 TeV
   2018-09-11T14:05:16:  Requested RoI .............: Centre(RA,DEC)=(83.633333333333, 22.514444444444) deg, Radius=2 deg
   2018-09-11T14:05:16:  RoI of data ...............: Centre(RA,DEC)=(83.633333333333, 22.514444444444) deg, Radius=2 deg
   2018-09-11T14:05:16:  Selected RoI ..............: Centre(RA,DEC)=(83.633333333333, 22.514444444444) deg, Radius=2 deg
   2018-09-11T14:05:16:  cfitsio selection .........: ENERGY >= 0.70794578 && ENERGY <= 100.00000000 && ANGSEP(83.633333,22.514444,RA,DEC) <= 2.000000
   2018-09-11T14:05:16:  FITS filename .............: /var/tmp/tmp.1.IoZj2b[EVENTS][ENERGY >= 0.70794578 && ENERGY <= 100.00000000 && ANGSEP(83.633333,22.514444,RA,DEC) <= 2.000000]
   2018-09-11T14:05:16: === HESS observation "Crab" (id=023559) ===
   2018-09-11T14:05:16:  Input filename ............: /project-data/hess/hess_dl3_dr1/data/hess_dl3_dr1_obs_id_023559.fits.gz
   2018-09-11T14:05:16:  Event extension name ......: EVENTS
   2018-09-11T14:05:16:  GTI extension name ........: GTI
   2018-09-11T14:05:16: === Event selection ===
   2018-09-11T14:05:16:  Selected energy range .....: 0.660693 - 100 TeV
   2018-09-11T14:05:16:  Requested RoI .............: Centre(RA,DEC)=(85.2533333381014, 22.014444444444) deg, Radius=2 deg
   2018-09-11T14:05:16:  RoI of data ...............: Centre(RA,DEC)=(85.2533333381014, 22.014444444444) deg, Radius=2 deg
   2018-09-11T14:05:16:  Selected RoI ..............: Centre(RA,DEC)=(85.2533333381014, 22.014444444444) deg, Radius=2 deg
   2018-09-11T14:05:16:  cfitsio selection .........: ENERGY >= 0.66069345 && ENERGY <= 100.00000000 && ANGSEP(85.253333,22.014444,RA,DEC) <= 2.000000
   2018-09-11T14:05:16:  FITS filename .............: /var/tmp/tmp.2.IAwBNH[EVENTS][ENERGY >= 0.66069345 && ENERGY <= 100.00000000 && ANGSEP(85.253333,22.014444,RA,DEC) <= 2.000000]
   2018-09-11T14:05:16: === HESS observation "Crab" (id=023592) ===
   2018-09-11T14:05:16:  Input filename ............: /project-data/hess/hess_dl3_dr1/data/hess_dl3_dr1_obs_id_023592.fits.gz
   2018-09-11T14:05:16:  Event extension name ......: EVENTS
   2018-09-11T14:05:16:  GTI extension name ........: GTI
   2018-09-11T14:05:16: === Event selection ===
   2018-09-11T14:05:16:  Selected energy range .....: 0.870964 - 100 TeV
   2018-09-11T14:05:16:  Requested RoI .............: Centre(RA,DEC)=(82.0133333285646, 22.014444444444) deg, Radius=2 deg
   2018-09-11T14:05:16:  RoI of data ...............: Centre(RA,DEC)=(82.0133333285646, 22.014444444444) deg, Radius=2 deg
   2018-09-11T14:05:16:  Selected RoI ..............: Centre(RA,DEC)=(82.0133333285646, 22.014444444444) deg, Radius=2 deg
   2018-09-11T14:05:16:  cfitsio selection .........: ENERGY >= 0.87096359 && ENERGY <= 100.00000000 && ANGSEP(82.013333,22.014444,RA,DEC) <= 2.000000
   2018-09-11T14:05:16:  FITS filename .............: /var/tmp/tmp.3.879wq4[EVENTS][ENERGY >= 0.87096359 && ENERGY <= 100.00000000 && ANGSEP(82.013333,22.014444,RA,DEC) <= 2.000000]

You may use the :ref:`csobsinfo` script to display a summary of the selected
events into the console:

.. code-block:: bash

   $ csobsinfo debug=yes
   Input event list, counts cube, or observation definition XML file [$HESSDATA/obs/obs_crab.xml] obs_crab_selected.xml
   Output DS9 region file [ds9.reg]

.. note::
   The ``debug=yes`` attribute instructs :ref:`csobsinfo` to direct the log
   file output also to the console. Duplication of log file output into
   the console using the ``debug=yes`` attribute works for any tool or script.

The output of :ref:`csobsinfo` is shown below.
:ref:`ctselect` selected 9675 events within an energy interval of
0.66 - 100 TeV.
This will be the energy range that you will from now on use for the data
analysis.

.. code-block:: none

   2018-09-11T20:30:16: +=========+
   2018-09-11T20:30:16: | Summary |
   2018-09-11T20:30:16: +=========+
   2018-09-11T20:30:16: === Observations ===
   2018-09-11T20:30:16:  Unbinned observations .....: 4
   2018-09-11T20:30:16:  Binned observations .......: 0
   2018-09-11T20:30:16: === Events ===
   2018-09-11T20:30:16:  Number of events ..........: 9675
   2018-09-11T20:30:16:  Number of bins ............: 0
   2018-09-11T20:30:16: === Pointings ===
   2018-09-11T20:30:16:  Mean offset angle .........: Unknown
   2018-09-11T20:30:16:  Mean zenith angle .........: 47.04 deg
   2018-09-11T20:30:16:  Mean azimuth angle ........: 13.76 deg
   2018-09-11T20:30:16: === Energy range ===
   2018-09-11T20:30:16:  Minimum energy ............: 660.693448007596 GeV
   2018-09-11T20:30:16:  Maximum energy ............: 100 TeV
   2018-09-11T20:30:16: === Time range ===
   2018-09-11T20:30:16:  MJD (days) ................: 53343.922 - 53347.933
   2018-09-11T20:30:16:  UTC .......................: 2004-12-04T22:07:06 - 2004-12-08T22:22:02
   2018-09-11T20:30:16:  Total ontime ..............: 6742.00 s = 112.37 min = 1.87 h
   2018-09-11T20:30:16:  Total livetime ............: 6313.81 s = 105.23 min = 1.75 h
