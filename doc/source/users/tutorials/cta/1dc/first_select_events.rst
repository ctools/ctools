.. _1dc_first_select_events:

Selecting the relevant events
-----------------------------

  .. admonition:: What you will learn

     You will learn how to **select a subset of events** from the selected
     observations for your analysis.

     The event selection is necessary if you want to restrict for example the
     energy range of the events or the time interval that should be covered,
     but also if you want to select only events within a given radial
     acceptance region (also called region of interest).

As next step, you have to select the relevant events from the selected
observations. You do this with the :ref:`ctselect` tool by typing:

.. code-block:: bash

   $ ctselect
   Input event list or observation definition XML file [events.fits] obs.xml
   Radius of ROI around pointing or specified RA/DEC (degrees) (0-180) [3.0] NONE
   Start time (UTC string, JD, MJD or MET in seconds) [NONE]
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Output event list or observation definition XML file [selected_events.fits] obs_selected.xml

In this example you selected events with energies comprised between 100 GeV
and 100 TeV from the observations. No selection was applied on the spatial
or temporal parameters of the events. :ref:`ctselect` will create new event
list files in the working directory that contain the selected events. In
addition, :ref:`ctselect` creates a new
:ref:`observation definition file <glossary_obsdef>`
``obs_selected.xml`` that contains the metadata of the new event list files,
and that will be used for the further analysis. Running

.. code-block:: bash

   $ csobsinfo debug=yes
   Input event list, counts cube, or observation definition XML file [obs.xml] obs_selected.xml
   Output DS9 region file [ds9.reg]

displays a summary of the content of this new
:ref:`observation definition file <glossary_obsdef>`:

.. code-block:: none

   2018-01-26T10:33:15: +=========+
   2018-01-26T10:33:15: | Summary |
   2018-01-26T10:33:15: +=========+
   2018-01-26T10:33:15: === Observations ===
   2018-01-26T10:33:15:  Unbinned observations .....: 57
   2018-01-26T10:33:15:  Binned observations .......: 0
   2018-01-26T10:33:15: === Events ===
   2018-01-26T10:33:15:  Number of events ..........: 3084595
   2018-01-26T10:33:15:  Number of bins ............: 0
   2018-01-26T10:33:15: === Pointings ===
   2018-01-26T10:33:15:  Mean offset angle .........: Unknown
   2018-01-26T10:33:15:  Mean zenith angle .........: 0.00 deg
   2018-01-26T10:33:15:  Mean azimuth angle ........: 0.00 deg
   2018-01-26T10:33:15: === Energy range ===
   2018-01-26T10:33:15:  Minimum energy ............: 100 GeV
   2018-01-26T10:33:15:  Maximum energy ............: 100 TeV
   2018-01-26T10:33:15: === Time range ===
   2018-01-26T10:33:15:  MJD (days) ................: 59235.500 - 59276.921
   2018-01-26T10:33:15:  UTC .......................: 2021-01-21T11:58:51 - 2021-03-03T22:04:51
   2018-01-26T10:33:15:  MET (seconds) .............: 664502400.000 - 668081160.000
   2018-01-26T10:33:15:  Total ontime ..............: 102600.00 s = 1710.00 min = 28.50 h
   2018-01-26T10:33:15:  Total livetime ............: 100548.00 s = 1675.80 min = 27.93 h

There are still 57 observations in the
:ref:`observation definition file <glossary_obsdef>`
but now the events are restricted to the energy interval 100 GeV - 100 TeV. In
total there are 3 084 595 events within this interval in the selected observations.

The content of ``obs_selected.xml`` will look similar to the content of
``obs.xml`` with the original event list file names replaced by the names of
the new event files:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <observation_list title="observation list">
     <observation name="GPS" id="110380" instrument="CTA">
       <parameter name="EventList" file="selected_gps_baseline_110380.fits" />
       <parameter name="Calibration" database="1dc" response="South_z20_50h" />
     </observation>
     <observation name="GPS" id="110381" instrument="CTA">
       <parameter name="EventList" file="selected_gps_baseline_110381.fits" />
       <parameter name="Calibration" database="1dc" response="South_z20_50h" />
     </observation>
     ...
     <observation name="GPS" id="111177" instrument="CTA">
       <parameter name="EventList" file="selected_gps_baseline_111177.fits" />
       <parameter name="Calibration" database="1dc" response="South_z20_50h" />
     </observation>
   </observation_list>
