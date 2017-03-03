.. _1dc_first_select_events:

Selecting the relevant events
-----------------------------

As next step, you have to select the relevant events from the selected
observations. You do this with the :ref:`ctselect` tool by typing:

.. code-block:: bash

   $ ctselect
   Input event list or observation definition XML file [events.fits] obs.xml
   RA for ROI centre (degrees) (0-360) [83.63] UNDEFINED
   Start time (CTA MET in seconds) [0.0] UNDEFINED
   Lower energy limit (TeV) [0.1]
   Upper energy limit (TeV) [100.0]
   Output event list or observation definition XML file [selected_events.fits] obs_selected.xml

In this example you selected events with energies comprised between 100 GeV
and 100 TeV from the observations. No selection was applied on the spatial
or temporal parameters of the events.

:ref:`ctselect` will create new event list files in the working directory
that contain the selected events. In addition, :ref:`ctselect` creates a new
:ref:`observation definition file <glossary_obsdef>`
``obs_selected.xml`` that contains the metadata of the new event list files,
and that will be used for the further analysis.
The content of the file will look similar to the content of ``obs.xml`` with
the original event list file names replaced by the names of the new files:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <observation_list title="observation list">
     <observation name="GC" id="000414" instrument="CTA">
       <parameter name="EventList" file="selected_gc_baseline_000414.fits" />
       <parameter name="Calibration" database="prod3b" response="South_z20_50h" />
     </observation>
     <observation name="GC" id="000415" instrument="CTA">
       <parameter name="EventList" file="selected_gc_baseline_000415.fits" />
       <parameter name="Calibration" database="prod3b" response="South_z20_50h" />
     </observation>
     <observation name="GC" id="000416" instrument="CTA">
       <parameter name="EventList" file="selected_gc_baseline_000416.fits" />
       <parameter name="Calibration" database="prod3b" response="South_z20_50h" />
     </observation>
     ...
     <observation name="GC" id="000581" instrument="CTA">
       <parameter name="EventList" file="selected_gc_baseline_000581.fits" />
       <parameter name="Calibration" database="prod3b" response="South_z20_50h" />
     </observation>
   </observation_list>
