.. _howto_display:

How to display the results?
---------------------------

  .. admonition:: What you will learn

     You will learn how to **display results** using some example scripts
     written in Python.

     Please not that these scripts need the ``matplotlib`` Python module
     installed.

The ctools package does not contain any tools or scripts for graphical
display of results since results are generally written into standard FITS
files that are readily displayed by existing astronomical tools.

Nevertheless, for your convenience several scripts for graphical display are
included in the ctools package that rely on the
`matplotlib <http://matplotlib.org>`_
Python module. You can find these scripts in the
``$CTOOLS/share/examples/python`` folder. The following scripts are available:

  +---------------------------------------+--------------------------------------+
  | Script                                | Usage                                |
  +=======================================+======================================+
  | ``show_butterfly.py``                 | Display butterfly diagram            |
  +---------------------------------------+--------------------------------------+
  | ``show_irf.py``                       | Display Instrument Response Function |
  +---------------------------------------+--------------------------------------+
  | ``show_lightcurve.py``                | Display light curve                  |
  +---------------------------------------+--------------------------------------+
  | ``show_obs.py``                       | Display observation summary          |
  +---------------------------------------+--------------------------------------+
  | ``show_pha.py``                       | Display Pulse High Analyser File     |
  +---------------------------------------+--------------------------------------+
  | ``show_phases.py``                    | Display event phases                 |
  +---------------------------------------+--------------------------------------+
  | ``show_pointings.py``                 | Display pointing directions          |
  +---------------------------------------+--------------------------------------+
  | ``show_residuals.py``                 | Display spectral residuals           |
  +---------------------------------------+--------------------------------------+
  | ``show_response.py``                  | Display response components          |
  +---------------------------------------+--------------------------------------+
  | ``show_rmf.py``                       | Display Redistribution Matrix File   |
  +---------------------------------------+--------------------------------------+
  | ``show_significance_distribution.py`` | Display significance distribution    |
  +---------------------------------------+--------------------------------------+
  | ``show_spectrum.py``                  | Display spectrum                     |
  +---------------------------------------+--------------------------------------+

**Do not hesitate to copy and adapt these scripts to your needs.**

Below some usage examples and the expected output.

show_butterfly.py
^^^^^^^^^^^^^^^^^

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_butterfly.py butterfly.txt

.. figure:: howto_display_butterfly.png
   :width: 500px
   :align: center

   *Butterfly diagram displayed with show_butterfly.py*


show_irf.py
^^^^^^^^^^^

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_irf.py prod2 South_50h

.. figure:: howto_display_irf.png
   :width: 600px
   :align: center

   *Instrument Response Function displayed with show_irf.py*


show_lightcurve.py
^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_lightcurve.py lightcurve.fits

.. figure:: howto_lightcurve.png
   :width: 500px
   :align: center

   *Light curve displayed with show_lightcurve.py*


show_obs.py
^^^^^^^^^^^

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_obs.py obs.xml

.. figure:: howto_display_obs.png
   :width: 600px
   :align: center

   *Observation summary displayed with show_obs.py*


show_pha.py
^^^^^^^^^^^

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_pha.py pha.fits

.. figure:: howto_display_pha.png
   :width: 500px
   :align: center

   *Pulse Hight Analyser File displayed with show_pha.py*


show_phases.py
^^^^^^^^^^^^^^

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_phases.py -n 50 events_phased.fits

.. figure:: howto_display_phases.png
   :width: 500px
   :align: center

   *Event phases displayed with show_phases.py*


show_pointings.py
^^^^^^^^^^^^^^^^^

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_pointings.py obs.xml

.. figure:: howto_display_pointings.png
   :width: 500px
   :align: center

   *Pointings displayed with show_pointings.py (zoomed in)*


show_residuals.py
^^^^^^^^^^^^^^^^^

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_residuals.py residual.fits

.. figure:: howto_display_residuals.png
   :width: 500px
   :align: center

   *Spectral residuals displayed with show_residuals.py*


show_response.py
^^^^^^^^^^^^^^^^

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_response.py

.. image:: howto_display_rsp_aeff.png
   :width: 32%
.. image:: howto_display_rsp_bkg.png
   :width: 32%
.. image:: howto_display_rsp_sens.png
   :width: 32%


show_rmf.py
^^^^^^^^^^^

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_rmf.py rmf.fits

.. figure:: howto_display_rmf.png
   :width: 500px
   :align: center

   *Redistribution Matrix File displayed with show_rmf.py*


show_significance_distribution.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_significance_distribution.py skymap.fits --exclude skymap.fits[EXCLUSION]

.. figure:: howto_display_significance.png
   :width: 500px
   :align: center

   *Significance distribution displayed with show_significance_distribution.py*


show_spectrum.py
^^^^^^^^^^^^^^^^

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_spectrum.py spectrum.fits

.. figure:: howto_display_spectrum.png
   :width: 500px
   :align: center

   *Spectrum displayed with show_spectrum.py*
