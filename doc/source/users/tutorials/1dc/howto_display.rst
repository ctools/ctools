.. _1dc_howto_display:

How to display the results?
---------------------------

.. admonition:: You will learn ...

   ... how to **display results**.

The ctools package does not contain any tools or scripts for graphical
display of results since results are generally written into standard FITS
files that are readily displayed by existing astronomical tools.

Nevertheless, for your convenience several scripts for graphical display are
included in the ctools package that rely on the
`matplotlib <http://matplotlib.org>`_
Python module. You can find these scripts in the
``$CTOOLS/share/examples/python`` folder. The following scripts are available:

+------------------------+-----------------------------+
| Script                 | Usage                       |
+========================+=============================+
| ``show_butterfly.py``  | Display a butterfly diagram |
+------------------------+-----------------------------+
| ``show_lightcurve.py`` | Display a light curve       |
+------------------------+-----------------------------+
| ``show_obs.py``        | Display observation summary |
+------------------------+-----------------------------+
| ``show_pointings.py``  | Display pointing directions |
+------------------------+-----------------------------+
| ``show_spectrum.py``   | Display a spectrum          |
+------------------------+-----------------------------+

Below some usage examples and the expected output:

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_butterfly.py butterfly_src001.txt

.. figure:: howto_display_butterfly.png
   :width: 600px
   :align: center

   *Butterfly diagram displayed with show_butterfly.py*

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_obs.py obs_selected.xml

.. figure:: howto_display_obs.png
   :width: 600px
   :align: center

   *Observation summary displayed with show_obs.py*

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_pointings.py obs_selected.xml

.. figure:: howto_display_pointings.png
   :width: 600px
   :align: center

   *Observation summary displayed with show_pointings.py (zoomed in)*

.. code-block:: bash

   $ $CTOOLS/share/examples/python/show_spectrum.py spectrum_src001.fits

.. figure:: howto_display_spectrum.png
   :width: 600px
   :align: center

   *Observation summary displayed with show_spectrum.py*


