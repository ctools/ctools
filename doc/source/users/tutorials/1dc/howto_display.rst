.. _1dc_howto_display:

How to display the results?
---------------------------

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
| ``show_spectrum.py``   | Display a spectrum          |
+------------------------+-----------------------------+
