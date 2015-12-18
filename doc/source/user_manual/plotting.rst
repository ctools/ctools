Plotting
========

Gammalib and ctools do not include classes and scripts to make plots.

This is intentional because creating and maintaining plotting tools that
work on most user's machines and cover most user's needs is a big project
in itself, and deciding to support a plotting package that is popular now
as part of the CTA science tools could become very problematic in the
future.

But we realise that most ctools users need to create plots for their
high-level science results, so on this page we provide a few pointers 
to external packages and tools.


General
-------

The most common way to create publication-quality astronomy plots at the
moment is to use the `matplotlib`_ Python plotting package. As illustrated
by the `matplotlib gallery`_, virtually any plot can be constructed,
but matplotlib doesn't have built-in easy-to-use classes or functions to
create the common plots for gamma-ray astronomy, which are sky maps and
spectra, so in the sections below we give some specific pointers for those.


Sky maps
--------

For creating publication quality sky map plots via a Python script a good
choice is to use `APLPy`_, `WCSAxes`_ or `kapteyn`_, which are all based on 
matplotlib.

You can also use interactive astronomy image viewers like `DS9`_, `Aladin`_
or `Ginga`_ which allow you to save the image to a PNG or PDF file after
customising it interactively, but also contain command line options to create
the plot, i.e. you can write a bash or Python script to make the plot
reproducible (DS9 and Ginga also feature direct communication with Python,
but that is a bit harder to use and usually used for interactive image
analysis, not plotting scripts).


Spectra and light curves
------------------------

To plot spectra or spectral energy distributions (including flux model curves,
flux model error bands, flux points and residuals) and light curves, you can
use matplotlib directly, or use the classes and scripts in `Gammapy`_.

Another option is to use virtual observatory interactive spectrum tools
like `IRIS`_, `CASSIS`_ or `VOSpec`_ which allow you to save the plot to a
PNG or PDF file after creating it interactively.

.. _DS9: http://ds9.si.edu/
.. _Aladin: http://aladin.u-strasbg.fr/
.. _Ginga: http://ejeschke.github.io/ginga/
.. _matplotlib: http://matplotlib.org/
.. _matplotlib gallery: http://matplotlib.org/gallery.html
.. _WCSAxes: http://wcsaxes.readthedocs.org/
.. _APLPy: http://aplpy.github.io/
.. _kapteyn: https://www.astro.rug.nl/software/kapteyn/
.. _Gammapy: https://gammapy.readthedocs.org/
.. _IRIS: http://cxc.cfa.harvard.edu/iris/
.. _CASSIS: http://cassis.irap.omp.eu/
.. _VOSpec: http://www.sciops.esa.int/index.php?project=SAT&page=vospec
