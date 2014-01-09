Introduction
============

ctools is a highly modular collection of utilities for processing
and analysing CTA reconstructed event data in the FITS (Flexible
Image Transport System) data format. Each utility presents itself
as a FTOOL (see http://heasarc.gsfc.nasa.gov/ftools/) and performs
a single simple task such as event binning, event selection or model
fitting. Individual utilities can easily be chained together in
scripts to achieve more complex operations, either by using the command
line interface, or by using the Python scripting language. The ctools
user interface is controlled by standard IRAF-style parameter files.
Software is written in C++ to provide portability across most computer
systems. The data format dependencies between hardware platforms are
isolated through the cfitsio library package from HEASARC 
(http://heasarc.gsfc.nasa.gov/fitsio/).

This User Manual describes the use of the ctools software.