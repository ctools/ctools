.. _about:

About
=====

ctools is a software package for the scientific analysis of astronomical
gamma-ray data. The package comprises an extensive set of tools for the analysis
of data from existing and future Cherenkov telescopes, including H.E.S.S.,
VERITAS, MAGIC and CTA. ctools supports also the analysis of data from
CGRO/COMPTEL, Fermi/LAT and INTEGRAL/SPI, enabling the exploration of the full
gamma-ray energy band, spanning from hundreds of keV to hundreds of TeV.

The ctools philosophy inherited from ftools consists of providing building
blocks that perform well-defined science data analysis tasks, including
observation and event selection, binning, sky map creation, source detection,
model fitting, spectra, phase curve and light curve generation, and observation
simulations. The building blocks are then combined by the user to create
flexible analysis workflows that can fit the needs of any scientist.

ctools can be used as command-line executables alike ftools to carry out simple
analyses even for users without any programming skills. Furthermore, the tools
are accessible through dedicated Python modules providing an alternative user
interface to the software, as well as the possibility to provide tutorials in
the form of Jupyter notebooks. Shell or Python scripts can be used by more
advanced users to build analysis pipelines to any degrees of complexity.

ctools are based on `GammaLib <http://cta.irap.omp.eu/gammalib>`_, a versatile
toolbox for the scientific analysis of astronomical gamma-ray data.

ctools are developed by a team of enthousiastic gamma-ray astronomers with
support from engineers.


Acknowledging or citing ctools
------------------------------

If you use ctools for work/research presented in a publication we ask you
to include the formal reference

   `J. Knoedlseder, M. Mayer, C. Deil, J.-B. Cayrou, E. Owen, N. Kelley-Hoskins,
   C.-C. Lu, R. Buehler, F. Forest, T. Louge, H. Siejkowski, K. Kosack,
   L. Gerard, A. Schulz, P. Martin, D. Sanchez, S. Ohm, T. Hassan, and
   S. Brau-Nogue, 2016, A&A, 593, A1 <https://www.aanda.org/articles/aa/pdf/2016/09/aa28822-16.pdf>`_

in your paper as well as the Astrophysics Source Code Library (ASCL) identifiers

.. image:: https://img.shields.io/badge/ascl-1601.005-blue.svg?colorB=262255
   :target: http://ascl.net/1601.005

for ctools and

.. image:: https://img.shields.io/badge/ascl-1110.007-blue.svg?colorB=262255
   :target: http://ascl.net/1110.007

for GammaLib. You may also reference the softwares by their Digital Object
Identifiers (DOI) on Zenodo, which are

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3265423.svg
   :target: https://doi.org/10.5281/zenodo.3265423

for ctools and

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3265404.svg
   :target: https://doi.org/10.5281/zenodo.3265404

for GammaLib.

In addition please add the following acknowledgment:

   *This research made use of ctools, a community-developed gamma-ray astronomy
   science analysis software. ctools is based on GammaLib, a community-developed
   toolbox for the scientific analysis of astronomical gamma-ray data.*

If the journal allows this, you can also include a link to
http://cta.irap.omp.eu/ctools/ in addition to the above text.

And to see who published an article using ctools you may check the `following link <http://cdsads.u-strasbg.fr/cgi-bin/nph-ref_query?bibcode=2016A%26A...593A...1K&amp;refs=CITATIONS&amp;db_key=AST>`_.


License
-------

ctools is free software distributed under the GNU GPL license version 3.
