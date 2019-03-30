.. _about:

About
=====

ctools is a software package developed for the scientific analysis of 
Cherenkov Telescope Array (CTA) data, as well as data from existing
Imaging Air Cherenkov Telescopes (such as H.E.S.S., MAGIC or VERITAS).

ctools comprises a set of ftools-like binary executables and Python scripts
with a command-line interface allowing for interactive step-wise data analysis.
ctools includes also Python modules to control all tools and scripts from
within Python.
Creation of shell or Python scripts and pipelines is supported as well.

ctools are based on GammaLib, a versatile toolbox for the high-level 
analysis of astronomical gamma-ray data. 
Besides CTA, GammaLib supports also the analysis of Fermi/LAT and COMPTEL 
data, and extensions to support further gamma-ray instruments are planned. 
An interface to virtual observatory ressources does also exist.
By making use of the GammaLib multi-instrument capabilities, ctools 
supports the joint analysis of CTA, H.E.S.S., MAGIC, VERITAS, Fermi/LAT and
COMPTEL data.

ctools are developed by a team of enthousiastic gamma-ray astronomers with
support from engineers. We regularily organise
`coding sprints <https://cta-redmine.irap.omp.eu/projects/ctools/wiki/Coding_sprints>`_
where key developers but also newcomers meet to discuss the developments 
and next steps, and advance with the coding of the software.


Acknowledging or citing ctools
------------------------------

If you use ctools for work/research presented in a publication we ask you
to include the formal reference

   J. Knoedlseder, M. Mayer, C. Deil, J.-B. Cayrou, E. Owen, N. Kelley-Hoskins,
   C.-C. Lu, R. Buehler, F. Forest, T. Louge, H. Siejkowski, K. Kosack,
   L. Gerard, A. Schulz, P. Martin, D. Sanchez, S. Ohm, T. Hassan, and
   S. Brau-Nogue, 2016, A&A, 593, A1 (`paper <http://www.aanda.org/articles/aa/pdf/2016/09/aa28822-16.pdf>`_)

in your paper and that you add the following acknowledgment:

   *This research made use of ctools, a community-developed analysis package
   for Imaging Air Cherenkov Telescope data. ctools is based on GammaLib,
   a community-developed toolbox for the high-level analysis of astronomical
   gamma-ray data.*

If you wish, you can also include a link to http://cta.irap.omp.eu/ctools/
(if the journal allows this) in addition to the above text, as well as the
identifier of the code in the Astrophysics Source Code Library (ASCL):
`ascl:1601.005 <http://ascl.net/1601.005>`_.

If you are giving a presentation or talk featuring work/research that makes
use of ctools, we suggest using this logo on your title slide:

.. figure:: ctools-logo.jpg
   :width: 150px
   :align: center

And to see who published an article using ctools you may check the `following link <http://cdsads.u-strasbg.fr/cgi-bin/nph-ref_query?bibcode=2016A%26A...593A...1K&amp;refs=CITATIONS&amp;db_key=AST>`_.


License
-------

ctools is free software distributed under the GNU GPL license version 3.
