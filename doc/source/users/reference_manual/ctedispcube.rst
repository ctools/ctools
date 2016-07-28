.. _ctedispcube:

ctedispcube
===========

Generate energy dispersion cube for a stacked analysis.


Synopsis
--------

This tool generates an energy dispersion cube for a stacked analysis. An
energy dispersion cube is a 4-dimensional cube spanned by Right Ascension or
Galactic longitude, Declination or Galactic latitude, energy, and migration
which is the ratio between reconstructed and true photon energy. The energy
binning of the cube may be either linear, logarithmic, or custom defined
using an input file

ctedispcube requires on input the event list or observation definition XML
file that has been used in the generation of the counts cube using :doc:`ctbin`.

It is not recommended to use the counts cube for the energy dispersion 
cube definition, although this is formally possible by specifying the counts 
cube as ``incube`` parameter. This leads however to a large FITS file on 
output since the number of bins in the counts cube will be multiplied by 
the number of migration (typically 100). Since the energy dispersion varies
only little over the field of view of the camera it is recommended to use a
rather coarse spatial binning to keep the file size manageable (with a typical
value of ``binsz=1.0``).

ctedispcube generates an energy dispersion cube FITS file comprising three
extensions. The primary extension contains a 4-dimensional image that contains
the energy disperison values. The next extension named ``ENERGIES`` contains
a binary table that defines the energies of the energy dispersion cube. The
last extension named ``MIGRAS`` contains a binary table that defines the
migration values of the cube which is the ratio between reconstructed and true
photon energy.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``incube [file]``
    Input counts cube file to extract energy dispersion cube definition.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrument response function.

``outcube [file]``
    Output energy dispersion cube file.

``ebinalg <FILE|LIN|LOG> [string]``
    Algorithm for defining energy bins.
 	 	 
``emin [real]``
    Lower energy value for first energy bin (in TeV).
 	 	 
``emax [real]``
    Upper energy value for last energy bin (in TeV).
 	 	 
``enumbins [integer]``
    Number of energy bins.
 	 	 
``ebinfile [file]``
    Name of the file containing the energy bin definition.
 	 	 
``(addbounds = no) [boolean]``
    Add energies to the energy dispersion cube at the observation energy boundaries?

``(usepnt = no) [boolean]``
    Use CTA pointing direction for cube centre instead of xref/yref parameters?
 	 	 
``nxpix [integer]``
    Number of cube bins in Right Ascension or Galactic longitude.
 	 	 
``nypix [integer]``
    Number of cube bins in Declination or Galactic latitude.
 	 	 
``binsz [real]``
    Cube bin size (in degrees/pixel).
 	 	 
``coordsys <CEL|GAL> [string]``
    Coordinate system (CEL - celestial, GAL - galactic).
 	 	 
``xref [real]``
    Right Ascension / Galactic longitude of cube centre (J2000, in degrees).
 	 	 
``yref [real]``
    Declination / Galactic latitude of cube centre (J2000, in degrees).
 	 	 
``proj <AIT|AZP|CAR|MER|MOL|STG|TAN> [string]``
    Projection method.

``(migramax = 2.0) [real]``
    Upper bound of ratio between reconstructed and true photon energy.

``(migrabins = 100) [integer]``
    Number of migration bins.
 	 	 

Standard parameters
-------------------

``(chatter = 2) [integer]``
    Verbosity of the executable:
     chatter = 0: no information will be logged
     
     chatter = 1: only errors will be logged
     
     chatter = 2: errors and actions will be logged
     
     chatter = 3: report about the task execution
     
     chatter = 4: detailed report about the task execution
 	 	 
``(clobber = yes) [boolean]``
    Specifies whether an existing output energy dispersion cube file should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctedispcube.log) [string]``
    Name of log file.


Related tools or scripts
------------------------

:doc:`ctbin`
