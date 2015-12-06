.. _ctpsfcube:

ctpsfcube
=========

Generate point spread function cube for a counts cube.


Synopsis
--------

This tool generates a point spread function cube for a counts cube. A point
spread function cube is a 4-dimensional cube spanned by Right Ascension or
Galactic longitude, Declination or Galactic latitude, energy, and offset 
angle between true and measured arrival direction of a photon. The energy
binning of the cube may be either linear, logarithmic, or custom defined
using an input file.

ctpsfcube requires on input the event list or observation definition file 
that has been used in the generation of the counts cube using :doc:`ctbin`.

It is not recommended to use the counts cube for the point spread function 
cube definition, although this is formally possibly by specifying the counts 
cube as ``incube`` parameter. This leads however to a large FITS file on 
output since the number of bins in the counts cube will be multiplied by 
the number of offset angle bins (typically 200). Since the point spread 
function varies only little over the field of fiew of the camera it is 
recommended to use a rather coarse spatial binning to keep the file size 
manageable (with a typical value of ``binsz=1.0``).

ctpsfcube generates a point spread function cube FITS file comprising three
extensions. The primary extension contains a 3-dimensional image that contains
the point spread function values. The energy and offset angle dimensions 
of the point spread function cube are folded into the 3rd dimension of the 
FITS image. The next extension named ``EBOUNDS`` contains a binary table
that defines the energy boundaries of the exposure cube. The last extension
named ``DELTAS`` contains a binary table that defines the offset angles 
between true and measured arrival direction of the photon.


General parameters
------------------

``inobs [file]``
    Input event list or observation definition XML file.

``incube [file]``
    Counts cube for point spread function cube definition.

``outcube [file]``
    Output point spread function cube file.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrument response function.

``(edisp = no) [boolean]``
    Apply energy dispersion for response computation.

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

``(amax = 0.3) [real]``
    Upper bound of angular separation between true and measued photon
    direction (in degrees).

``(anumbins = 200) [integer]``
    Number of angular separation bins.
 	 	 

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
    Specifies whether an existing output counts cube should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = ctpsfcube.log) [string]``
    Name of log file.


Related tools
-------------

:doc:`ctbin`
