.. _csbkgmodel:

csbkgmodel
==========

Generates background model for 3D analysis.


Synopsis
--------

This script generates a background model for a 3D analysis. The background
model is composed of spatial and spectral components, controlled through the
parameters ``spatial`` and ``spectral``. The spatial component can be either
an effective area model (``AEFF``), a template extracted from the instrument
response function (``IRF``) or an analytical 2D Gaussian function (``GAUSS``).
For the analytical 2D Gaussian function, a bilinear multiplicative background
gradient can be added using ``gradient=yes``. The spectral component can be
either a power law (``PLAW``) or a nodes function (``NODES``). The number of
nodes for the node function is controlled through the ``enumbins`` parameter.
Specific node energies can be specified through the ``ebinalg`` parameter.

If ``runwise=yes``, a background model component will be added for each observation
in the input observation definition XML file. The script fits the background
model to the input observation(s) to preset the background model parameters
with resonable values.

On output, :ref:`csbkgmodel` writes an model definition XML file that can be
used for model fitting. Source model components have to be added as needed.


General parameters
------------------

``inobs [file]``
    Event list, counts cube or input observation definition XML file.

``caldb [string]``
    Calibration database.

``irf [string]``
    Instrument response function.

``outmodel [file]``
    Output model XML file.

``instrument [string]``
    Instrument name for which background models should be generated. If the
    input observation definition XML contains only observations for a given
    instrument, the instrument name will be extracted from the observation
    definition XML file.

``spatial <IRF|AEFF|GAUSS> [string]``
    Spatial model component. So far IRF template, effective area and radial
    Gaussian models are supported.

``gradient [boolean]``
    Allow for a spatial gradient in the background event distribution?
    This option only applies to the ``GAUSS`` spatial model component.

``spectral <PLAW|NODES> [string]``
    Spectral model component. ``PLAW`` specifies a simple power law model,
    ``NODES`` specifies a piecewise broken power-law.

``ebinalg <FILE|LIN|LOG> [string]``
    Algorithm for defining energy nodes. For ``FILE``, the energy nodes are
    defined in a FITS file that is specified by the ``ebinfile`` parameter,
    for ``LIN`` and ``LOG`` there will be ``enumbins`` energy nodes spaced
    linearly or logarithmically between ``emin`` and ``emax``, respectively.

``emin [real]``
    Lower energy limit (in TeV).

``emax [real]``
    Upper energy limit (in TeV).

``enumbins [integer]``
    Number of energy nodes if ``LIN`` or ``LOG`` energy algorithms are used.

``ebinfile [file]``
    Name of the file containing the energy node definition if ``ebinalg=FILE``.

``runwise [boolean]``
    Generate runwise background model? If ``yes`` is specified a background
    model component for each observation in the observation definition XML
    file will be added to the output model XML file.

``(rad = 2.0) [real]``
    Radius for event selection. Only events within this radius around the
    pointing direction will be used for model fitting.


Standard parameters
-------------------

``(chatter = 2) [integer]``
    Verbosity of the executable:
     ``chatter = 0``: no information will be logged

     ``chatter = 1``: only errors will be logged

     ``chatter = 2``: errors and actions will be logged

     ``chatter = 3``: report about the task execution

     ``chatter = 4``: detailed report about the task execution

``(clobber = yes) [boolean]``
    Specifies whether an existing energy boundaries output file should be overwritten.

``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.

``(mode = ql) [string]``
    Mode of automatic parameters (default is ``ql``, i.e. "query and learn").

``(logfile = csbkgmodel.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`ctlike`
:doc:`csmodelmerge`
