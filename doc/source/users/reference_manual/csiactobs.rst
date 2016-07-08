.. _csiactobs:

csiactobs
===========

Creates an observation XML file for IACTs using a list of observation IDs as
input.


Synopsis
--------

This script creates an observation XML file from a list of observation IDs
provided in an ASCII file. Such an ASCII file can for instance be created with
:doc:`csfindobs`. The script operates on IACT data that needs to be present on
the user machine. `See here <http://gamma-astro-data-formats.readthedocs.org/en/latest/index.html>`__
how IACT data should be structured.

The parameter ``datapath`` is only queried if the environment variable $VHEFITS
is not set.

Besides creating an observation XML file, this script also stores a model XML
file that contains a background model for each observation. The number of free
parameters per background model is steered via the ``bkgpars`` parameter. 

The user can also specify as ``inmodel`` parameter an XML model file with e.g.
sky models. The output model will then contain these models as well.

The script further allows to specify hierarchies for the data formats and
background models. The script will pick the first available option and fall back
to the next one if the first choice is unavailable. If e.g. the IRF background
model is not available, the script will assign the next background model in the
hierarchy ``bkg_mod_hiera``. If the user wants to use a Gaussian background model
as first choice, the parameter ``bkg_mod_hiera`` should be specified to "gauss". 
In addition, the user can specify start parameters for the background parameters
that will be in the output XML model. 

If an observation that given in the input ASCII file is not available on the user
machine, the script dumps a warning into the logfile and the observation is
skipped. 
 

General parameters
------------------

``datapath [string]``
    Path were data are located.

``prodname [string]``
    Name of FITS production (Run :doc:`csiactdata` to view your options).
    
``infile [file]``
    Input runlist ASCII file.

``(inmodel = NONE) [file]``
    Input model XML file.
    
``outobs [file]``
    Output Observation Definition XML file.
    
``outmodel [file]``
    Output model XML file.

``bkgpars [integer]``
    Number of free parameters for each background model.

``(master_indx = master.json) [file]``
    Name of master index file.
    
``(bkg_scale = yes) [boolean]``
    Specifies whether the background scaling factor from the observation index
    file should be applied if available. 

``(ev_hiera = events) [string]``
    Hierarchy of event file formats.

``(aeff_hiera = aeff_2d) [string]``
    Hierarchy of effective area formats.

``(psf_hiera = psf_king|psf_3gauss) [string]``
    Hierarchy of point spread function formats.

``(edisp_hiera = edisp_2d) [string]``
    Hierarchy of energy disperson formats.
    
``(bkg_hiera = bkg_3d) [string]``
    Hierarchy of background formats.
 
``(bkg_mod_hiera = irf|aeff|gauss) [string]``
    Hierarchy of background models.

``(bkg_gauss_norm = 1e-8) [real]``
    Input normalisation for Gaussian background.

``(bkg_gauss_index = -2.0) [real]``
    Input spectral index for Gaussian background.

``(bkg_gauss_sigma = 1e-8) [real]``
    Input sigma for Gaussian background.

``(bkg_aeff_norm = 1e-14) [real]``
    Input normalisation for effective area background.

``(bkg_aeff_index = -2.0) [real]``
    Input spectral index for effective area background.

``(bkg_range_factor = 100.0) [real]``
    Factor to determine range of background normalisation.

    
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
    Specifies whether existing output files should be overwritten.
 	 	 
``(debug = no) [boolean]``
    Enables debug mode. In debug mode the executable will dump any log file output to the console.
 	 	 
``(mode = ql) [string]``
    Mode of automatic parameters (default is "ql", i.e. "query and learn").

``(logfile = csiactobs.log) [filename]``
    Log filename.


Related tools or scripts
------------------------

:doc:`csiactobs`