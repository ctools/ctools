.. _hess_dr1_data:

Getting the data
================

The data for the H.E.S.S. first public release can be downloaded
from `this web site <https://www.mpi-hd.mpg.de/hfm/HESS/pages/dl3-dr1>`_.
The web site also contains release notes for the data and a link to the
description of the data formats.

Once downloaded, you should have a tarball named ``hess_dl3_dr1.tar.gz``
on your disk.

Now uncompress the tarball at any place by typing

.. code-block:: bash

   $ tar xfvz hess_dl3_dr1.tar.gz

You should now have a folder named ``hess_dl3_dr1`` in your current working
directory with the following structure:

.. code-block:: bash

   hess_dl3_dr1/
   hess_dl3_dr1/data
   hess_dl3_dr1/hdu-index.fits.gz
   hess_dl3_dr1/obs-index.fits.gz
   hess_dl3_dr1/README.txt

You will only need the files in the ``data`` directory. Each file corresponds
to one H.E.S.S observing run, which typically lasts about 28 minutes. The
folder contains observations of the Crab nebula, PKS 2155-304, MSH 15-52 and
RX J1713.7-3946. In addition, there are observations corresponding to empty
fields that can be used to estimate the instrumental background.

In addition to the data files you will need
:ref:`observation definition XML files <glossary_obsdef>`
that comprise all observations for a given target.
Please download the following tarball
:download:`hess_dl3_dr1_obs.tar.gz`
and uncompress it by typing

.. code-block:: bash

   $ tar xfvz hess_dl3_dr1_obs.tar.gz

You should now have a folder named ``obs`` in your current working directory
with the following content:

.. code-block:: bash

   hess_dl3_dr1/obs/obs_crab.xml
   hess_dl3_dr1/obs/obs_msh.xml
   hess_dl3_dr1/obs/obs_off.xml
   hess_dl3_dr1/obs/obs_pks.xml
   hess_dl3_dr1/obs/obs_pks_flare.xml
   hess_dl3_dr1/obs/obs_pks_steady.xml
   hess_dl3_dr1/obs/obs_rx.xml

Before continuing, please set the environment variable ``HESSDATA`` to the
``hess_dl3_dr1`` directory, e.g.

.. code-block:: bash

   $ export HESSDATA=$PWD/hess_dl3_dr1

To check whether the data are properly installed you can use the
:ref:`csobsinfo` script:

.. code-block:: bash

   $ csobsinfo debug=yes
   Input event list, counts cube, or observation definition XML file [obs.xml] $HESSDATA/obs/obs_crab.xml
   Output DS9 region file [ds9.reg]
   2018-09-11T13:30:28: +============+
   2018-09-11T13:30:28: | Parameters |
   2018-09-11T13:30:28: +============+
   2018-09-11T13:30:28:  inobs .....................: $HESSDATA/obs/obs_crab.xml
   2018-09-11T13:30:28:  outds9file ................: ds9.reg
   2018-09-11T13:30:28:  offset ....................: no
   2018-09-11T13:30:28:  ra ........................: [not queried]
   2018-09-11T13:30:28:  dec .......................: [not queried]
   2018-09-11T13:30:28:  chatter ...................: 2
   2018-09-11T13:30:28:  clobber ...................: yes
   2018-09-11T13:30:28:  debug .....................: yes
   2018-09-11T13:30:28:  mode ......................: ql
   2018-09-11T13:30:28:  logfile ...................: csobsinfo.log
   2018-09-11T13:30:28:
   2018-09-11T13:30:28: +==============+
   2018-09-11T13:30:28: | Observations |
   2018-09-11T13:30:28: +==============+
   2018-09-11T13:30:28:
   2018-09-11T13:30:28: +=========+
   2018-09-11T13:30:28: | Summary |
   2018-09-11T13:30:28: +=========+
   2018-09-11T13:30:28: === Observations ===
   2018-09-11T13:30:28:  Unbinned observations .....: 4
   2018-09-11T13:30:28:  Binned observations .......: 0
   2018-09-11T13:30:28: === Events ===
   2018-09-11T13:30:28:  Number of events ..........: 30129
   2018-09-11T13:30:28:  Number of bins ............: 0
   2018-09-11T13:30:28: === Pointings ===
   2018-09-11T13:30:28:  Mean offset angle .........: Unknown
   2018-09-11T13:30:28:  Mean zenith angle .........: 47.04 deg
   2018-09-11T13:30:28:  Mean azimuth angle ........: 13.76 deg
   2018-09-11T13:30:28: === Energy range ===
   2018-09-11T13:30:28:  Minimum energy ............: undefined
   2018-09-11T13:30:28:  Maximum energy ............: undefined
   2018-09-11T13:30:28: === Time range ===
   2018-09-11T13:30:28:  MJD (days) ................: 53343.922 - 53347.933
   2018-09-11T13:30:28:  UTC .......................: 2004-12-04T22:07:06 - 2004-12-08T22:22:02
   2018-09-11T13:30:28:  Total ontime ..............: 6742.00 s = 112.37 min = 1.87 h
   2018-09-11T13:30:28:  Total livetime ............: 6313.81 s = 105.23 min = 1.75 h
   2018-09-11T13:30:28:
   2018-09-11T13:30:28: +============================+
   2018-09-11T13:30:28: | Save pointings in DS9 file |
   2018-09-11T13:30:28: +============================+
   2018-09-11T13:30:28:  DS9 filename ..............: ds9.reg
   2018-09-11T13:30:28:
   2018-09-11T13:30:28: Application "csobsinfo" terminated after 11 wall clock seconds, consuming 0.169556 seconds of CPU time.

