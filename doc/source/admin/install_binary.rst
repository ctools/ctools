.. _sec_install_binary:

Installing as binary package
============================

Mac OS X
--------

After :ref:`downloading the Mac OS X binary disk image <sec_download>`
double-click on the image to mount the disk

.. image:: macosx-disk-image.jpg
   :height: 300px
   :alt: Mac OS X binary disk image
   :align: center

Now open the installer package by right-clicking (ctrl-click) on the installer
package icon. This allows to install a software package with an invalid
certificate (ctools has so far no Apple certificate)

.. image:: macosx-open-package.jpg
   :height: 300px
   :alt: Open installer package
   :align: center

A window will show up, warning that the package comes from an unidentified
developer. Click on ``Open`` to continue.

.. image:: macosx-unidentified.jpg
   :height: 200px
   :alt: Unidentified developer warning
   :align: center

Now the installer opens, warning that the package has been signed with an
invalid certificate. Click on ``Continue`` to proceed.

.. image:: macosx-invalid-cert.jpg
   :height: 300px
   :alt: Invalid certificate warning
   :align: center

Follow now the installer instructions by clicking on ``Continue``. As the
software will be installed in ``/usr/local/gamma`` the installer will ask
you for the ``root`` password on your machine.

.. image:: macosx-installer.jpg
   :height: 300px
   :alt: Installer introduction
   :align: center

At the end you should see a window like this:

.. image:: macosx-installed.jpg
   :height: 300px
   :alt: Installer summary
   :align: center

