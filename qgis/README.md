PyRATBridge for QGIS
====================

Overview
--------

This Plugin imports most of the functionality of PyRAT into QGIS.

Installation
============

You need a local installation of PyRAT to install PyRATBridge.
There are multiple ways of installing PyRATBridge.

Installation (simple)
---------------------

Compress the PyRATBridge folder into a zip file.
In QGIS open the Plugin Browser (Plugins -> Manage and install plugins).
Click on the tab Install from Zip, select the zipped file and install the plugin.

Installation on Linux (alternative, with anaconda)
--------------------------------------------------

Switch into the PyRATBridge directory and run the setup script:

    ./setup.sh

Then you can start qgis in your (new) conda environment.
You only need to enable the plugin in the QGIS Plugin Browser (Plugins -> Manage and install plugins, tab: Installed, check PyRATBridge)

Installation on Windows (alternative, with anaconda, untested)
-----------------------------------------------

Create a conda environment with the requirements:

    anaconda env create -n ENVIONMENT_NAME -f requirements.yml

It's possible that you need to install pyreadline additionally:

    anaconda install -n ENVIRONMENT_NAME pyreadline

After that, you need to copy the folder pyrat inside the PyRAT repository into the 
environments site-packages folder. You'll find the site-packages folder under:

    C:\Users\YOUR_USER_NAME\AppData\Local\Continuum\anaconda3\envs\ENVIRONMENT_NAME\Lib\site-packages

Zip the contents of this folder (PyRATBridge).

Select in QGIS Plugins -> Manage and install plugins...

Go to the tab 'Install from ZIP' and select the zip file of the tool.

Issues
======

QGIS crashing under Linux in conda environment
----------------------------------------------

PyRAT requires the module readline, which calls a system library.
QGIS calls the system library too, but might use a different version (conda version and system version).
This leads to a segmentation fault.
If you're affected, start QGIS with the following command:

    LD_PRELOAD=PATH_TO_CONDA_ENVIRONMENT/lib/libedit.so qgis

the plugin will itself offer you an alias for this command after a crash on startup.
