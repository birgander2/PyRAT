#!/bin/bash

echo ""
echo "This script installs a conda environment with qgis and pyrat"
read -p "Name of the (new) conda environment to install PyRAT and QGIS: " envname

echo ""
echo "To install PyRATBridge, you need PyRAT version 0.6 or above."
echo "WARNING: Please type a absolute path (beginning with /)!"
read -p "Where is the PyRAT directory located (with the python-module folder 'pyrat' inside)? " PyRAT

echo ""
echo "Starting installation..."
echo ""

if conda env create --file requirements.yml -n $envname; then
    echo "Environment sucessfully created."
else
    echo "Environment already exists, updating it:"
    conda env update -f=requirements.yml -n $envname
fi

#Copy/Install the plugin
mkdir -p ~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/PyRATBridge
cp -r . ~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/PyRATBridge

cp -r $PyRAT/pyrat ~/.conda/envs/$envname/lib/python3.6/site-packages/pyrat
