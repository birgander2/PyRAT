#!/bin/bash

echo ""
echo "This script installs a conda environment with qgis and pyrat"
read -p "Name of the (new) conda environment to install PyRAT and QGIS: " envname

echo ""
echo "To install PyRATBridge, you need PyRAT version 0.6 or above."
echo ""
echo "Starting installation..."
echo ""

if conda env create --file requirements.yml -n "$envname"; then
    echo "Environment sucessfully created."
else
    echo "Environment already exists, updating it:"
    conda env update -f=requirements.yml -n "$envname"
fi

#Copy/Install the plugin
mkdir -p ~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/PyRATBridge
cp -r . ~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/PyRATBridge

curdir=$(dirname "$0")
envpath=$(conda env list | grep qgis | cut -d' ' -f2- | sed 's/^ *//g')

cp -r "$curdir"/../../pyrat "$envpath"/lib/python3.6/site-packages/pyrat
