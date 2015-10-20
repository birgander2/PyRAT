# PyRAT (Python Radar Analysis Tools)

PyRat is a flexible framework for postprocessing synthetic aperture radar (SAR) data. It
is made for both airborne and spaceborne data and especially focused on providing an
easy plugin-based programming interface. 

Technically, PyRat is implemented in Python (supported by some Cython) and uses HDF5 based 
disc containers for temporary storage. It features automatic multithreaded block 
processing for speed and memory efficiency, a powerful batch system and a Qt-based GUI. 
It is expandable by plugins without deep knowledge of framework itself.

## Installation

Compile / build:
    
    python setup.py build_ext --inplace

Install (with root rights)

    python setup.py install

Install (as user)

    python setup.py install --user

## Usage

CLI Interface:

    ./pyrat.py -b [rat filename]
    ./pyrat.py --batch [rat filename]
    
GUI Interface:

    ./pyrat.py [rat filename]

Current modules:
* load:      Importing of data
* save:      Exporting of data
* filter:    Various image manipulations
* transform: Geometrical transformations
* insar:     Interferometric processing
* polar:     Polarimetric processing

More information about modules and contents (replace 'module' by correct name, e.g. 'filter):
    
    >>> module.info()

## Example batch usage

    ./ pyrat.py -b
    >>> x1 = load.rat(filename='abc.rat')
    >>> x2 = filter.lee(looks=3)
    >>> x3 = filter.boxcar(layer=x1)
    >>> save.pixmap(filename='abc.jpg', layer=x2)
    >>> var = getdata(layer=x1)
    >>> show()

## Implementing your own modules

PyRat has a quite simple programming interface. Have a look at the file 'pyrat/filter/Template.py',
this should explain at least the basics of programming own modules. Put your own code
in the 'plugins' directory, it is automatically scanned at startup. PyRat will automatically
attach your code to the GUI and run it using parallel processing.

For more detailed questions, please contact the authors directly.
