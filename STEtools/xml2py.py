import numpy as np
import xml.etree.ElementTree as ET

class Xml2Py:
    """Turns a STEP XML document (e.g. processing parameters) into a structure.

    Example usage:
    pp = Xml2Py(<path to pp file>)
    pp.v0
    (Python prints 89.1)
    pp.v0?
    (Type is numpy double)
    pp.r?
    (Numpy ndarray of doubles)

    Currently does not support complex data, pointer arrays and structure
    arrays.

     :author: Marc Jaeger

     :param fileName: Path to XML document.
     :type fileName: string
     
     :param elem: For internal use only! Do not use.
     :type elem: XML element object
     
    """

    def __init__(self, fileName, elem=None):
        if (elem == None):
            tree = ET.parse(fileName)
            root = tree.getroot()
            elem = root[0]

        self.xml2struct(elem)


    def xml2struct(self, elem):
        for f in elem:
            self.xml2field(f)


    def xml2field(self, elem, name=None):
        typElem = elem.find('datatype')
        dType = typElem.text
        dDim = typElem.attrib['length']
        dDim = np.asarray([np.long(d) for d in dDim.split()])[::-1]
        dLen = np.prod(dDim)

        if (name == None):
            name = elem.attrib['name']

        valElem = elem.find('value')
        if (dType == 'pointer'):
            self.xml2field(valElem.find('parameter'),elem.attrib['name'])
            return

        if (dType == 'struct'):
            o = Xml2Py(None, valElem[0])
            setattr(self,name,o)
            return

        val = elem.find('value').text    
        if (dLen > 1):
            val = val.strip('[]').split(',')

        conv = {'int': np.int, 'long': np.long, 'float': np.float, 'double': np.double, 'string': lambda s:s}
        try:
            if (dLen > 1):
                val = np.asarray([conv[dType](v) for v in val]).reshape(dDim)
            else:
                val = conv[dType](val)
        except KeyError:
            print 'WARNING: Unsupported data type {} in field {}! Ignoring...'.format(dType, name)
            return

        setattr(self,name,val)
