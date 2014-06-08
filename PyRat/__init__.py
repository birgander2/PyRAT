# PyRat __init__

from LayerData import *
from FilterWorker import *
from ImportWorker import *
from ExportWorker import *
from LayerWorker import *
import Layer
import Filter
import Import
import Export
import Viewer
import os, logging
import multiprocessing
import sys

#logging.basicConfig(format='  %(asctime)s %(levelname)s: %(message)s', level=logging.DEBUG)
logging.basicConfig(format='  %(levelname)s: %(message)s', level=logging.INFO)

def init(filename = 'PyRat.hd5'):
    global Data, MP_Pool
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    if os.path.exists(filename):
        os.remove(filename)
    Data = LayerData(filename)
    MP_Pool = multiprocessing.Pool(processes=8)


def exit(filename = 'PyRat.hd5'):
    global MP_Pool
    if os.path.exists(filename):
        os.remove(filename)
    MP_Pool.close()

def cleanup():
    global Data
    Data.close()
    del Data
    os.system("h5repack PyRAT.hd5 PyRAT.hd5.temp")
    os.system("mv PyRAT.hd5.temp PyRAT.hd5")
    Data = LayerData('PyRAT.hd5')

