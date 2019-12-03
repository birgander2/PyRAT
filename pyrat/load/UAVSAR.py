"""
Import UAVSAR data (SLC's and Kz's) to perform TomoSAR

As in : https://github.com/mdenbina/kapok

Author: Gustavo Daniel Martín del Campo Becerra
E-Mail: Gustavo.MartindelCampoBecerra@dlr.de
"""

import logging
import os
import os.path
import time
import pyrat
import numpy as np

from glob import glob
from scipy.ndimage.interpolation import zoom
from pyrat.tools import bcolors



def findsegment(az, azsize):
    """For a given azimuth index, return the segment number and
    the azimuth index within the given segment.

    Arguments:
        az: Azimuth index of interest.
        azsize: List containing the azimuth size of each segment.

    Returns:
        seg: Segment number.
        azoff: Azimuth index within the segment.

    """
    azstart = np.insert(np.cumsum(azsize), 0, 0)[0:-1]

    if az <= np.sum(azsize):
        seg = np.max(np.where(azstart <= az))
        azoff = az - azstart[seg]
    else:
        seg = azstart.shape[0] - 1
        azoff = azsize[seg]
        logging.info('')
        logging.info(bcolors.WARNING + 'SLC row index of ' + str(az) + ' is larger than the size of the data. Returning maximum index.' + bcolors.ENDC)

    return seg, azoff



def mlook(data, mlwin):
    """Multilook/rebin image to smaller image by averaging.

    Arguments:
        data: Array (up to 4D) containing data to multilook.
        mlwin (tuple, int): Tuple of ints containing the smoothing window
            sizes in each dimension.

    Returns:
        mldata: Array containing multilooked data.

    """
    if mlwin == (1, 1):
        return data

    data = np.asarray(data)
    n_dim = len(data.shape)

    nshape = np.array(data.shape) // np.array(list(mlwin) + [1] * (n_dim - len(mlwin)))

    sh = np.array([[nshape[i], data.shape[i] // nshape[i]] for i, x in enumerate(nshape)]).flatten()

    if n_dim == 2:
        if not any(np.mod(data.shape, nshape)):
            return data.reshape(sh).mean(-1).mean(1)
        else:
            return data[0:sh[0] * sh[1], 0:sh[2] * sh[3]].reshape(sh).mean(-1).mean(1)
    elif n_dim == 3:
        if not any(np.mod(data.shape, nshape)):
            return data.reshape(sh).mean(1).mean(2).mean(3)
        else:
            return data[0:sh[0] * sh[1], 0:sh[2] * sh[3], 0:sh[4] * sh[5]] \
                .reshape(sh).mean(1).mean(2).mean(3)
    elif n_dim == 4:
        if not any(np.mod(data.shape, nshape)):
            return data.reshape(sh).mean(1).mean(2).mean(3).mean(4)
        else:
            return data[0:sh[0] * sh[1], 0:sh[2] * sh[3], 0:sh[4] * sh[5], 0:sh[6] * sh[7]] \
                .reshape(sh).mean(1).mean(2).mean(3).mean(4)
    else:
        logging.info('')
        logging.info(bcolors.FAIL + 'Given number of dimensions not considered. Aborting.' + bcolors.ENDC)

    return



def getslcblock(file, rngsize, azstart, azend, rngbounds=None, file2=None,
                azsize=None):
    """Load SLC data into a NumPy array buffer.  If the file2 argument is
    specified in the arguments, this function will treat the two files as
    consecutive segments, and will join them.

    Arguments:
        file (str): Filename of the first SLC.
        rngsize (int): Number of columns (range bins).  Same for both SLCs.
        azstart (int): Azimuth index at which to start the buffer, in the
            first SLC.
        azend (int): Azimuth index at which to end the buffer.  If file2
            is specified, this is an azimuth index in the second SLC.  If
            file2 is not specified, this is an azimuth index in the first SLC.
            The row specified by azend is not actually included in the buffer,
            as in the Python range() function.  (azend-1) is the last line
            included in the buffer.  To load the entire SLC, azend should be
            equal to the number of rows in SLC.
        rngbounds (tuple, int): Starting and ending range bounds, if range
            subsetting is desired.
        file2 (str): Filename of the second SLC, if combining multiple
            segments is desired.
        azsize (int): Number of rows (azimuth bins) for the first SLC.  Only
            required if file2 is specified.  Otherwise we only load in the
            lines of the SLC before azend.

    Returns:
        block: NumPy array of complex64 datatype, containing the loaded SLC
        data between the specified azimuth bounds.

    """
    if file2 is None:
        byteoffset = rngsize * 8 * azstart
        slc = np.memmap(file, dtype='complex64', mode='c', offset=byteoffset, shape=(azend - azstart, rngsize))
        if rngbounds is not None:
            slc = slc[:, rngbounds[0]:rngbounds[1]]
        return slc
    elif azsize is None:
        logging.info('')
        logging.info(bcolors.FAIL + '"file2" argument specified, but "azsize" argument missing. Aborting.' + bcolors.ENDC)
        return
    else:
        byteoffset = rngsize * 8 * azstart
        slca = np.memmap(file, dtype='complex64', mode='c', offset=byteoffset, shape=(azsize - azstart, rngsize))
        slcb = np.memmap(file2, dtype='complex64', mode='c', shape=(azend, rngsize))
        if rngbounds is not None:
            slca = slca[:, rngbounds[0]:rngbounds[1]]
            slcb = slcb[:, rngbounds[0]:rngbounds[1]]
        return np.vstack((slca, slcb))



class Ann(object):
    """Class for loading and interacting with a UAVSAR annotation file."""

    def __init__(self, file):
        """Load in the specified .ann file as a list of strings and
            initialize.

        Arguments:
            file (str): UAVSAR annotation filename to load.

        """
        self.file = file
        fd = open(self.file, 'r')
        self.ann = fd.read().split('\n')
        return

    def query(self, keyword):
        """Query the annotation file for the specified annotation keyword.

        Arguments:
            keyword (str): The keyword to query.

        Returns:
            value: The value of the specified keyword.

        """
        for n in range(len(self.ann)):
            if self.ann[n].startswith(keyword):
                try:
                    val = self.ann[n].rsplit('=')[-1].split(';')[0].split()[0]
                    val = np.array(val, dtype='float')  # if we can convert the string to a number, do so
                    if (val - np.floor(val)) == 0:
                        val = np.array(val,
                                       dtype='int')  # if it's an integer, convert it to one (e.g., number of samples)
                    return val
                except ValueError:  # if we can't convert the string to a number, leave it as a string
                    val = self.ann[n].split('=', maxsplit=1)[-1].split(';')[0].strip()
                    return val

        return None



class Scene(pyrat.Worker):
    """Scene object for reading the UAVSAR dataset composed of only one segment.
        Arguments:
            infile (str): Path and filename of a UAVSAR .ann file containing the
                metadata for a UAVSAR SLC stack.
            azbounds: List containing starting and ending SLC azimuth index for
                desired subset.  Data outside these bounds will not be loaded.
                Default is to import all data.
            rngbounds: List containing starting and ending SLC range index for
                desired subset.  Data outside these bounds will not be loaded.
                Default is to import all data.
            select: List containing desired track indices to import.  For example,
                tracks=[0,1] will import only the first two tracks listed in the
                annotation file.  Default: All tracks imported.
            master: Master track index.
    """
    para = [
        {'var': 'infile', 'value': ''},
        {'var': 'azbounds', 'value': None},
        {'var': 'rngbounds', 'value': None},
        {'var': 'select', 'value': None},
        {'var': 'master', 'value': 0}
    ]

    def __init__(self, *args, **kwargs):
        """Scene initialization method."""
        super(Scene, self).__init__(*args, **kwargs)
        self.name = "UAVSAR COV AND KZ IMPORT"
        self.blockprocess = False
        return



    ###################
    # Getters
    def get_npols(self):
        return self.__num_pol

    def get_ntracks(self):
        return self.__num_tracks

    def get_site(self):
        return self.__site

    def get_url(self):
        return self.__url

    def get_sensor(self):
        return self.__sensor

    def get_stack_name(self):
        return self.__stack_name



    ###################
    # Loading data
    def run(self, *args, **kwargs):

        # Load the annotation file.
        try:
            ann = Ann(self.infile)
        except:
            logging.info('')
            logging.info(bcolors.FAIL + "Cannot load UAVSAR annotation file. Aborting." + bcolors.ENDC)
            return

        # Get SLC dimensions
        rngsize_slc = int(ann.query('slc_1_1x1 Columns'))
        azsize_slc = int(ann.query('slc_1_1x1 Rows'))

        # Get track filenames and number of tracks.
        temp = ann.query('stackline1')
        num = 1
        if temp is not None:
            tracknames = [temp.split('_L090')[0]]
            num += 1
            temp = ann.query('stackline' + str(num))
            while temp is not None:
                tracknames.append(temp.split('_L090')[0])
                num += 1
                temp = ann.query('stackline' + str(num))
        else:
            logging.info('')
            logging.info(bcolors.FAIL + "Cannot find track names in UAVSAR annotation file. Aborting." + bcolors.ENDC)
            return

        # Subset track names if desired tracks were specified:
        tracknames = np.array(tracknames)
        if self.select is not None:
            tracks = np.array(self.select, dtype='int')
            tracknames = tracknames[tracks]

        logging.info('')
        logging.info(bcolors.OKBLUE + "Importing metadata" + bcolors.ENDC)

        self.__tracks = np.array(tracknames, dtype='S')
        self.__num_tracks = len(tracknames)  # Number of Tracks
        self.__num_baselines  = int(self.__num_tracks * (self.__num_tracks - 1) / 2)  # Number of Baselines
        self.__num_pol = int(3)  # Number of Polarizations (HH, sqrt(2)*HV, VV)
        self.__num_elements = self.__num_tracks * self.__num_pol

        self.__stack_name = ann.query('Stack Name')
        self.__site = ann.query('Site Description')
        self.__url = ann.query('URL')
        self.__sensor = 'UAVSAR'

        logging.info('')
        logging.info(bcolors.OKBLUE + 'Stack ID: ' + self.__stack_name + bcolors.ENDC)
        logging.info(bcolors.OKBLUE + 'Site Description: ' + self.__site + bcolors.ENDC)
        logging.info(bcolors.OKBLUE + 'URL: ' + self.__url + bcolors.ENDC)

        self.__average_altitude = float(ann.query('Average Altitude'))
        self.__image_starting_slant_range = float(ann.query('Image Starting Slant Range')) * 1000
        self.__slc_azimuth_pixel_spacing = float(ann.query('1x1 SLC Azimuth Pixel Spacing'))
        self.__slc_slant_range_pixel_spacing = float(ann.query('1x1 SLC Range Pixel Spacing'))

        # Get SCH Peg from Annotation File
        self.__peglat = float(ann.query('Peg Latitude'))
        self.__peglon = float(ann.query('Peg Longitude'))
        self.__peghdg = float(ann.query('Peg Heading'))

        if (self.azbounds is not None) or (self.rngbounds is not None):
            self.__subset = True
        else:
            self.__subset = False

        # Check azimuth bounds for validity.
        if self.azbounds is None:
            self.azbounds = [0, np.sum(azsize_slc)]

        if self.azbounds[1] <= self.azbounds[0]:
            logging.info('')
            logging.info(bcolors.FAIL + 'Invalid azimuth bounds. Must be ascending. Aborting.' + bcolors.ENDC)
            return

        if self.azbounds[0] < 0:
            logging.info('')
            logging.info(bcolors.WARNING + 'Lower azimuth bound (' + str(self.azbounds[0]) + ') is less than zero. Setting lower azimuth bound to zero.' + bcolors.ENDC)
            self.azbounds[0] = 0

        if self.azbounds[1] > np.sum(azsize_slc):
            logging.info('')
            logging.info(bcolors.WARNING + 'Upper azimuth bound (' + str(self.azbounds[1]) + ') greater than number of SLC lines (' + str(azsize_slc) + ').' + bcolors.ENDC)
            logging.info(bcolors.WARNING + 'Setting upper azimuth bound to ' + str(azsize_slc) + '.' + bcolors.ENDC)
            self.azbounds[1] = np.sum(azsize_slc)

        # Check range bounds for validity.
        if self.rngbounds is None:
            self.rngbounds = [0, rngsize_slc]

        if self.rngbounds[1] <= self.rngbounds[0]:
            logging.info('')
            logging.info(bcolors.FAIL + 'Invalid range bounds. Must be ascending. Aborting.' + bcolors.ENDC)
            return

        if self.rngbounds[0] < 0:
            logging.info('')
            logging.info(bcolors.WARNING + 'Lower range bound (' + str(self.rngbounds[0]) + ') is less than zero.  Setting lower azimuth bound to zero.' + bcolors.ENDC)
            self.rngbounds[0] = 0

        if self.rngbounds[1] > rngsize_slc:
            logging.info('')
            logging.info(bcolors.WARNING + 'Upper range bound (' + str(self.rngbounds[1]) + ') greater than number of SLC columns (' + str(rngsize_slc) + ').' + bcolors.ENDC)
            logging.info(bcolors.WARNING + 'Setting upper range bound to ' + str(rngsize_slc) + '.' + bcolors.ENDC)
            self.rngbounds[1] = rngsize_slc

        # Image dimensions:
        azsize = (self.azbounds[1] - self.azbounds[0])
        rngsize = (self.rngbounds[1] - self.rngbounds[0])
        self.__dim = (azsize, rngsize)
        self.__dim_slc = (int(np.sum(azsize_slc)), int(rngsize_slc))
        self.__dim_slc_stack = (self.__num_pol, azsize, rngsize)
        self.__dim_kz_stack = (self.__num_tracks, azsize, rngsize)

        # Path containing SLCs and other files (assume in same folder as .ann):
        datapath = os.path.dirname(self.infile)
        if datapath != '':
            datapath = datapath + '/'

        # Get filenames of SLCs for each track, in polarization order HH, HV, VV.
        slcfiles = []
        for tr in range(self.__num_tracks):
            for pol in ['HH', 'HV', 'VV']:
                file = glob(datapath + tracknames[tr] + '*' + pol + '_*_s' + str(1) + '_1x1.slc')

                if len(file) == 1:
                    slcfiles.append(file[0])
                elif len(file) > 1:
                    logging.info('')
                    logging.info(bcolors.FAIL + 'Too many SLC files matching pattern: "' + datapath + tracknames[tr] + '*_' + pol + '_*_1x1.slc' + '". Aborting.' + bcolors.ENDC)
                    return
                else:
                    logging.info('')
                    logging.info(bcolors.FAIL + 'Cannot find SLC file matching pattern: "' + datapath + tracknames[tr] + '*_' + pol + '_*_1x1.slc' + '".  Aborting.' + bcolors.ENDC)
                    return

        slcstack = np.empty(self.__dim_slc_stack, dtype='complex64')
        slc_layer = []
        index = 0
        for slcnum in range(0, self.__num_elements):
            logging.info('')
            logging.info(bcolors.OKBLUE + 'Loading SLCs: ' + str(slcnum + 1) + '/' + str(self.__num_elements) + '. (' + time.ctime() + ')' + bcolors.ENDC)

            file = slcfiles[slcnum]
            slc = np.memmap(file, dtype='complex64', mode='c', shape=(self.__dim_slc))

            if self.rngbounds is not None:
                slc = slc[:, self.rngbounds[0]:self.rngbounds[1]]
            if self.azbounds is not None:
                slc = slc[self.azbounds[0]:self.azbounds[1], :]

            if (slcnum % 3) == 1:  # HV Polarization
                slcstack[index,...] = np.sqrt(2) * slc
            else:  # HH or VV Polarization
                slcstack[index,...] = slc

            index += 1
            if index == 3:
                slc_layer.append(pyrat.adddata(slcstack))
                index = 0

        # LLH Subset Bounds and Offsets for Trimming Multilooked Arrays
        mlwin_lkv = (int(ann.query('Number of Azimuth Looks in 2x8 SLC')), int(ann.query('Number of Range Looks in 2x8 SLC')))
        azllhstart = self.azbounds[0] // mlwin_lkv[0]
        azllhend = self.azbounds[1] // mlwin_lkv[0]
        azllhoffset = self.azbounds[0] % mlwin_lkv[0]
        rngllhstart = self.rngbounds[0] // mlwin_lkv[1]
        rngllhend = self.rngbounds[1] // mlwin_lkv[1]
        rngllhoffset = self.rngbounds[0] % mlwin_lkv[1]

        # Import .kz files

        self.__kz_units = 'radians/meter'
        self.__kz_description = 'Interferometric Vertical Wavenumber'
        self.__kz_indexing = 'track'

        kz_stack = np.empty(self.__dim_kz_stack, dtype='complex64')
        for tr in range(0, self.__num_tracks):
            logging.info('')
            logging.info(bcolors.OKBLUE + 'Importing kz between track ' + str(tr) + ' and reference track. (' + time.ctime() + ')' + bcolors.ENDC)

            kz_temp = None
            file = glob(datapath + tracknames[tr] + '*s' + str(1) + '*.kz')
            file_alt = glob(datapath + tracknames[tr] + '*.kz')

            if len(file) >= 1:
                kz_rows = int(ann.query('lkv_' + str(1) + '_2x8 Rows'))
                kz_cols = int(ann.query('lkv_' + str(1) + '_2x8 Columns'))
                if kz_temp is None:
                    kz_temp = np.memmap(file[0], dtype='float32', mode='r', shape=(kz_rows, kz_cols))
                else:
                    kz_temp = np.vstack(
                        (kz_temp, np.memmap(file[0], dtype='float32', mode='r', shape=(kz_rows, kz_cols))))
            elif (len(file_alt) >= 1):
                kz_rows = int(ann.query('lkv_' + str(1) + '_2x8 Rows'))
                kz_cols = int(ann.query('lkv_' + str(1) + '_2x8 Columns'))
                if kz_temp is None:
                    kz_temp = np.memmap(file_alt[0], dtype='float32', mode='r', shape=(kz_rows, kz_cols))
                else:
                    kz_temp = np.vstack((kz_temp, np.memmap(file_alt[0], dtype='float32', mode='r', shape=(kz_rows, kz_cols))))
            else:
                break

            # kz files appear to be kz_i0 rather than kz_0i.  Multiplying by -1 so that in the Scene object, kz_ij = kz[j] - kz[i].
            kz_stack[tr,...] = -1*zoom(kz_temp[azllhstart:azllhend,rngllhstart:rngllhend],mlwin_lkv)[azllhoffset:(azllhoffset+azsize),rngllhoffset:(rngllhoffset+rngsize)]

        kz = np.empty(self.__dim_kz_stack, dtype=np.complex64)
        for i in range(self.__num_tracks):
            kz[i, ...] = (kz_stack[self.master,...] - kz_stack[i,...])

        kz_layer = pyrat.adddata(kz)

        logging.info('')
        logging.info(bcolors.OKBLUE + 'Complete. (' + time.ctime() + ')' + bcolors.ENDC)

        return kz_layer, slc_layer


    def quicklook(self, tr=0, pol='hh', mlwin=(40, 10), savefile=None, crop=None):
        """Display a quick look intensity image for a given UAVSAR SLC stack.

        Arguments:
            infile (str): Input annotation file of the UAVSAR stack.
            tr (int): Track index of the desired image.  Default: 0
                (first track in .ann file).
            pol (str): Polarization str of the desired image.  Options are 'hh',
                'hv', or 'vv'.  Default: 'hh'
            mlwin (tuple): Multi-looking window size to use for the quick look
                image.  Note:  The original SLC azimuth indices will be displayed
                on the axes of the image, so that the image can be used as a guide
                for suitable values of the azbounds and rngbounds keywords in
                uavsar.load.  These multi-looking windows are in terms
                of the original 1x1 SLC image size, not the 8x2 or 4x1
                image sizes.
            savefile (str): Output path and filename to save the displayed image.

        """
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches

        # Load the annotation file.
        try:
            ann = Ann(self.infile)
        except:
            logging.info('')
            logging.info(bcolors.FAIL + "Cannot load UAVSAR annotation file. Aborting." + bcolors.ENDC)
            return

        # Get track filenames and number of tracks.
        temp = ann.query('stackline1')
        num = 1
        if temp is not None:
            tracknames = [temp.split('_L090')[0]]
            num += 1
            temp = ann.query('stackline' + str(num))
            while temp is not None:
                tracknames.append(temp.split('_L090')[0])
                num += 1
                temp = ann.query('stackline' + str(num))
        else:
            logging.info('')
            logging.info(
                bcolors.FAIL + "Cannot find track names in UAVSAR annotation file. Aborting." + bcolors.ENDC)
            return

        # Path containing SLCs and other files (assume in same folder as .ann):
        datapath = os.path.dirname(self.infile)
        if datapath != '':
            datapath = datapath + '/'

        pol = pol.upper()

        slcfiles = []
        file = glob(datapath + tracknames[tr] + '*' + pol + '_*_s' + str(1) + '_2x8.slc')

        if len(file) == 1:
            slcfiles.append(file[0])
            slcwindow = (8, 2)
            rngsize_slc = ann.query('slc_1_2x8 Columns')
            azsize_slc = ann.query('slc_1_8x2 Rows')
        elif len(file) > 1:
            logging.info('')
            logging.info(
                bcolors.FAIL + 'Too many SLC files matching pattern: "' + datapath + tracknames[tr] + '*_' + pol + '_*_2x8.slc' + '".  Aborting.' + bcolors.ENDC)
            return
        else:
            file = glob(datapath + tracknames[tr] + '*' + pol + '_*_s' + str(1) + '_4x1.slc')

            if len(file) == 1:
                slcfiles.append(file[0])
                slcwindow = (4, 1)
                rngsize_slc = ann.query('slc_1_1x4 Columns')
                azsize_slc = ann.query('slc_1_1x4 Rows')
            elif len(file) > 1:
                logging.info('')
                logging.info(bcolors.FAIL + 'Too many SLC files matching pattern: "' + datapath + tracknames[
                    tr] + '*_' + pol + '_*_1x4.slc' + '". Aborting.' + bcolors.ENDC)
                return
            else:
                file = glob(datapath + tracknames[tr] + '*' + pol + '_*_s' + str(1) + '_1x1.slc')

                if len(file) == 1:
                    slcfiles.append(file[0])
                    slcwindow = (1, 1)
                    rngsize_slc = ann.query('slc_1_1x1 Columns')
                    azsize_slc = ann.query('slc_1_1x1 Rows')
                elif len(file) > 1:
                    logging.info('')
                    logging.info(bcolors.FAIL + 'Too many SLC files matching pattern: "' + datapath + tracknames[
                        tr] + '*_' + pol + '_*_1x1.slc' + '". Aborting.' + bcolors.ENDC)
                    return
                else:
                    logging.info('')
                    logging.info(bcolors.FAIL + 'Cannot find SLC file matching pattern: "' + datapath + tracknames[
                            tr] + '*_' + pol + '_*_1x1.slc' + '". Aborting.' + bcolors.ENDC)
                    return

        mlwin = (mlwin[0] // slcwindow[0], mlwin[1] // slcwindow[1])
        azsize = np.sum(azsize_slc) // mlwin[0]
        rngsize = rngsize_slc // mlwin[1]

        # Get SLC dimensions for the 1x1 SLCs -- these get displayed on the plot axes.
        rngsize_slc1x1 = ann.query('slc_1_1x1 Columns')
        azsize_slc1x1 = ann.query('slc_1_1x1 Rows')

        # Load and multilook quicklook intensity image.
        qlimage = np.zeros((azsize, rngsize), dtype='float32')

        num_blocks = int(np.sum(azsize_slc) * rngsize_slc * 8 / 1e9)
        if num_blocks < 2:
            num_blocks = 2
        az_vector = np.round(np.linspace(0, azsize, num=num_blocks + 1)).astype('int')
        az_vector[num_blocks] = azsize

        for n, azstart in enumerate(az_vector[0:-1]):
            azend = az_vector[n + 1]
            azstart_slc = azstart * mlwin[0]
            azend_slc = azend * mlwin[0]
            seg_start, azoffset_start = findsegment(azstart_slc, azsize_slc)
            seg_end, azoffset_end = findsegment(azend_slc, azsize_slc)

            file = slcfiles[seg_start]

            if seg_start == seg_end:
                slc = getslcblock(file, rngsize_slc, azoffset_start, azoffset_end)
            else:
                file2 = slcfiles[seg_end]
                slc = getslcblock(file, rngsize_slc, azoffset_start, azoffset_end, file2=file2,
                                  azsize=azsize_slc[seg_start])

            qlimage[azstart:azend, :] = mlook(np.real(slc * np.conj(slc)), mlwin)

        qlimage = np.real(qlimage)
        qlimage[qlimage <= 1e-10] = 1e-10
        qlimage = 10 * np.log10(qlimage)

        fig1, ax1 = plt.subplots(1)

        cs1 = ax1.imshow(qlimage, vmin=-25, vmax=0, cmap='gray', aspect=0.25, interpolation='nearest',
                         extent=(0, rngsize_slc1x1, np.sum(azsize_slc1x1), 0))

        if crop is not None:
            # Create a Rectangle patch
            rect1 = patches.Rectangle((crop[2], crop[0]), crop[3]- crop[2], crop[1]- crop[0],
                                      linewidth=2, edgecolor='r', facecolor='none')
            # Add the patch to the Axes
            ax1.add_patch(rect1)

        cbar = fig1.colorbar(cs1)
        cbar.set_label(label=pol + ' Backscatter (dB)', fontsize=20)
        cbar.ax.tick_params(labelsize=18)

        ax1.set_xlabel('Range Index', fontsize=20)
        ax1.set_ylabel('Azimuth Index', fontsize=20)
        ax1.tick_params(labelsize=18)
        fig1.tight_layout()

        if savefile is not None:
            fig1.savefig(savefile, dpi=125, bbox_inches='tight', pad_inches=0.1)

        return


@pyrat.docstringfrom(Scene)
def scene(*args, **kwargs):
    return Scene(*args, **kwargs).run(*args, **kwargs)

@pyrat.docstringfrom(Scene)
def sceneObject(*args, **kwargs):
    return Scene(*args, **kwargs)