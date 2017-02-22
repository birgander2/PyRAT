import pyrat
import h5py


class HDF5(pyrat.Worker):
    """
    Generic HDF5 writer (experimental)
    It should work, but might not be compatible with some other code we use. Probably needs to be improved!!!

    :author: Andreas Reigber
    """
    gui = {'menu': 'File|Export to', 'entry': 'HDF5'}
    para = [{'var': 'file', 'value': '', 'type': 'savefile', 'text': 'Save as HDF5', 'extensions': 'HDF5 (*.hd5)'}]

    def __init__(self, *args, **kwargs):
        super(HDF5, self).__init__(*args, **kwargs)
        self.name = "HDF5 EXPORT"
        self.nthreads = 1
        if len(args) == 1:
            self.file = args[0]

    def run(self, *args, **kwargs):
        if isinstance(self.file, tuple):                                       # remove file type if present
            self.file = self.file[0]

        if isinstance(self.layer, list):
            layers = self.layer
        else:
            layers = [self.layer]

        self.file = h5py.File(self.file, 'w')
        for k, layer in enumerate(layers):
            query = pyrat.data.queryLayer(layer)
            self.dset = self.file.create_dataset("D" + str(k + 1), query['shape'], dtype=query['dtype'])
            self.layer_extract(self.block_writer, silent=False, layer=layer, **kwargs)
            meta = pyrat.data.getAnnotation(layer=layer)
            for key, val in meta.items():
                self.dset.attrs[key] = val
        self.file.close()
        return layers


    def block_writer(self, array, **kwargs):
        self.dset[..., kwargs['block'][0] + kwargs['valid'][0]:kwargs['block'][0] + kwargs['valid'][1], :] \
            = array[..., kwargs['valid'][0]:kwargs['valid'][1], :]
        return True


@pyrat.docstringfrom(HDF5)
def hdf5(*args, **kwargs):
    return HDF5(*args, **kwargs).run(*args, **kwargs)
