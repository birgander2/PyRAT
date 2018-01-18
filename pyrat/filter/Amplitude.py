import pyrat
import numpy as np


class FlattenAmp(pyrat.Worker):
    """
    Removes amplitude trends

    :author: Andreas Reigber
    """
    gui = {'menu': 'SAR|Amplitude', 'entry': 'Remove trends'}
    para = [
        {'var': 'order', 'value': 1, 'type': 'int', 'range': [1, 9], 'text': 'polynomial order'},
        {'var': 'axis', 'value': 'range', 'type': 'list', 'range': ['range', 'azimuth', 'both'], 'text': 'adjust'}
    ]

    def run(self, *args, **kwargs):
        layer = pyrat.data.active

        # STEP1: Estimate profiles
        azprof, rgprof = self.layer_accumulate(self.estimate_profiles, combine=self.combine_profiles)

        # STEP2: Fit correction
        rgprof /= np.mean(rgprof, axis=-1, keepdims=True)
        azprof /= np.mean(azprof, axis=-1, keepdims=True)

        # todo: from here on adapt to nd-data sets
        rgaxis = np.arange(rgprof.shape[-1])
        azaxis = np.arange(azprof.shape[-1])
        rgcorr = np.empty_like(rgprof)
        azcorr = np.empty_like(azprof)
        if rgprof.ndim == 1:
            rgcorr = np.polyval(np.polyfit(rgaxis, rgprof, self.order), rgaxis)
            azcorr = np.polyval(np.polyfit(azaxis, azprof, self.order), azaxis)
        elif rgprof.ndim == 2:
            for k in range(rgprof.shape[0]):
                rgcorr[k, :] = np.polyval(np.polyfit(rgaxis, rgprof[k, :], self.order), rgaxis)
                azcorr[k, :] = np.polyval(np.polyfit(azaxis, azprof[k, :], self.order), azaxis)
        elif rgprof.ndim == 3:
            for k in range(rgprof.shape[0]):
                for l in range(rgprof.shape[1]):
                    rgcorr[k, l, :] = np.polyval(np.polyfit(rgaxis, rgprof[k, l, :], self.order), rgaxis)
                    azcorr[k, l, :] = np.polyval(np.polyfit(azaxis, azprof[k, l, :], self.order), azaxis)

        # STEP3: Apply correction
        outlayer = self.layer_process(self.applyfix, axis=self.axis, correction=(azcorr, rgcorr), siltent=False,
                                      **kwargs)
        pyrat.activate(outlayer)
        return outlayer

    @staticmethod
    def applyfix(data, axis='range', correction=(1.0, 1.0), block=(0, 0, 0, 0), **kwargs):
        rgcorr = correction[1]
        azcorr = correction[0][..., block[0]:block[1]]
        if axis == 'range' or axis == 'both':
            data /= rgcorr[..., np.newaxis, :]
        if axis == 'azimuth' or axis == 'both':
            data /= azcorr[..., np.newaxis]
        return data

    @staticmethod
    def estimate_profiles(data, **kwargs):
        valid = kwargs['valid']
        rgprof = np.mean(np.abs(data), axis=-2)
        azprof = np.mean(np.abs(data), axis=-1)[..., valid[0]:valid[1]]
        return azprof, rgprof

    @staticmethod
    def combine_profiles(values, **kwargs):
        azshp = list(values[0][0].shape)
        azshp[-1] = 0
        rgprof = np.zeros(values[0][1].shape)
        azprof = np.zeros(azshp)
        for val in values:
            rgprof += val[1]
            azprof = np.append(azprof, val[0], axis=-1)
        return azprof, rgprof


@pyrat.docstringfrom(FlattenAmp)
def flattenamp(*args, **kwargs):
    return FlattenAmp(*args, **kwargs).run(*args, **kwargs)
