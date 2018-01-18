import numpy as np

def sarscale(img, factor=2.5):
    """
    Returns a SAR-bytscaled version of an array for visualisation, clipped
    between zero and 2.5 * mean(array)
    """
    scl = np.nanmean(img)
    if scl == 0:
        return np.zeros_like(img, dtype=np.uint8)
    else:
        return np.uint8(np.clip(255.0 * img / factor / scl, 0, 255))


def phascale(img):
    """
    Returns a byte version of an phase array (scaled between -!pi to !pi)
    """
    return np.uint8(np.clip(img / np.pi * 128 + 127, 0.0, 255.0))


def cohscale(img):
    """
    Returns a byte version of an coherence array (scaled between 0.0 and 1.0)
    """
    return np.uint8(np.clip(img * 255.0, 0.0, 255.0))


def histscale(img, factor=0.05):
    """
    Returns a histgram scaled version of array. Cuts the percentage given in factor.
    """
    hist = np.histogram(img, bins=1000, range=[0.0, 10 * np.mean(img)])
    hist[0][0] = 0
    mh = np.sum(hist[0]) * factor
    st = 0
    while np.sum(hist[0][:st]) < mh:
        st += 1
    en = len(hist[0])
    while np.sum(hist[0][en:]) < mh:
        en -= 1
    return np.uint8(np.clip((img - hist[1][st]) / (hist[1][en] - hist[1][st]) * 255, 0, 255))


def npabs(array, **kwargs):
    """
    Helper routines to simulate behaviour of cython multithreaded version. Identical to np.abs()
    """
    return np.abs(array)


def subsample(args):
    """
    Rebin / Congrid variant
    """
    arr, shape, mode = args  # unpack arguments

    if mode == 'phase' and np.iscomplexobj(arr):
        arr = np.angle(arr)
    if np.iscomplexobj(arr):
        arr = np.abs(arr)

    if arr.shape == shape:
        return arr

    oshap = arr.shape

    for d in range(arr.ndim):
        n1 = shape[d]
        n2 = oshap[d]
        if n1 < n2:
            s = list(arr.shape)
            s.insert(d + 1, n2 // n1)
            s[d] = n1
            if mode == 'phase':
                arr = np.angle(np.exp(1j * arr.reshape(s)).mean(d + 1))
            elif mode == 'lables':
                arr = np.take(arr.reshape(s), 0, d + 1)
            else:
                arr = arr.reshape(s).mean(d + 1)
        else:
            arr = arr.repeat(n1 // n2, axis=d)
    return arr
