"""
Provides :func:`loadmatbunch` for more convenient access to the
data in Matlab .mat files.

Also *datenum_to_day* and *datenum_to_mpl* for
converting MATLAB datenums to decimal days and matplotlib
datenums, respectively.
"""
import datetime

import numpy as np
from scipy.io import loadmat

from misc import Bunch

def _unicode(arr):
    """
    loadmat seems to be mishandling strings when there is a difference
    in byte order between the machine that wrote the file and the one
    that is reading the file.  The result is, e.g.,
    
    """
    try:
        return str(arr)
    except UnicodeEncodeError:
        dt = arr.dtype.newbyteorder('S')
        return str(arr.view(dt))

def crunch(arr, masked=True):
    arr = arr.squeeze()
    if arr.ndim == 0:
        kind = arr.dtype.kind
        if kind == 'f':
            return float(arr)
        if kind in 'ui':
            return int(arr)
        if kind == 'U':
            try:
                return _unicode(arr)
            except UnicodeDecodeError:
                return "Could not decode."
        if kind == 'S':
            return str(arr)
        if kind == 'O':
            return arr
        return arr  # warn?  Other kinds need to be handled?
    if masked and arr.dtype.kind == 'f':  # check for complex also
        arrm = np.ma.masked_invalid(arr)
        if arrm.count() < arrm.size:
            arr = arrm
        else:
            arr = np.array(arr) # copy to force a read
    else:
        arr = np.array(arr)
    return arr

def structured_to_bunch(arr, masked=True):
    if arr.dtype.kind == 'V' and arr.shape == (1,1):
        b = Bunch()
        x = arr[0,0]
        for name in x.dtype.names:
            b[name] = structured_to_bunch(x[name], masked=masked)
        return b
    return crunch(arr, masked=masked)


def loadmatbunch(fname, masked=True):
    """
    Wrapper for loadmat that dereferences (1,1) object arrays,
    converts floating point arrays to masked arrays, and uses
    nested Bunch objects in place of the matlab structures.
    """
    out = Bunch()
    fobj = open(fname, 'rb')
    xx = loadmat(fobj)
    keys = [k for k in list(xx.keys()) if not k.startswith("__")]
    for k in keys:
        out[k] = structured_to_bunch(xx[k], masked=masked)
    fobj.close()
    return out

def datenum_to_mpl(dnum):
    "Convert MATLAB datenum(s) to matplotlib datenum(s)."
    return dnum - 366

def datenum_to_day(dnum, yearbase):
    "Convert MATLAB datenum(s) to decimal day relative to a yearbase."
    return dnum - (366 + datetime.date(yearbase, 1, 1).toordinal())

