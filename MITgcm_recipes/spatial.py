#----- spatial transformation (averages,...)
import numpy as _np
import xarray as _xr

def avg3d_llc(DataArray, darea=None, dz=None):
    ''' return the 3dimensional average of a xarray.DataArray in the LLC geometry '''
    # MITgcm uses O as missing value
    data  = DataArray.where(DataArray !=0).to_masked_array()  # convert to numpy.masked_array
    mask  = _np.ones(data.shape)         # init the mask to one
    if darea is None:
        darea = DataArray['rA'].values   # get cell area
    if dz is None:
        dz = DataArray['drF'].values     # get layer thickness

    # fill values of binary mask with logical mask
    # where values are masked in np.masked_array is becoming zero
    mask[data.mask] = 0

    # we assume time axis is present and in first dimension
    # sum first on layers (i,j,face) then along depth
    out = ((data * darea * mask).sum(axis=4).sum(axis=3).sum(axis=2) * dz).sum(axis=1) / \
          ((darea * mask).sum(axis=4).sum(axis=3).sum(axis=2) * dz).sum(axis=1)

    return out

def avg3d_llc_xr(DataArray):
    ''' return the 3dimensional average of a xarray.DataArray in the LLC geometry 
        xarray version '''
    wa = (DataArray * DataArray['rA'] * DataArray['drF']).where(DataArray !=0).sum(dim=['face','k','j','i'])
    norm = (DataArray['rA'] * DataArray['drF']).where(DataArray !=0).sum(dim=['face','k','j','i'])
    out = wa / norm
    return out

def avg2d_llc(DataArray, darea=None):
    ''' return the 2dimensional average of a xarray.DataArray in the LLC geometry '''
    # MITgcm uses O as missing value
    data  = DataArray.where(DataArray !=0).to_masked_array()  # convert to numpy.masked_array
    mask  = _np.ones(data.shape)         # init the mask to one
    if darea is None:
        darea = DataArray['rA'].values   # get cell area

    # fill values of binary mask with logical mask
    # where values are masked in np.masked_array is becoming zero
    mask[data.mask] = 0

    # we assume time axis is present and in first dimension
    # sum on layers (i,j,face) 
    out = (data * darea * mask).sum(axis=3).sum(axis=2).sum(axis=1) / \
          (darea * mask).sum(axis=3).sum(axis=2).sum(axis=1) 

    return out
