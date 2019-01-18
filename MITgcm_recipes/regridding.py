import xesmf as xe
import gridfill
import xarray as xr
from MITgcm_recipes import mod_drown_sosie
import numpy as np
import matplotlib.pylab as plt
import scipy.interpolate as spint
from MITgcm_recipes import akima1d

def drown_field(dataarray, dims_drown=['lat', 'lon'], mask=None, periodicity=0):

    # get the dimensions on which to loop
    dims_da = list(dataarray.dims)
    dims_loop = dims_da.copy()
    for d in dims_drown:
        dims_loop.remove(d)

    nloopdims = len(dims_loop)

    # init to original to have the structure right
    data_drowned = np.empty(dataarray.shape)

    if nloopdims == 0:
        data = dataarray.values
        drowned = drown_2d_field_sosie(data) #+options
    # I can't figure out a smarter way of doing this
    # but I'd love to make embedded loops invisible
    elif nloopdims == 1:
        index0 = dataarray.get_axis_num(dims_loop[0])
        nloop0 = dataarray.shape[index0]
        for kloop0 in range(nloop0):
            dataslice = dataarray.isel({dims_loop[0]: kloop0}).values
            drownedslice = drown_2d_field_sosie(dataslice) #+options
            if index0 == 0:
                data_drowned[kloop0,:,:] = drownedslice
            elif index0 == 1:
                data_drowned[:,kloop0,:] = drownedslice
            elif index0 == 2:
                data_drowned[:,:,kloop0] = drownedslice
    elif nloopdims == 2:    
        index0 = dataarray.get_axis_num(dims_loop[0])
        nloop0 = dataarray.shape[index0]
        index1 = dataarray.get_axis_num(dims_loop[1])
        nloop1 = dataarray.shape[index1]
        for kloop0 in range(nloop0):
            for kloop1 in range(nloop1):
                dataslice = dataarray.isel({dims_loop[0]: kloop0,
                                            dims_loop[1]: kloop1,}).values
                drownedslice = drown_2d_field_sosie(dataslice) #+options
                # there must be a more elegant way
                if index0 == 0:
                    if index1 == 1:
                        data_drowned[kloop0,kloop1,:,:] = drownedslice
                    if index1 == 2:
                        data_drowned[kloop0,:,kloop1,:] = drownedslice
                    if index1 == 3:
                        data_drowned[kloop0,:,:,kloop1] = drownedslice
                elif index0 == 1:
                    if index1 == 0:
                        data_drowned[kloop1,kloop0,:,:] = drownedslice
                    if index1 == 2:
                        data_drowned[:,kloop0,kloop1,:] = drownedslice
                    if index1 == 3:
                        data_drowned[:,kloop0,:,kloop1] = drownedslice
                elif index0 == 2:
                    if index1 == 0:
                        data_drowned[kloop1,:,kloop0,:] = drownedslice
                    if index1 == 1:
                        data_drowned[:,kloop1,kloop0,:] = drownedslice
                    if index1 == 3:
                        data_drowned[:,:,kloop0,kloop1] = drownedslice
                elif index0 == 3:
                    if index1 == 0:
                        data_drowned[kloop1,:,:,kloop0] = drownedslice
                    if index1 == 1:
                        data_drowned[:,kloop1,:,kloop0] = drownedslice
                    if index1 == 2:
                        data_drowned[:,:,kloop1,kloop0] = drownedslice
                else:
                    raise ValueError('not possible')
    
    dataarray_out = xr.DataArray(data_drowned, coords=dataarray.coords, dims=dataarray.dims)

    return dataarray_out


def drown_2d_field_sosie(array, mask=None, spval=None, periodicity=0, nb_inc=200, nb_smooth=2):
    """ fill land values from array, with optional mask, periodicity, ... """

    if len(array.shape) != 2:
        raise ValueError('array should be 2d')

    if mask is not None:
        if mask.shape != data.shape:
            raise ValueError('mask should be 2d')
        if mask.min() < 0 or mask.max() > 1:
            raise ValueError('mask values should be 0/1')
    else:
        mask = int_mask_from_missing_value(array, spval=spval)

    # set expected fill value where mask is land
    tmpin = array.copy()
    tmpin[np.where(mask == 0)] = -9999.

    drowned = mod_drown_sosie.mod_drown.drown(periodicity, tmpin.transpose(),\
                                              mask.transpose(), nb_inc=nb_inc,\
                                              nb_smooth=nb_smooth)

    out = drowned.transpose()
    return out



def drown_2d_field_gridfill(array, mask=None, spval=None, periodic=True, itermax=100, relax=.6):
    """ fill land values from array, with optional mask, periodicity, ... """

    if len(array.shape) != 2:
        raise ValueError('array should be 2d')

    if mask is not None:
        if mask.shape != data.shape:
            raise ValueError('mask should be 2d')
        if mask.min() < 0 or mask.max() > 1:
            raise ValueError('mask values should be 0/1')
    else:
        mask = logic_mask_from_missing_value(array, spval=spval)

    grids = np.ma.masked_array(data=array, mask=mask)

    xdim=1
    ydim=0
    eps=1e-4
    
    drowned, conv = gridfill.fill(grids, xdim, ydim, eps, relax=.6, itermax=100, initzonal=True,\
                                  cyclic=periodic, verbose=False)

    return drowned


def int_mask_from_missing_value(array, spval=None):
    ''' create binary 1/0 mask for masked values in array '''

    mask = np.ones(array.shape, dtype=np.int)
    if spval is None: # nan by default
        mask[np.isnan(array)] = 0
    else:
        mask[np.where(array == spval)] = 0

    return mask


def logic_mask_from_missing_value(array, spval=None):
    ''' create True/False mask for masked values in array '''

    if spval is None: # nan by default
        mask = np.isnan(array)
    else:
        mask = np.where(array == spval)

    return mask













