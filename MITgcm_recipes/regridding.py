import xesmf as xe
import gridfill
import xarray as xr
from MITgcm_recipes import mod_drown_sosie
import numpy as np
import matplotlib.pylab as plt
import scipy.interpolate as spint
from MITgcm_recipes import akima1d

def drown_field(dataset, variable, dims_drown=['lat', 'lon'], mask=None, periodicity=0):

    dataarray = dataset[variable]
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
    
    out = xr.Dataset({dataarray.name: (dataarray.dims,  data_drowned)},
                      coords = dataarray.coords)

    return out


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

def regrid_2_mitgcm_llc(input_dataset, mitgcm_grid, list_variables, 
                        lonname='lon', latname='lat', method='bilinear',
                        faces2blend=[], periodic=True, reuse_weights=True):

    # rename mitgcm grid variables
    target_grid = mitgcm_grid.rename({'XC': 'lon', 'YC': 'lat'})

    iface = target_grid['lon'].get_axis_num('face')
    nfaces = target_grid['lon'].shape[iface]

    # create  a dictionary of interpolators
    regridders = {}
    backup_regridders = {}
    for face in range(nfaces):
        lonmin, lonmax, latmin, latmax = get_bounds(target_grid, face)
        tmpds = input_dataset.sel({lonname: slice(lonmin, lonmax),
                                   latname: slice(latmin, latmax)}) 
        periodic_face = restrict_periodicity(periodic, input_dataset, tmpds, lonname)
        regridders.update({face: xe.Regridder(tmpds,
                                              target_grid.sel(face=face),
                                              method, periodic=periodic_face, 
                                              filename='regrid_face'+str(face)+'.nc',
                                              reuse_weights=reuse_weights)})
        if face in faces2blend:
            backup_regridders.update({face: xe.Regridder(tmpds, 
                                                         target_grid.sel(face=face), 
                                                         'nearest_s2d', periodic=periodic,
                                                         filename='backup_regrid_face'+str(face)+'.nc',
                                                         reuse_weights=reuse_weights)})


    hremapped = xr.Dataset()
    hremapped.update({'XC':target_grid['lon']})
    hremapped.update({'YC':target_grid['lat']})

    for variable in list_variables:
        for face in range(nfaces):
            lonmin, lonmax, latmin, latmax = get_bounds(target_grid, face)
            tmpds = input_dataset[variable].sel({lonname: slice(lonmin, lonmax),
                                                 latname: slice(latmin, latmax)})

            dataface = regridders[face](tmpds)
            if face in faces2blend:
                backup_dataface = backup_regridders[face](tmpds)
                dataface = blend(dataface, backup_dataface, missing=0)

            dataface = dataface.rename({'lon': 'XC', 'lat': 'YC'})
            #stack/concatenate
            if face == 0:
                data_all = dataface.expand_dims(dim='face')
            else:
                data_all = xr.concat([data_all, dataface], dim='face')

        hremapped.update({variable: data_all})

        #hremapped.rename({'lon': 'XC', 'lat': 'YC'}, inplace=True)
    return hremapped


def blend(da1, da2, missing=None):
    tmp1 = da1.values
    tmp2 = da2.values
    mask = np.ones(tmp1.shape)
    mask[np.where(tmp1 == missing)] = 0
    tmp3 = (mask * tmp1) + ((1-mask) * tmp2)
    da3 = xr.DataArray(tmp3, name=da1.name, dims=da1.dims,
                       coords = da1.coords)
    return da3


def get_bounds(target_grid, face):
    lon = get_true_coords(target_grid, 'lon', face)
    lat = get_true_coords(target_grid, 'lat', face)
    lonmin = lon.min() -1
    lonmax = lon.max() +1
    latmin = lat.min() -1
    latmax = lat.max() +1
    return lonmin, lonmax, latmin, latmax


def get_true_coords(grid, coord, face):
    tmp = grid.sel(face=face)[coord].values
    out = np.ma.masked_values(tmp,0)
    return out


def restrict_periodicity(periodicity, full, subset, lonname='lon'):
    indxlonfull = full[lonname].get_axis_num(lonname)
    indxlonsubset = subset[lonname].get_axis_num(lonname)
    sizelonfull = full[lonname].shape[indxlonfull]
    sizelonsubset = subset[lonname].shape[indxlonsubset]
    lonmax_subset = subset[lonname].max()
    lonmin_subset = subset[lonname].min()
    extend = np.abs(lonmax_subset-lonmin_subset)
    res=(full[lonname][1] - full[lonname][0])
    if extend < 360-res:
        periodicity=False
    return periodicity








