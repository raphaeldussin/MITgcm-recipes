import xesmf as xe
import gridfill
import xarray as xr
from MITgcm_recipes import mod_drown_sosie, akima1d
import numpy as np
import matplotlib.pylab as plt
import scipy.interpolate as spint
from MITgcm_recipes import akima1d

#--------------------
# DROWN
#--------------------

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
        data_drowned = drown_2d_field_sosie(data) #+options
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


#--------------------
# REGRID
#--------------------

def regrid_2_mitgcm_llc(input_dataset, mitgcm_grid, list_variables, point='T',
                        lonname='lon', latname='lat', method='bilinear',
                        faces2blend=[], blend_mask={},
                        periodic=True, reuse_weights=True,
                        regridname='regrid_face'):

    # rename mitgcm grid variables
    if point == 'T':
        target_grid = mitgcm_grid.rename({'XC': 'lon', 'YC': 'lat'})
    elif point == 'U':
        target_grid = mitgcm_grid.rename({'XG': 'lon', 'YC': 'lat'})
        target_grid['lon'] = target_grid['lon'].rename({'j_g':'j'})
    elif point == 'V':
        target_grid = mitgcm_grid.rename({'XC': 'lon', 'YG': 'lat'})
        target_grid['lon'] = target_grid['lon'].rename({'j':'j_g'}) # needed by xesmf to get dim right
        target_grid['lat'] = target_grid['lat'].rename({'i_g':'i'})


    coords = create_coords_regridded(mitgcm_grid, point=point)
    nfaces = len(coords['face'].values)

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
                                              filename=regridname+str(face)+'.nc',
                                              reuse_weights=reuse_weights)})
        if face in faces2blend:
            print('Creating nearest neighbor weights for face', face)
            backup_regridders.update({face: xe.Regridder(tmpds, 
                                                         target_grid.sel(face=face),
                                                         'nearest_s2d', periodic=periodic_face,
                                                         filename='backup_'+regridname+str(face)+'.nc',
                                                         reuse_weights=reuse_weights)})


    hremapped = xr.Dataset()

    for variable in list_variables:
        for face in range(nfaces):
            lonmin, lonmax, latmin, latmax = get_bounds(target_grid, face)
            tmpds = input_dataset[variable].sel({lonname: slice(lonmin, lonmax),
                                                 latname: slice(latmin, latmax)})

            dataface = regridders[face](tmpds)
            if face in faces2blend:
                print('running nn regridding for face', face)
                backup_dataface = backup_regridders[face](tmpds)
                #dataface = blend(dataface, backup_dataface, missing=0)
                dataface = blend_using_mask(dataface, backup_dataface, blend_mask[face])

            # rewrite the coordinates
            dataface = rewrite_coords_regridded(dataface, coords, point=point)

            #stack/concatenate
            if face == 0:
                data_all = dataface.expand_dims(dim='face')
            else:
                data_all = xr.concat([data_all, dataface], dim='face')

        # once concat has happened, add faces
        data_all.coords['face'] = coords['face']
        # add to dataset
        hremapped.update({variable: data_all})

    return hremapped


#--------------------
# VERTICAL interpol
#--------------------

def vertical_interpolation(inputds, depth_out, list_variables,
                           lonvar='lon', latvar='lat',
                           depth_varin='depth',
                           timevar=None):

    """
    depth_out: data_array

    """
    # We need to create a new depth coordinate for the dataset
    outputds = xr.Dataset()
    outputds.update({lonvar:inputds[lonvar]})
    outputds.update({latvar:inputds[latvar]})
    if timevar is not None:
        outputds.update({timevar:inputds[timevar]})

    outputds.update({'Z':depth_out})
    # input depth vector
    z_in = inputds[depth_varin].values
    # output depth vector
    z_out = depth_out.values
    nz = len(z_out)

    # test for direction (xmitgcm has depth negative)
    if z_out.min() < 0:
        z_out = - z_out

    for variable in list_variables:
        idx = inputds[variable].get_axis_num(lonvar)
        idy = inputds[variable].get_axis_num(latvar)
        nx = inputds[variable].shape[idx]
        ny = inputds[variable].shape[idy]
        if timevar is not None:
            idt = inputds[variable].get_axis_num(timevar)
            nt = inputds[variable].shape[idt]
            data_out = np.empty((nt, nz, ny, nx))
            for kt in range(nt):
                data_in = inputds[variable].isel({timevar: kt}).values
                data_out[kt,:,:,:] = akima1d.mod_akima_1d.vertical_interpolation(z_in, data_in, z_out)
        else:
            data_out = np.empty((nz, ny, nx))
            data_in = inputds[variable].values
            data_out[:] = akima1d.mod_akima_1d.vertical_interpolation(z_in, data_in, z_out)

        if timevar is not None:
            dims_out = [timevar, 'k', latvar, lonvar]
        else:
            dims_out = ['k', latvar, lonvar]

        outputds.update({variable: xr.DataArray(data_out, dims=dims_out)})

    return outputds

#--------------------
# Masking
#--------------------

def mask_output(dataarray, mitgcm_grid, point='T'):
    ''' mask the output from the regridding with mask from the model '''
    if point == 'T':
        fac = 'hFacC'
    elif point == 'U':
        fac = 'hFacW'
    elif point == 'V':
        fac = 'hFacS'
    if 'k' in dataarray.dims:
        out = dataarray.where(mitgcm_grid[fac] != 0)
    else:
        out = dataarray.where(mitgcm_grid[fac].sel(k=0) != 0)
    return out

def mask_output_zeros(dataarray, mitgcm_grid, point='T'):
    ''' mask the output from the regridding with mask from the model '''
    tmp = mask_output(dataarray, mitgcm_grid, point=point)
    out = tmp.fillna(0)
    return out



#--------------------
# helping functions
#--------------------

def blend(da1, da2, missing=None):
    ''' blend dataarray1 with values from dataarray2 where
    dataarray1 has missing value '''
    tmp1 = da1.values
    tmp2 = da2.values
    mask = np.ones(tmp1.shape)
    mask[np.where(tmp1 == missing)] = 0
    tmp3 = (mask * tmp1) + ((1-mask) * tmp2)
    da3 = xr.DataArray(tmp3, name=da1.name, dims=da1.dims,
                       coords = da1.coords)
    return da3

def blend_using_mask(da1, da2, mask):
    ''' blend dataarray1 with values from dataarray2
    mask should be zero where values needed to be blended, 1 otherwise'''
    tmp1 = da1.values
    tmp2 = da2.values
    tmp3 = (mask * tmp1) + ((1-mask) * tmp2)
    da3 = xr.DataArray(tmp3, name=da1.name, dims=da1.dims,
                       coords = da1.coords)
    return da3

def get_bounds(target_grid, face):
    ''' returns lon/lat min/max for llc face '''
    lon = get_true_coords(target_grid, 'lon', face)
    lat = get_true_coords(target_grid, 'lat', face)
    lonmin = lon.min() -1
    lonmax = lon.max() +1
    latmin = lat.min() -1
    latmax = lat.max() +1
    return lonmin, lonmax, latmin, latmax


def get_true_coords(grid, coord, face):
    ''' returns lon/lat while filtering zeros = missing value'''
    tmp = grid.sel(face=face)[coord].values
    out = np.ma.masked_values(tmp,0)
    return out


def restrict_periodicity(periodicity, full, subset, lonname='lon'):
    ''' decide if original periodicity has been broken by subsetting array '''
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


def geo_roll_360_to_180(ds, londim):
    ind = (np.abs(ds[londim]-180.)).argmin().values
    ds = ds.roll(**{londim: -ind})
    datalon=ds[londim].values
    datalon[:ind] -=360
    ds[londim] = datalon
    return ds


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


def create_coords_regridded(mitgcm_grid, point='T'):
    ''' generate new set of coords in the xmitgcm style
    for regridded array '''

    xdim, ydim = dims_of_pointtype(point)
    lon, lat = geo_of_pointtype(point)

    # create horizontal coordinates "i,j"
    ix = mitgcm_grid[lon].get_axis_num(xdim)
    nx = mitgcm_grid[lon].shape[ix]

    iy = mitgcm_grid[lat].get_axis_num(ydim)
    ny = mitgcm_grid[lat].shape[iy]

    x = xr.DataArray(np.arange(nx), dims=[xdim])
    y = xr.DataArray(np.arange(ny), dims=[ydim])

    # create face coordinate 
    iface = mitgcm_grid[lon].get_axis_num('face')
    nfaces = mitgcm_grid[lon].shape[iface]
    
    faces = xr.DataArray(np.arange(nfaces), dims=['face'])

    coords = {xdim: x, ydim: y, 'face': faces}

    return coords


def rewrite_coords_regridded(da, coords, point='T'):
    ''' rewrite the dataarray with xmitgcm coords '''

    xdim, ydim = dims_of_pointtype(point)

    # rewrite the coordinates
    da.rename({'lon': xdim})
    da.rename({'lat': ydim})
    da.coords[xdim] = coords[xdim]
    da.coords[ydim] = coords[ydim]
    da = da.drop('lon')
    da = da.drop('lat')

    return da


def dims_of_pointtype(point):
    ''' returns the dims pair for each point type (e.g. T, U, V) '''

    if point == 'T':
        xdim='i' ; ydim='j'
    elif point == 'U':
        xdim='i_g' ; ydim='j'
    elif point == 'V':
        xdim='i' ; ydim='j_g'
    else:
        raise ValueError('no such point, correct values are T, U and V')
    return xdim, ydim

def geo_of_pointtype(point):
    ''' returns the geo (lon/lat) pair for each point type (e.g. T, U, V) '''

    if point == 'T':
        lon = 'XC' ; lat = 'YC'
    elif point == 'U':
        lon = 'XG' ; lat = 'YC'
    elif point == 'V':
        lon = 'XC' ; lat = 'YG'
    else:
        raise ValueError('no such point, correct values are T, U and V')
    return lon, lat



