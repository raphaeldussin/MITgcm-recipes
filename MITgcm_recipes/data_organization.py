import xarray as _xr
import numpy as _np

def facets2faces_aste270(target_da, variable, facet1=None, facet2=None, facet3=None,
                         facet4=None, facet5=None, xdim='x', ydim='y'):

    """ create a xarray.DataArray using the template target_da for coords and dimensions,
        then add the different facets to the dataarray 

    target_da : data array we want to mimick the shape
    variable   : name of variable in facets netcdf

    """
    
    zeros = _np.zeros(target_da.shape)

    # create a dataarray and initialize with zeros
    da = _xr.DataArray(zeros, coords=target_da.coords, 
                       dims=target_da.dims, name=variable, attrs=None)
    
    if facet1 is not None:
        dsobs1 = _xr.open_dataset(facet1)
        # southeast atlantic 
        da.loc[dict(face=0, j=slice(90,270))] = dsobs1[variable].sel(y=dsobs1[ydim][slice(0,180)])
        # northeast atlantic
        da.loc[dict(face=1)] = dsobs1[variable].sel(y=dsobs1[ydim][slice(180,450)])

    if facet3 is not None:
        dsobs3 = _xr.open_dataset(facet3)
        # arctic
        da.loc[dict(face=2)] = dsobs3[variable]

    if facet4 is not None:
        dsobs4 = _xr.open_dataset(facet4)
        # North Pacific
        da.loc[dict(face=3, i=slice(0,180-1))] = dsobs4[variable]

    if facet5 is not None:
        dsobs5 = _xr.open_dataset(facet5)
        # northwest atlantic
        da.loc[dict(face=4)] = dsobs5[variable].sel(x=dsobs5[xdim][slice(0,270)])
        # southwest atlantic
        da.loc[dict(face=5, i=slice(0,180-1))] = dsobs5[variable].sel(x=dsobs5[xdim][slice(270,450)])
        
    return da

def build_ACgrid_aste270(facet1_grid='facet1.nc', facet5_grid='facet5.nc'):
    """ build an Atlantic Common grid for ASTE based on facet 1 and 5 + extrapolation """

    lon_out = _np.zeros((540, 540))
    lat_out = _np.zeros((540, 540))

    grid1 = _xr.open_dataset(facet1_grid)
    grid5 = _xr.open_dataset(facet5_grid)

    # NorthWest Atlantic
    lon_out[90:,:270] = grid5['XC'].values.transpose()[::-1,:]
    lat_out[90:,:270] = grid5['YC'].values.transpose()[::-1,:]
    # NorthEast Atlantic
    lon_out[90:,270:] = grid1['XC'].values
    lat_out[90:,270:] = grid1['YC'].values
    # Southward expansion
    # for longitude we assume = to the first 90 rows
    # for latitude we assume same dy
    # note that there is no data so this does not matter much
    # as long as we don't give zeros (f up the plots)
    lon_out[:90,:] = lon_out[90:180,:]
    dlat = lat_out[91,:] - lat_out[90,:]
    _, incr = _np.meshgrid(dlat, _np.arange(-90,0))
    lat_out[:90,:] = lat_out[90,:] + dlat * incr

    # create dataset
    ds = _xr.Dataset({},
                     coords={'XC': (['y','x'], lon_out),
                             'YC': (['y','x'], lat_out)} )
    return ds

def regrid_to_AC(dataarray, dsAC):
    """ take a xarray.dataarray in the llc ASTE format and regrid
    to the atlantic common grid """

    data_out = _np.zeros((540, 540))

    # Face 0
    data_out[:270,270:] = dataarray.sel(face=0).values
    # Face 1
    data_out[270:,270:] = dataarray.sel(face=1).values
    # Face 4
    data_out[270:,:270] = dataarray.sel(face=4).values.transpose()[::-1,:]
    # Face 5
    data_out[:270,:270] = dataarray.sel(face=5).values.transpose()[::-1,:]

    out = _xr.Dataset({dataarray.name: (('y', 'x'), data_out)},
                     coords={'XC': (['y','x'], dsAC['XC'].values),
                             'YC': (['y','x'], dsAC['YC'].values)} )
    return out

def build_PACgrid_aste270(facet4_grid='facet4.nc'):
    """ build a grid for the Pacific, expanding lon, lat and patching zeros """

    lon_out = _np.zeros((270,270))
    lat_out = _np.zeros((270,270))

    grid4 = _xr.open_dataset(facet4_grid)

    # North Pacific
    lon_out[90:,:] = grid4['XC'].values.transpose()[::-1,:]
    lat_out[90:,:] = grid4['YC'].values.transpose()[::-1,:]
    # Southward expansion
    lon_out[:90,:] = lon_out[90:180,:]
    dlat = lat_out[91,:] - lat_out[90,:]
    _, incr = _np.meshgrid(dlat, _np.arange(-90,0))
    lat_out[:90,:] = lat_out[90,:] + dlat * incr

    # remove discontinuity by having all lon negative
    lon_out[_np.where(lon_out >=0)] -= 360

    # create dataset
    ds = _xr.Dataset({},
                     coords={'XC': (['y','x'], lon_out),
                             'YC': (['y','x'], lat_out)} )
    return ds

def regrid_to_PAC(dataarray,dsPAC):
    """ take a ASTE dataarray and regrid for the Pacific """

    data_out = _np.zeros((270, 270))

    # Face 3
    data_out[:,:] = dataarray.sel(face=3).values.transpose()[::-1,:]

    out = _xr.Dataset({dataarray.name: (('y', 'x'), data_out)},
                     coords={'XC': (['y','x'], dsPAC['XC'].values),
                             'YC': (['y','x'], dsPAC['YC'].values)} )
    return out

def build_ARCgrid_aste270(facet3_grid='facet3.nc'):
    """ build an Arctic grid (basically to patch zeros in XC, YC) """

    lon_out = _np.zeros((270,270))
    lat_out = _np.zeros((270,270))

    grid3 = _xr.open_dataset(facet3_grid)

    # Arctic
    lon_out[:,:] = grid3['XC'].values
    lat_out[:,:] = grid3['YC'].values

        # create dataset
    ds = _xr.Dataset({},
                     coords={'XC': (['y','x'], lon_out),
                             'YC': (['y','x'], lat_out)} )
    return ds

def regrid_to_ARC(dataarray,dsARC):
    """ take a ASTE dataarray and regrid for the Arctic """

    data_out = _np.zeros((270, 270))

    # Face 2
    data_out[:,:] = dataarray.sel(face=2).values

    out = _xr.Dataset({dataarray.name: (('y', 'x'), data_out)},
                     coords={'XC': (['y','x'], dsARC['XC'].values),
                             'YC': (['y','x'], dsARC['YC'].values)} )
    return out
