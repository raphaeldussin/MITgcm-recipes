#----- spatial transformation (averages,...)
import numpy as _np
import xarray as _xr
import xgcm

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
    DataArray = DataArray.fillna(0)
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

def avg2d_llc_xr(DataArray):
    ''' return the 3dimensional average of a xarray.DataArray in the LLC geometry 
        xarray version '''
    DataArray = DataArray.fillna(0)
    wa = (DataArray * DataArray['rA'] ).where(DataArray !=0).sum(dim=['face','j','i'])
    norm = (DataArray['rA'] ).where(DataArray !=0).sum(dim=['face','j','i'])
    out = wa / norm
    return out

#----- hovmullers ----

def hov_time_depth_llc(DataArray):
    ''' return the hovmuller of llc data with time vs depth dims '''
    wa = (DataArray * DataArray['rA']).where(DataArray !=0).sum(dim=['face','j','i'])
    norm = (DataArray['rA']).where(DataArray !=0).sum(dim=['face','j','i'])
    out = wa / norm
    return out

#----- velocity rotation ----

def rotate_to_llc(u, v, grid):
    """ rotate eastward/northward u/v velocities to the i/j grid axes
    u : data array for zonal veloicty
    v : data array for meridional velocity
    grid : model grid, must contain CS and SN
    """

    return None

def rotate_llc_to_geo(u, v, grid, face_connections, boundary='extend'):
    """ interp velocities to cell center and rotate
    to geographical axes
    u : data array for zonal velocity
    v : data array for meridional velocity
    grid : model grid, must contain CS and SN
    """
    
    xgrid = xgcm.Grid(grid, face_connections=face_connections)
    uv_center = xgrid.interp_2d_vector({'X': u, 'Y': v}, boundary=boundary)

    u_geo = uv_center['X'] * grid['CS'] - uv_center['Y'] * grid['SN']
    v_geo = uv_center['Y'] * grid['CS'] + uv_center['X'] * grid['SN']

    # this is a wrong but works
    #u_geo = u.rename({'i_g': 'i'}) * grid['CS'] - v.rename({'j_g': 'j'}) * grid['SN']
    #v_geo = v.rename({'j_g': 'j'}) * grid['CS'] + u.rename({'i_g': 'i'}) * grid['SN']

    return u_geo, v_geo

#------- subsetting,...

def extract_box(ds, lonmin=-180, lonmax=180,
                    latmin=-90, latmax=90):

    tmp = ds.where(ds.XC >= lonmin).where(ds.XC <= lonmax)
    out = tmp.where(ds.YC >= latmin).where(ds.YC <= latmax)

    return out
