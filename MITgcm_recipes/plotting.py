#-------- plotting functions -------------------
import matplotlib.pylab as _plt
import numpy as _np
import cartopy as _cart
from matplotlib import colors as _colors
from gridfill import fill as _fill

# plot ASTE on the best possible projection
def plot_ASTE(dataarray, dict_plt):
    ''' make a plot for ASTE from a xmitgcm dataarray
    '''
    
    # explicit list of keys
    figsize = dict_plt['figsize']
    vmin    = dict_plt['vmin']
    vmax    = dict_plt['vmax']
    cmap    = dict_plt['cmap']

    fig = _plt.figure(figsize=figsize)
    m = _plt.axes(projection=_cart.crs.Orthographic(central_longitude=-35, central_latitude=40))
    m.add_feature(_cart.feature.LAND, facecolor='0.75')
    gl = m.gridlines(draw_labels=False)
    norm = _colors.Normalize(vmin=vmin, vmax=vmax)
    #contours = _np.arange(vmin,vmax,0.2)
    

    kw = dict(eps=1e-2, relax=0.6, itermax=1e4, initzonal=False,
              cyclic=False, verbose=False)
    
    # first plot the faces without padding
    for f in [1,2,4]:
        lon = dataarray['XC'].sel(face=f).values
        lat = dataarray['YC'].sel(face=f).values
        # fill missing values in lon/lat arrays
        lonplt, converged = _fill(_np.ma.masked_values(lon,0), 1, 0, **kw)
        latplt, converged = _fill(_np.ma.masked_values(lat,0), 1, 0, **kw)
        # mask data where needed
        data = dataarray.sel(face=f).values
        depth= dataarray['Depth'].sel(face=f).values
        mask = _np.full(depth.shape, False)
        mask[_np.where(depth == 0)] = True
        dataplt = _np.ma.masked_array(data=data,mask=mask)
        
        C = m.pcolormesh(lonplt, latplt,
                         dataplt, norm=norm, cmap=cmap,\
                         transform=_cart.crs.PlateCarree())
        
    for f in [3]:
        lon = dataarray['XC'].sel(face=f, i=slice(0,125)).values
        lat = dataarray['YC'].sel(face=f, i=slice(0,125)).values
        # mask data where needed
        data = dataarray.sel(face=f, i=slice(0,125)).values
        depth= dataarray['Depth'].sel(face=f, i=slice(0,125)).values
        mask = _np.full(depth.shape, False)
        mask[_np.where(depth == 0)] = True
        dataplt = _np.ma.masked_array(data=data,mask=mask)
        
        C = m.pcolormesh(lon, lat,
                         dataplt, norm=norm, cmap=cmap,\
                         transform=_cart.crs.PlateCarree())
        
    for f in [5]:
        lon = dataarray['XC'].sel(face=f, i=slice(0,155)).values
        lat = dataarray['YC'].sel(face=f, i=slice(0,155)).values
        # mask data where needed
        data = dataarray.sel(face=f, i=slice(0,155)).values
        depth= dataarray['Depth'].sel(face=f, i=slice(0,155)).values
        mask = _np.full(depth.shape, False)
        mask[_np.where(depth == 0)] = True
        dataplt = _np.ma.masked_array(data=data,mask=mask)
        
        C = m.pcolormesh(lon, lat,
                         dataplt, norm=norm, cmap=cmap,\
                         transform=_cart.crs.PlateCarree())
        
    for f in [0]:
        lon = dataarray['XC'].sel(face=f, j=slice(115,270)).values
        lat = dataarray['YC'].sel(face=f, j=slice(115,270)).values
        # mask data where needed
        data = dataarray.sel(face=f, j=slice(115,270)).values
        depth= dataarray['Depth'].sel(face=f, j=slice(115,270)).values
        mask = _np.full(depth.shape, False)
        mask[_np.where(depth == 0)] = True
        dataplt = _np.ma.masked_array(data=data,mask=mask)
        
        C = m.pcolormesh(lon, lat,
                         dataplt, norm=norm, cmap=cmap,\
                         transform=_cart.crs.PlateCarree())
        _plt.colorbar(C, norm=norm, shrink=0.75)
        
    gl = m.gridlines(draw_labels=False)

    return fig
        


