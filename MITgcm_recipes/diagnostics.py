import xarray as xr
import numpy as np

def mld_rho_003(ds, saltname='SALT', tempname='THETA'):
    ''' compute the mixed layer depth with the drho = 0.03 kg/m3 criterion '''
    # get density
    dens = compute_potdens(ds, saltname=saltname, tempname=tempname)
    # compute anomaly to 10 meters density
    anom_dens = dens - dens.swap_dims({'k': 'Z'}).interp(Z=-10)
    # return depth of last level satisfying condition
    mld = - ds['Z'].where(dens !=0).where(anom_dens>=0).where(anom_dens < 0.03).min(dim='k')
    return mld

def mld_temp_02(ds, tempname='THETA'):
    ''' compute the mixed layer depth with the dT = 0.2 degC criterion '''
    #temp = ds[tempname].load().chunk({'time':1, 'face':1})
    #dT = temp - temp.swap_dims({'k': 'Z'}).interp(Z=-10)
    #mld = - ds['Z'].where(dT !=0).where(dT<=0).where(dT > -0.2).min(dim='k')
    return mld

#def compute_dens(ds, saltname='SALT', tempname='THETA'):
#    import seawater
#    # needs to recompute in-situ temperature
#    temp = seawater.temp(ds[saltname], ds[tempname], ds['Z'])
#    dens = seawater.dens0(ds[saltname], temp)
#    dens[dens < 1000] = 0
#    out = xr.DataArray(dens, dims=('time','k','face','j','i') )
#    out = out.assign_coords(XC=ds['XC'], YC=ds['YC'], Z=ds['Z'])
#    return out

def compute_potdens(ds, saltname='SALT', tempname='THETA'):
    import gsw
    """ compute the potential density
    """
    # compute the Conservative Temperature from the model's potential temperature
    temp = ds[tempname].transpose(*('time','k','face','j','i'))
    salt = ds[tempname].transpose(*('time','k','face','j','i'))
    CT = gsw.CT_from_pt(salt, temp)
    z, lat = xr.broadcast(ds['Z'], ds['YC'])
    z = z.transpose(*('k','face','j','i'))
    lat = lat.transpose(*('k','face','j','i'))
    # compute pressure from depth
    p = gsw.p_from_z(z, lat)
    # compute in-situ temperature
    T = gsw.t_from_CT(salt, CT, p)
    # compute potential density
    rho = gsw.pot_rho_t_exact(salt, T, p, 0.)
    # create new dataarray
    darho = xr.full_like(temp,0.)
    darho = darho.load().chunk({'time':1, 'face':1})
    darho.name = 'RHO'
    darho.attrs['long_name'] = 'Potential Density ref at 0m'
    darho.attrs['standard_name'] = 'RHO'
    darho.attrs['units'] = 'kg/m3'
    darho.values = rho
    # filter special value
    darho = darho.where(darho > 1000)
    darho = darho.assign_coords(XC=ds['XC'], YC=ds['YC'], Z=ds['Z'])
    return darho
