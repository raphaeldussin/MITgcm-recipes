import xarray as xr
import numpy as np
from .regridding import dims_of_pointtype

def define_obc_aste(point='T'):
    """ dumb definition of obc position in ASTE 270 model
    """

    xdim, ydim = dims_of_pointtype(point)
    if point == 'V':
        yoffset=1
    else:
        yoffset=0
    # south
    obcs = []
    # obc in south atlantic (24 + 90 of padding -1 C/F)
    obcs.append({'name': 'south', 'face':0, 'saxis': xdim, 'naxis': ydim, 
                 'smin':0, 'smax':270, 'obcpoint': 113+yoffset,
                 'obc_start': 0, 'obc_end': 270})

    # east
    obce = []
    # gibraltar strait
    obce.append({'name': 'east', 'face':1, 'saxis': ydim, 'naxis': xdim,
                 'smin':81, 'smax':83, 'obcpoint': 95,
                 'obc_start': 261, 'obc_end': 263})
    # pacific boundary
    obce.append({'name': 'east', 'face':3, 'saxis': ydim, 'naxis': xdim,
                 'smin':0, 'smax':270, 'obcpoint': 124,
                 'obc_start': 720, 'obc_end': 990})
    # south atlantic boundary
    obce.append({'name': 'east', 'face':5, 'saxis': ydim, 'naxis': xdim,
                 'smin':0, 'smax':270, 'obcpoint': 156,
                 'obc_start': 990, 'obc_end': 1260})

    return obcs, obce

def extract_obc_raw(da, obc_dicts, nxy, nz, nt):
    """ extract the open boundary values from dataarray array according
    to info from list of dicts obc_dicts 

    obc_dict = [{obc}, {obc}]
    with obc = {'name': south, 
                'face':0, 'axis': 'i', 'smin':0, 'smax':270, 'obcpoint': 113, # in xarray world
                'obc_start': 0, 'obc_end': 270} # in obc world

    """

    # create the output array
    out = np.zeros((nt, nz, nxy))

    for obc_dict in obc_dicts:
        extract = da.isel({obc_dict['saxis']: slice(obc_dict['smin'], obc_dict['smax']), 
                           obc_dict['naxis']: obc_dict['obcpoint'],
                           'face': obc_dict['face']}).values
        out[:,:,obc_dict['obc_start']:obc_dict['obc_end']] = extract

    return out
    
