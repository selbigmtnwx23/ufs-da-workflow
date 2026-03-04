#!/usr/bin/env python3

import os, sys
import logging
import netCDF4 as nc
import numpy as np
from scipy.ndimage import gaussian_filter, distance_transform_edt
import warnings
import yaml


# Main part (will be called at the end) ============================= CHJ =====
def main():

    yaml_file = "calc_scales4parameter.yaml"
    with open(yaml_file, 'r') as f:
        yaml_data = yaml.load(f, Loader=yaml.FullLoader)
    f.close()

    bkg_fp = yaml_data['bkg_fp']
    gridspec_fn = yaml_data['gridspec_fn']
    HZ_MAX = yaml_data['HZ_MAX']
    HZ_MIN_GRID_MULT = yaml_data['HZ_MIN_GRID_MULT']
    HZ_ROSSBY_MULT = yaml_data['HZ_ROSSBY_MULT']
    mld_fp = yaml_data['mld_fp']
    mld_vn = yaml_data['mld_vn']
    output_fn = yaml_data['output_fn']
    output_variable_hz = yaml_data['output_variable_hz']
    output_variable_vt = yaml_data['output_variable_vt']
    PY_LOG_LEVEL = yaml_data['PY_LOG_LEVEL']
    VT_MIN = yaml_data['VT_MIN']
    VT_MAX = yaml_data['VT_MAX']
    work_dir = yaml_data['work_dir']

    # Set logging config
    log_level_str = PY_LOG_LEVEL.upper()
    try:
        log_level = getattr(logging, log_level_str)
    except AttributeError:
        log_level_str = "INFO"
        log_level = logging.INFO
        print(f''' WARNING: Invalid log level "{PY_LOG_LEVEL.upper()}", set to INFO.''')
    print(f''' Python Log Level = str: {log_level_str}, attr: {log_level}''')
    logging.basicConfig(format='%(levelname)s::%(pathname)s::L%(lineno)d::%(message)s', level=log_level)
    logging.info(f''' YAML Data: {yaml_data}''')

    if isinstance(HZ_MAX, str):
        HZ_MAX = float(HZ_MAX)

    # read input
    gridspec_fp = os.path.join(work_dir,gridspec_fn)
    with nc.Dataset(gridspec_fp, 'r') as src:
        rossbyRadius = src.variables['rossby_radius'][0]
        dx = src.variables['dx'][0]
        dy = src.variables['dy'][0]
        area = src.variables['area'][0]
        mask = src.variables['mask2d'][0]

    with nc.Dataset(bkg_fp, 'r') as src:
        h = src.variables['h'][0]
        varType = src.variables['h'].datatype

    with nc.Dataset(mld_fp, 'r') as src:
        mld = src.variables[mld_vn][0]
    nl3 = mld.shape
  
    nz = h.shape[0]
    ny = h.shape[1]
    nx = h.shape[2]
  
    # calculate horizontal scales
    hz_scales = rossbyRadius * HZ_ROSSBY_MULT
    hz_scales = np.clip(hz_scales, a_min=dx*HZ_MIN_GRID_MULT, a_max= HZ_MAX)
    hz_scales = np.clip(hz_scales, a_min=dy*HZ_MIN_GRID_MULT, a_max= HZ_MAX)
  
    # calculate hz smoothing scales for use inside this script.
    # The scales (in units of # of grid cells) is zonally averaged, we assume hz scales
    # do not vary as much with longitude as they do with latitude
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        smoothingScale = np.nanmean(np.where(mask, hz_scales/np.sqrt(area), np.NaN), axis=1)
    smoothingScale[0] = smoothingScale[1]
    smoothingScale[-1] = smoothingScale[-2]
  
    # smooth the horizontal scales WITH the horizontal scales.
    hz_scales_orig = hz_scales.copy()
    for j in range(hz_scales.shape[0]): # zonally for efficiency (otherwise it would be am i/j set of for loops)
        hz_scales[j,:] = gaussian_filter(hz_scales_orig, sigma=smoothingScale[j], mode='nearest')[j,:]
    hz_scales = np.where(mask, hz_scales, 0)
    del hz_scales_orig
  
    # calculate vertical scales.
    # ---------------------------------------------------------
  
    # initialize with zeros
    vtScales = np.zeros_like(h)
  
    # smooth the MLD with the horizontal scales.
    mld_orig = mld.copy()
    # fill missing
    mld_orig = mld[tuple(distance_transform_edt(mask == 0, return_distances=False, return_indices=True))]
    for j in range(mld.shape[0]): # zonally for efficiency (otherwise it would be am i/j set of for loops)
        mld[j,:] = gaussian_filter(mld_orig, sigma=smoothingScale[j], mode='nearest')[j,:]
    mld = np.where(mask, mld, 0)
    del mld_orig
  
    # calculate layer depths for use later
    layerDepth = np.cumsum(h, axis=0) - h/2.0
    maxLevels = np.sum(np.where(h > 0.01, 1, 0), axis=0) # number of lvls before we hit bottom
  
    # find index of last level in the mixed layer
    mlLastLevel = np.where(layerDepth < np.stack([mld]*nz, axis=0), 1, 0)
    mlLastLevel[h <= 0.01] = 0 # don't count small bottom layers
    mlLastLevel = np.clip(np.sum(mlLastLevel, axis=0)-1, a_min=0, a_max=nz-2)
  
    # add to that a fraction of layer at the bottom of the mixed layer
    # for a total fractional number of levels in the ML
    y_idx, x_idx = np.indices(mlLastLevel.shape)
    depth1 = layerDepth[mlLastLevel, y_idx, x_idx]
    depth2 = layerDepth[mlLastLevel+1, y_idx, x_idx]
    mlLayers = mlLastLevel + (mld - depth1) / (depth2-depth1)
    mlLayers = np.clip(mlLayers, a_min=1,a_max=maxLevels)
  
    # set the top to the number of levels in the ML, interpolate down to bottom of ML
    # AND enforce a min/max value. The max value is important in order
    # to keep the number of diffusion iterations in check. The resulting diracs
    # are not quite the same... but close enough.
    layer, _, _ = np.indices(h.shape)
    mlLayers3D = np.stack([mlLayers]*nz, axis=0)
    vtScales[:] = np.clip(mlLayers3D - layer, a_min=VT_MIN, a_max=VT_MAX)
  
    # ignore thin layers at the bottom
    vtScales[h <= 0.01] = 0
  
    # write output file
    output_fp = os.path.join(work_dir,output_fn)
    with nc.Dataset(output_fp, 'w') as dst:
        dst.createDimension('Time', 1)
        dst.createDimension('zaxis_1', vtScales.shape[0] )
        dst.createDimension('yaxis_1', vtScales.shape[1] )
        dst.createDimension('xaxis_1', vtScales.shape[2] )
  
        time = dst.createVariable('Time', varType, ('Time',))
        var_vt = dst.createVariable(output_variable_vt, varType, ('Time','zaxis_1','yaxis_1','xaxis_1'))
        var_hz = dst.createVariable(output_variable_hz, varType, ('Time','zaxis_1','yaxis_1','xaxis_1'))
  
        var_vt[:] = vtScales
        var_hz[0,0:nz-1,:,:] = hz_scales[:,:]


# Main call ========================================================= CHJ =====
if __name__ == "__main__":
    main()
