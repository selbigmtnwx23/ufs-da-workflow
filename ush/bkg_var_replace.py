#!/usr/bin/env python3

###################################################################### CHJ #####
## Name		  : bkg_var_replace.py
## Usage	  : Replace variables of background files with those of JEDI output
## NOAA/EPIC
## History ===============================
## V000: 2025/07/09: Chan-Hoo Jeon : Preliminary version
## V001: 2025/11/13: Chan-Hoo Jeon : Add single file option
## V002: 2026/03/04: Chan-Hoo Jeon : Add grid-mask
###################################################################### CHJ #####

import os
import sys
import logging
import yaml
import xarray as xr
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt


# Main part (will be called at the end) ============================= CHJ =====
def main():

    global bkg_data_fn_suffix,fn_data_base,fn_grid_mask,jedi_out_fn_prefix,jedi_out_fn_suffix
    global new_bkg_data_fn_suffix,work_dir

    yaml_file = "bkg_var_replace.yaml"
    with open(yaml_file, 'r') as f:
        yaml_data = yaml.load(f, Loader=yaml.FullLoader)
    f.close()

    bkg_data_fn_suffix = yaml_data['bkg_data_fn_suffix']
    fn_data_base = yaml_data['fn_data_base']
    fn_grid_mask = yaml_data['fn_grid_mask']
    jedi_out_fn_prefix = yaml_data['jedi_out_fn_prefix']
    jedi_out_fn_suffix = yaml_data['jedi_out_fn_suffix']
    JEDI_TYPE_SOCA = yaml_data['JEDI_TYPE_SOCA']
    new_bkg_data_fn_suffix = yaml_data['new_bkg_data_fn_suffix']
    num_tiles = yaml_data['num_tiles']
    PY_LOG_LEVEL = yaml_data['PY_LOG_LEVEL']
    replace_opt = yaml_data['replace_opt']
    work_dir = yaml_data['work_dir']

    # Set logging config
    log_level_str = PY_LOG_LEVEL.upper()
    try:
        log_level = getattr(logging, log_level_str)
    except AttributeError:
        log_level_str = "INFO"
        log_level = logging.INFO
        print(f''' WARNING: Invalid log level "{PY_LOG_LEVEL.upper()}", set to INFO.''')
    print(f''' Python Log Level= str: {log_level_str}, attr: {log_level}''')
    logging.basicConfig(format='%(levelname)s::%(pathname)s::L%(lineno)d::%(message)s', level=log_level)
    logging.info(f''' YAML Data: {yaml_data}''')
   
    if JEDI_TYPE_SOCA == "YES":
        var_list = ["Salt", "Temp", "ave_ssh", "h"]
    else:
        var_list = ["smc"]

    if num_tiles == 0:
        # Input/output files
        bkg_data_fn = f'''{fn_data_base}{bkg_data_fn_suffix}'''
        jedi_out_fn = f'''{jedi_out_fn_prefix}{jedi_out_fn_suffix}'''
        new_bkg_data_fn = f'''{fn_data_base}{new_bkg_data_fn_suffix}'''
        # Path to input files
        bkg_data_fp = os.path.join(work_dir, bkg_data_fn)
        jedi_out_fp = os.path.join(work_dir, jedi_out_fn)
        new_bkg_data_fp = os.path.join(work_dir, new_bkg_data_fn)
        logging.info(f''' File 1: {bkg_data_fp}''')
        logging.info(f''' File 2: {jedi_out_fp}''')
        replace_var_file(var_list,bkg_data_fp,jedi_out_fp,new_bkg_data_fp,replace_opt)
        if replace_opt != "mask":
            compare_vars_two_files(jedi_out_fp,new_bkg_data_fp, var_list)
    else:
        replace_var_tile(var_list, num_tiles)


# Replace variables for single file ============================= CHJ =====
def replace_var_file(var_list,bkg_data_fp,jedi_out_fp,new_bkg_data_fp,replace_opt):

    # Open mask file
    ds_mask = Dataset(fn_grid_mask, 'r')
    mask2d_raw = ds_mask.variables['mask2d'][:]
    mask2d = np.squeeze(mask2d_raw)

    # Open NetCDF datasets
    dsA = Dataset(bkg_data_fp, 'r')
    dsB = Dataset(jedi_out_fp, 'r')
    # Create new file
    ds_out = Dataset(new_bkg_data_fp, 'w', format='NETCDF4')

    # Copy dimensions
    for name, dim in dsA.dimensions.items():
        ds_out.createDimension(name, (len(dim) if not dim.isunlimited() else None))
    
    # Copy variables
    for name, varin in dsA.variables.items():
        outVar = ds_out.createVariable(name, varin.datatype, varin.dimensions)
        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    
        if name in var_list:
            logging.info(f''' ===== Variable: {name} ===== ''')
            # Check variable existence
            if name not in dsA.variables:
                raise ValueError(f'''Variable '{name}' not found in {bkg_data_fp}''')
            if name not in dsB.variables:
                raise ValueError(f'''Variable '{name}' not found in {jedi_out_fp}''')

            varA = dsA.variables[name]
            varB = dsB.variables[name]

            logging.info(f''' File A: {varA.shape}''')
            logging.info(f''' File B: {varB.shape}''')

            if varA.shape != varB.shape:
                raise ValueError(f'''Dimension mismatch (excluding Time): {varA.shape} vs {varB.shape}''')

            # Replace data (excluding 'Time' dimension)
            if replace_opt == "mask":
                dataA = varA[:]
                dataB = varB[:]                
                if 'Time' in varA.dimensions:
                    # assume Time is first dimension
                    for t in range(dataA.shape[0]):
                        if dataA.ndim == 3:  # (Time, y, x)
                            dataA[t, mask2d == 1] = dataB[t, mask2d == 1]
                        elif dataA.ndim == 4:  # (Time, z, y, x)
                            for k in range(dataA.shape[1]):
                                dataA[t, k, mask2d == 1] = dataB[t, k, mask2d == 1]
                else:
                    if dataA.ndim == 2:  # (y, x)
                        dataA[mask2d == 1] = dataB[mask2d == 1]
                    elif dataA.ndim == 3:  # (z, y, x)
                        for k in range(dataA.shape[0]):
                            dataA[k, mask2d == 1] = dataB[k, mask2d == 1]
                
                outVar[:] = dataA
            else:
                if 'Time' in varA.dimensions and 'Time' in varB.dimensions:
                    # Use time-mean from B
                    dataB = np.mean(varB[:], axis=0)
                    dataA = varA[:]
                    for t in range(dataA.shape[0]):
                        dataA[t, ...] = dataB
                    outVar[:] = dataA
                else:
                    outVar[:] = varB[:]
    
            logging.info(f'''Replaced variable '{name}' in {bkg_data_fp} using {jedi_out_fp} (excluding 'Time').''')
    
        else:
            # Keep original variable
            outVar[:] = varin[:]
    
    # Copy global attributes
    ds_out.setncatts({k: dsA.getncattr(k) for k in dsA.ncattrs()})
    
    # Close all files
    dsA.close()
    dsB.close()
    ds_out.close()
    logging.info(f''' Variables saved to "{new_bkg_data_fp}" successfully.''')


# Check replaced variables in two files ========================= CHJ =====
def compare_vars_two_files(fileA, fileB, vars_to_check):

    ncA = Dataset(fileA, "r")
    ncB = Dataset(fileB, "r")
    logging.info(f'''Comparing variables: {vars_to_check}''')

    for var in vars_to_check:
        if var not in ncA.variables or var not in ncB.variables:
            logging.error(f'''FATAL ERROR: {var}: NOT found in both files''')
            sys.exit(1)

        varA = ncA.variables[var]
        varB = ncB.variables[var]
        # Check dimension count
        if len(varA.dimensions) != len(varB.dimensions):
            logging.error(f'''FATAL ERROR: {var}: dimension mismatch ({len(varA.dimensions)} vs {len(varB.dimensions)})''')
            sys.exit(1)

        # Determine first dimension name and its length
        first_dim = varA.dimensions[0]
        lenA0 = ncA.dimensions[first_dim].size
        lenB0 = ncB.dimensions[first_dim].size
        if lenA0 != lenB0:
            logging.error(f'''FATAL ERROR: {var}: first dimension length mismatch ({lenA0} vs {lenB0})''')
            sys.exit(1)

        # Slice excluding the first dimension
        try:
            # Build slice: first dimension 1:end, others full range
            slices = [slice(1, None)] + [slice(None)]*(len(varA.dimensions)-1)

            dataA = varA[tuple(slices)]
            dataB = varB[tuple(slices)]
        except Exception as e:
            logging.warning(f'''{var}: error slicing first dimension → {e}''')
            continue
        # Compare shapes
        if dataA.shape != dataB.shape:
            logging.error(f'''FATAL ERROR: {var}: shape mismatch after excluding first dimension''')
            sys.exit(1)

        # Compare values
        identical = np.allclose(dataA, dataB, equal_nan=True)
        if identical:
            logging.info(f'''Good Job: {var}: values identical excluding the first dimension''')
        else:
            logging.error(f'''FATAL ERROR: {var}: values differ excluding the first dimension''')

    ncA.close()
    ncB.close()


# Replace variables for tiled files ============================= CHJ =====
def replace_var_tile(var_list, num_tiles):

    for it in range(num_tiles):
        itp = it+1
        # Input and output file name
        bkg_data_fn = f'''{fn_data_base}{itp}{bkg_data_fn_suffix}'''
        jedi_out_fn = f'''{jedi_out_fn_prefix}{fn_data_base}{itp}{jedi_out_fn_suffix}'''
        new_bkg_data_fn = f'''{fn_data_base}{itp}{new_bkg_data_fn_suffix}'''
        # Path to input files
        bkg_data_fp = os.path.join(work_dir, bkg_data_fn)
        jedi_out_fp = os.path.join(work_dir, jedi_out_fn)
        logging.info(f''' File 1: {bkg_data_fp}''')
        logging.info(f''' File 2: {jedi_out_fp}''')
        # Open the NetCDF datasets
        try:
            ds1 = xr.open_dataset(bkg_data_fp)
            ds2 = xr.open_dataset(jedi_out_fp)
        except FileNotFoundError:
            logging.error(f'''Error: One of both files not found at {work_dir}''')
        except Exception as e:
            logging.error(f'''An error occurred: {e}''')

        # Check if the values of slmsk are identical in two netcdf files
        slmsk1 = ds1['slmsk'].squeeze('Time')
        slmsk2 = ds2['slmsk'].squeeze('Time')
        slmsk1_val = slmsk1.values
        slmsk2_val = slmsk2.values
        nan_count = np.sum(np.isnan(slmsk2_val))
        logging.info(f''' slmsk2: {slmsk2_val.shape}: number of NaN elements: {nan_count}''')

        chk_slmsk = np.array_equal(slmsk1_val, slmsk2_val, equal_nan=True)
        if chk_slmsk:
            logging.info(f''' The values of the Sea-Land Mask (slmsk) are identical in two files.''')
        else:
            logging.warning(f''' WARNING: The values of the Sea-Land Mask (slmsk) are NOT identical in two files !!!''')
            slmsk_diff = slmsk1_val - slmsk2_val
            non_zero_count = np.count_nonzero(slmsk_diff)
            logging.info(f''' Number of non-identical elements: {non_zero_count}''')
            plot_comp_var_tile('slmsk', slmsk1_val, slmsk2_val, itp, 0, 'msk')

        # Check the target variables and replace them with JEDI output
        for var in var_list:
            logging.info(f''' ===== Tile #: {itp}, Variable: {var} =====''')
            # Check if the variable exists in both datasets
            if var in ds1.variables and var in ds2.variables:
                # Exclude 1st dimension (Time)
                var1_3d = ds1[var].squeeze('Time')
                var2_3d = ds2[var].squeeze('Time')
              
                num_zaxis = var1_3d.shape[0]
                for iz in range(num_zaxis):
                    izp = iz+1
                    plot_comp_var_tile(var, var1_3d[iz,:,:], var2_3d[iz,:,:], itp, izp, 'var')

            else:
                logging.error(f''' Variable "{var}" not found in one or both datasets.''')
    
            # Replace the variable values in ds1 with those from ds2 (excluding 1st dimension)
            ds1[var].values[..., :, :, :] = var2_3d

            # Plot the replaced variable
            for iz in range(num_zaxis):
                izp = iz+1
                plot_comp_var_tile(var, var2_3d[iz,:,:], ds1[var].values[0,iz,:,:], itp, izp, 'chk')

        # Save the modified dataset to a new NetCDF file
        ds1.to_netcdf(new_bkg_data_fn)

        logging.info(f''' Variable "{var}" replaced and saved to "{new_bkg_data_fn}" successfully.''')
        ds1.close()
        ds2.close()
    

# Plot var values in two files for comparison ======================== CHJ =====
def plot_comp_var_tile(var_nm, var1, var2, tile_num, lyr_num, opt):

    if opt == 'msk':
        out_fn = f'''plot_comp_bkg_{var_nm}_tile{tile_num}'''
        fig1_title = f'''BKG_DATA :: {var_nm} :: Tile {tile_num}'''
        fig2_title = f'''JEDI_Output :: {var_nm} :: Tile {tile_num}'''
    elif opt == 'chk':
        out_fn = f'''plot_comp_chk_{var_nm}_layer{lyr_num}_tile{tile_num}'''        
        fig1_title = f'''JEDI_Output :: {var_nm} :: Layer {lyr_num} :: Tile {tile_num}'''
        fig2_title = f'''Replaced SFC :: {var_nm} :: Layer {lyr_num} :: Tile {tile_num}'''
    else:
        out_fn = f'''plot_comp_bkg_{var_nm}_layer{lyr_num}_tile{tile_num}'''        
        fig1_title = f'''BKG_DATA :: {var_nm} :: Layer {lyr_num} :: Tile {tile_num}'''
        fig2_title = f'''JEDI_Output :: {var_nm} :: Layer {lyr_num} :: Tile {tile_num}'''

    var_all = np.concatenate((var1, var2))
    var_max = np.nanmax(var_all)
    var_min = np.nanmin(var_all)
    logging.info(f''' {opt}:: {var_nm}, Layer: {lyr_num}, Max: {var_max}, Min: {var_min}''')
 
    cs_map = 'plasma'
    cs_max = var_max
    cs_min = var_min
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=4

    fig, axs = plt.subplots(1, 2, figsize=(7,3))
    cs1 = axs[0].pcolormesh(var1, cmap=cs_map, rasterized=True, vmin=cs_min, vmax=cs_max)
    axs[0].set_title(fig1_title, fontsize=tlb_sz+1)
    axs[0].tick_params(direction='out',length=tick_ln,width=tick_wd,labelsize=tlb_sz) 
    
    cs2 = axs[1].pcolormesh(var2, cmap=cs_map, rasterized=True, vmin=cs_min, vmax=cs_max)
    axs[1].set_title(fig2_title, fontsize=tlb_sz+1)
    axs[1].tick_params(direction='out',length=tick_ln,width=tick_wd,labelsize=tlb_sz) 

    cbar = fig.colorbar(cs2, ax=axs, orientation='vertical', fraction=0.046, pad=0.04)
    cbar.set_label(var_nm, fontsize=tlb_sz+0.5)
    cbar.ax.tick_params(length=tick_ln,width=tick_wd,labelsize=tlb_sz)

    # Output figure
    ndpi = 300
    out_file(work_dir, out_fn, ndpi)


# Background plot ==================================================== CHJ =====
def back_plot(ax):
    # Resolution of background natural earth data ('50m' or '110m')
    back_res='50m'

    fline_wd=0.5  # line width
    falpha=0.7 # transparency

    # natural_earth
    land=cfeature.NaturalEarthFeature('physical','land',back_res,
                      edgecolor='face',facecolor=cfeature.COLORS['land'],
                      alpha=falpha)
    lakes=cfeature.NaturalEarthFeature('physical','lakes',back_res,
                      edgecolor='blue',facecolor='none',
                      linewidth=fline_wd,alpha=falpha)
    coastline=cfeature.NaturalEarthFeature('physical','coastline',
                      back_res,edgecolor='black',facecolor='none',
                      linewidth=fline_wd,alpha=falpha)
    states=cfeature.NaturalEarthFeature('cultural','admin_1_states_provinces',
                      back_res,edgecolor='green',facecolor='none',
                      linewidth=fline_wd,linestyle=':',alpha=falpha)
    borders=cfeature.NaturalEarthFeature('cultural','admin_0_countries',
                      back_res,edgecolor='red',facecolor='none',
                      linewidth=fline_wd,alpha=falpha)

#    ax.add_feature(land)
#    ax.add_feature(lakes)
#    ax.add_feature(states)
#    ax.add_feature(borders)
    ax.add_feature(coastline)


# Output file ======================================================= CHJ =====
def out_file(work_dir,out_file,ndpi):
    # Output figure
    fp_out=os.path.join(work_dir,out_file)
    plt.savefig(fp_out+'.png',dpi=ndpi,bbox_inches='tight')
    plt.close('all')


# Main call ========================================================= CHJ =====
if __name__=='__main__':
    main()
