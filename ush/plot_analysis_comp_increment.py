#!/usr/bin/env python3

###################################################################### CHJ #####
## Name		: plot_analysis_comp_increment.py
## Usage	: Plot comparison of JEDI increment by analysis task
## NOAA/EPIC
## History ===============================
## V000: 2024/12/10: Chan-Hoo Jeon : Preliminary version
## V001: 2025/10/23: Chan-Hoo Jeon : Add options for SOCA C-test
## V002: 2026/03/04: Chan-Hoo Jeon : Add SOCA mask option
###################################################################### CHJ #####

import os, sys
import logging
import yaml
import numpy as np
import netCDF4 as nc
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
import matplotlib as mpl
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable


# Main part (will be called at the end) ============================= CHJ =====
def main():

    global num_tiles,c_lon
    num_tiles=6
    # center of map
    c_lon=-77.0369

    yaml_file="plot_analysis_comp_increment.yaml"
    with open(yaml_file, 'r') as f:
        yaml_data=yaml.load(f, Loader=yaml.FullLoader)
    f.close()

    cartopy_ne_path = yaml_data['cartopy_ne_path']
    TYPE_ANAL_FCST = yaml_data['TYPE_ANAL_FCST']
    fn_ice_data = yaml_data['fn_ice_data']
    fn_ice_incr = yaml_data['fn_ice_incr']
    fn_ocn_data = yaml_data['fn_ocn_data']
    fn_ocn_incr = yaml_data['fn_ocn_incr']
    fn_sfc_data = yaml_data['fn_sfc_data']
    fn_sfc_incr = yaml_data['fn_sfc_incr']
    JEDI_ALGORITHM = yaml_data['JEDI_ALGORITHM']
    JEDI_TYPE_SNOW = yaml_data['JEDI_TYPE_SNOW']
    JEDI_TYPE_SOCA = yaml_data['JEDI_TYPE_SOCA']
    JEDI_TYPE_SOIL_MOISTURE = yaml_data['JEDI_TYPE_SOIL_MOISTURE']
    orog_path = yaml_data['orog_path']
    orog_fn_base = yaml_data['orog_fn_base']
    out_fn_base_prefix = yaml_data['out_fn_base_prefix']
    PDY = yaml_data['PDY']
    PY_LOG_LEVEL = yaml_data['PY_LOG_LEVEL']
    snowdepth_vn = yaml_data['snowdepth_vn']
    work_dir = yaml_data['work_dir']
    zlvl = yaml_data['zlevel_number']

    zlvlm1 = int(zlvl)-1
    grid_soca_fn = "soca_gridspec.nc"
    grid_soca_ctest_fn = "soca_gridspec.72x35x25.nc"

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

    # Set the path to Natural Earth dataset
    cartopy.config['data_dir'] = cartopy_ne_path

    list_jedi_type = []
    if JEDI_TYPE_SNOW == "YES":
        list_jedi_type.append("snow")
    if JEDI_TYPE_SOIL_MOISTURE == "YES":
        list_jedi_type.append("soil_moisture")
    if JEDI_TYPE_SOCA == "YES":
        if TYPE_ANAL_FCST == "ctest":
            list_jedi_type.append("soca_ctest")
        else:
            list_jedi_type.append("soca")
    logging.info(f''' list of JEDI types: {list_jedi_type}''')

    var_list_sfc = []
    var_list_ocn = []
    var_list_ice = []
    for jtype in list_jedi_type:
        if jtype == "snow":
            var_list_sfc.append(snowdepth_vn)
        elif jtype == "soil_moisture":
            var_list_sfc.append("smc")
        elif jtype == "soca":
            var_list_ocn += ["Salt", "Temp", "ave_ssh", "h"]
        elif jtype == "soca_ctest":
            var_list_ocn += ["Salt", "Temp", "ave_ssh"]
            if JEDI_ALGORITHM == "3dvar":
                var_list_sfc += ["sw_rad", "latent_heat", "fric_vel"]
                var_list_ice += ["hi_h", "hs_h"]
        logging.info(f''' list of vars for sfc: {var_list_sfc}''')
        logging.info(f''' list of vars for ocn: {var_list_ocn}''')
        logging.info(f''' list of vars for ice: {var_list_ice}''')

        # sfc file
        if var_list_sfc:
            # Set output file name and title base
            out_title_base=f'''UFS-DA::COMP::SFC::{jtype}::{JEDI_ALGORITHM}::{PDY}'''
            out_fn_base=f'''{out_fn_base_prefix}sfc_{jtype}_{JEDI_ALGORITHM}_{PDY}'''
            if jtype == "snow" or jtype == "soil_moisture":
                # get lon, lat from orography
                get_geo_tile(orog_path,orog_fn_base)
                for var_nm in var_list_sfc:
                    # get data before analysis
                    var_data1 = get_data_tile(work_dir,fn_sfc_data,var_nm,zlvlm1,jtype,
                                              out_title_base,out_fn_base,'before',True)
                    # get data after analysis
                    var_data2 = get_data_tile(work_dir,fn_sfc_data,var_nm,zlvlm1,jtype,
                                              out_title_base,out_fn_base,'after',True)
                    # get increment data of analysis
                    var_data_inc = get_data_tile(work_dir,fn_sfc_incr,var_nm,zlvlm1,jtype,
                                                 out_title_base,out_fn_base,'inc',True)
                    # compare data1 and data2
                    compare_data(var_data1,var_data2,var_nm,zlvlm1,jtype,
                                 out_title_base,out_fn_base,work_dir,True)
            elif jtype == "soca_ctest":
                # get lon, lat from grid file
                grid_path = os.path.join(work_dir,"data_output") 
                get_geo_grd(grid_path,grid_soca_ctest_fn,jtype)
                for var_nm in var_list_sfc:
                    # get data before analysis
                    var_data1 = get_data(work_dir,fn_sfc_data,var_nm,zlvlm1,jtype,
                                         out_title_base,out_fn_base,'before',False)
                    # get data after analysis
                    var_data2 = get_data(work_dir,fn_sfc_data,var_nm,zlvlm1,jtype,
                                        out_title_base,out_fn_base,'after',False)
                    # get increment data of analysis
                    var_data_inc = get_data(work_dir,fn_sfc_incr,var_nm,zlvlm1,jtype,
                                            out_title_base,out_fn_base,'inc',False)
                    # compare data1 and data2
                    compare_data(var_data1,var_data2,var_nm,zlvlm1,jtype,
                                 out_title_base,out_fn_base,work_dir,False)

        # ocn file
        if var_list_ocn:
            # Set output file name and title base
            out_title_base=f'''UFS-DA::COMP::OCN::{jtype}::{JEDI_ALGORITHM}::{PDY}'''
            out_fn_base=f'''{out_fn_base_prefix}ocn_{jtype}_{JEDI_ALGORITHM}_{PDY}'''
            if jtype == "soca":
                # get lon, lat from grid file
                get_geo_grd(work_dir,grid_soca_fn,jtype)
                for var_nm in var_list_ocn:
                    # get data before analysis
                    var_data1 = get_data(work_dir,fn_ocn_data,var_nm,zlvlm1,jtype,
                                         out_title_base,out_fn_base,'before',False)
                    # get data after analysis
                    var_data2 = get_data(work_dir,fn_ocn_data,var_nm,zlvlm1,jtype,
                                        out_title_base,out_fn_base,'after',False)
                    # compare data1 and data2
                    compare_data(var_data1,var_data2,var_nm,zlvlm1,jtype,
                                 out_title_base,out_fn_base,work_dir,False)
            elif jtype == "soca_ctest":
                # get lon, lat from grid file
                grid_path = os.path.join(work_dir,"data_output")
                get_geo_grd(grid_path,grid_soca_ctest_fn,jtype)
                for var_nm in var_list_ocn:
                    # get data before analysis
                    var_data1 = get_data(work_dir,fn_ocn_data,var_nm,zlvlm1,jtype,
                                         out_title_base,out_fn_base,'before',False)
                    # get data after analysis
                    var_data2 = get_data(work_dir,fn_ocn_data,var_nm,zlvlm1,jtype,
                                        out_title_base,out_fn_base,'after',False)
                    # get increment data of analysis
                    var_data_inc = get_data(work_dir,fn_ocn_incr,var_nm,zlvlm1,jtype,
                                            out_title_base,out_fn_base,'inc',False)
                    # compare data1 and data2
                    compare_data(var_data1,var_data2,var_nm,zlvlm1,jtype,
                                 out_title_base,out_fn_base,work_dir,False)

        # ice file
        if var_list_ice:
            # Set output file name and title base
            out_title_base=f'''UFS-DA::COMP::ICE::{jtype}::{JEDI_ALGORITHM}::{PDY}'''
            out_fn_base=f'''{out_fn_base_prefix}ice_{jtype}_{JEDI_ALGORITHM}_{PDY}'''
            if jtype == "soca_ctest":
                # get lon, lat from grid file
                grid_path = os.path.join(work_dir,"data_output")
                get_geo_grd(grid_path,grid_soca_ctest_fn,jtype)
                for var_nm in var_list_ice:
                    # get data before analysis
                    var_data1 = get_data(work_dir,fn_ice_data,var_nm,zlvlm1,jtype,
                                         out_title_base,out_fn_base,'before',False)
                    # get data after analysis
                    var_data2 = get_data(work_dir,fn_ice_data,var_nm,zlvlm1,jtype,
                                        out_title_base,out_fn_base,'after',False)
                    # get increment data of analysis
                    var_data_inc = get_data(work_dir,fn_ice_incr,var_nm,zlvlm1,jtype,
                                            out_title_base,out_fn_base,'inc',False)
                    # compare data1 and data2
                    compare_data(var_data1,var_data2,var_nm,zlvlm1,jtype,
                                 out_title_base,out_fn_base,work_dir,False)


# geo lon/lat from grid file ======================================== CHJ =====
def get_geo_grd(grid_path,grid_fn,jtype):
    
    global glon,glat
    logging.info(f''' ===== geo data files ==============================================''')
    fp_data = os.path.join(grid_path,grid_fn)
    try: data_raw = nc.Dataset(fp_data)
    except: raise Exception('Could NOT find the file',fp_data)
    logging.info(f''' Variables: {list(data_raw.variables)}''')   
    # Extract geo data
    if jtype == "soca" or jtype == "soca_ctest":
        glon_o = np.ma.masked_invalid(data_raw.variables['lon'])
        glat_o = np.ma.masked_invalid(data_raw.variables['lat'])
        glon = np.squeeze(glon_o,axis=0)
        glat = np.squeeze(glat_o,axis=0)

    if jtype == "soca":
        mask_o = np.ma.masked_invalid(data_raw.variables['mask2d'])
        mask2d = np.squeeze(mask_o,axis=0)
        plot_data(mask2d,"mask","orig",0,jtype,"SOCA::GRIDSPEC","soca_gridspec",grid_path,False)

    data_raw.close()
    logging.info(f''' glon = {glon.shape}''')
    logging.info(f''' glat = {glat.shape}''')


# geo lon/lat from orography ======================================== CHJ =====
def get_geo_tile(orog_path,orog_fn_base):

    global glon,glat
    logging.info(f''' ===== geo data files ==============================================''')
    cres=orog_fn_base.split('_')[0]
    glon_all=[]
    glat_all=[]
    for it in range(num_tiles):
        itp=it+1
        fn_orog = f'''{orog_fn_base}.tile{itp}.nc'''
        fp_orog = os.path.join(orog_path,fn_orog)
        try: orog = xr.open_dataset(fp_orog)
        except: raise Exception('Could NOT find the file',fp_orog)
        # Extract longitudes, and latitudes
        geolon = np.ma.masked_invalid(orog['geolon'].data)
        geolat = np.ma.masked_invalid(orog['geolat'].data)
        glon_all.append(geolon[None,:])
        glat_all.append(geolat[None,:])

    glon = np.vstack(glon_all)
    glat = np.vstack(glat_all)
    logging.info(f''' glon = {glon.shape}''')
    logging.info(f''' glat = {glat.shape}''')


# Get data from file and plot ======================================= CHJ =====
def get_data(path_data,fn_data_base,var_nm,zlvl,jtype,out_title_base,
             out_fn_base,data_opt,opt_tile):

    logging.info(f''' ===== {jtype} :: {var_nm} :: {data_opt} ===============================''')
    var_data_all=[]
    if data_opt == 'before':
        fn_data_ext=f'''_{jtype}_before_inc'''
    elif data_opt == 'after':
        fn_data_ext=f'''_{jtype}_after_inc'''
    else:
        fn_data_ext=""
    fn_data = f'''{fn_data_base}{fn_data_ext}'''
    fp_data = os.path.join(path_data,fn_data)
    try: data_raw = nc.Dataset(fp_data)
    except: raise Exception('Could NOT find the file',fp_data)
    # Extract valid variable
    var_orig = data_raw.variables[var_nm]
    var_data = np.ma.masked_invalid(var_orig)
    var_data_dim = var_data.ndim
    logging.info(f''' {var_nm} :: dimensions = {var_data_dim}''')
    if var_data_dim == 4:
        var_data3d = np.squeeze(var_data,axis=0)
        var_data2d = var_data3d[zlvl,:,:]
    else:
        var_data2d = np.squeeze(var_data,axis=0)
    data_var = var_data2d

    if data_opt == 'inc':
        plot_increment(data_var,var_nm,data_opt,zlvl,jtype,out_title_base,
                       out_fn_base,path_data,opt_tile)
    else:
        plot_data(data_var,var_nm,data_opt,zlvl,jtype,out_title_base,
                  out_fn_base,path_data,opt_tile)

    return data_var


# Get data tiles from files and plot ================================ CHJ =====
def get_data_tile(path_data,fn_data_base,var_nm,zlvl,jtype,out_title_base,
                  out_fn_base,data_opt,opt_tile):

    logging.info(f''' ===== sfc files: {var_nm} :: {data_opt} ===============================''')
    var_data_all=[]
    if data_opt == 'before':
        fn_data_ext=f'''.nc_{jtype}_before_inc'''
    elif data_opt == 'after':
        fn_data_ext=f'''.nc_{jtype}_after_inc'''
    else:
        fn_data_ext=".nc"

    for it in range(num_tiles):
        itp=it+1
        fn_data = f'''{fn_data_base}{itp}{fn_data_ext}'''
        fp_data = os.path.join(path_data,fn_data)
        try: ds = xr.open_dataset(fp_data)
        except: raise Exception('Could NOT find the file',fp_data)
        var_data = np.ma.masked_invalid(ds[var_nm].data)
        if var_nm == 'stc' or var_nm == 'smc' or var_nm == 'slc':
            var_data3d = np.squeeze(var_data,axis=0)
            var_data2d = var_data3d[zlvl,:,:]
        else:
            var_data2d = np.squeeze(var_data,axis=0)
        var_data_all.append(var_data2d[None,:])

    data_var = np.vstack(var_data_all)

    if data_opt == 'inc':
        plot_increment(data_var,var_nm,data_opt,zlvl,jtype,out_title_base,
                       out_fn_base,path_data,opt_tile)
    else:
        plot_data(data_var,var_nm,data_opt,zlvl,jtype,out_title_base,
                  out_fn_base,path_data,opt_tile)
   
    return data_var


# Compare two data set and plot ===================================== CHJ =====
def compare_data(var_data1,var_data2,var_nm,zlvl,jtype,out_title_base,
                 out_fn_base,work_dir,opt_tile):

    logging.info(f''' ===== compare files ===============================================''')
    logging.info(f''' data 1= {var_data1.shape}''')
    logging.info(f''' data 2= {var_data2.shape}''')
    diff_data = var_data2 - var_data1
    logging.info(f''' diff. data= {diff_data.shape}''')
    plot_increment(diff_data,var_nm,'diff',zlvl,jtype,out_title_base,
                   out_fn_base,work_dir,opt_tile)


# increment/difference plot ========================================== CHJ =====
def plot_increment(plt_var,plt_var_nm,plt_out_txt,zlvl,jtype,out_title_base,
                   out_fn_base,work_dir,opt_tile):

    var_max=np.nanmax(plt_var)
    var_min=np.nanmin(plt_var)
    logging.info(f''' {plt_var_nm}: {plt_out_txt} : var_max= {var_max}''')
    logging.info(f''' {plt_var_nm}: {plt_out_txt} : var_min= {var_min}''')

    if var_max == var_min:
        cs_max = max(abs(var_max),abs(var_min))+0.1
        cs_min = cs_max*-1.0
    else:
        cs_max = max(abs(var_max),abs(var_min))
        cs_min = cs_max*-1.0

    cs_cmap='seismic'
    nm_svar='\u0394'+plt_var_nm
    n_rnd=0
    cbar_extend='neither'

    if plt_var_nm == 'snodl' or plt_var_nm == 'snwdph':
        cs_max = 150
        cs_min = -150
        cbar_extend='both'

    if jtype == 'soil_moisture':
        out_title=f'''{out_title_base}::{plt_var_nm}::L{zlvl+1}::{plt_out_txt}'''
        out_fn=f'''{out_fn_base}_{plt_var_nm}_z{zlvl+1}_{plt_out_txt}'''
    else:
        out_title=f'''{out_title_base}::{plt_var_nm}::{plt_out_txt}'''
        out_fn=f'''{out_fn_base}_{plt_var_nm}_{plt_out_txt}'''

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_title(out_title, fontsize=6)
    # Call background plot
    back_plot(ax)
    if opt_tile:
        for it in range(num_tiles):
            cs=ax.pcolormesh(glon[it,:,:],glat[it,:,:],plt_var[it,:,:],cmap=cs_cmap,
                             rasterized=True,vmin=cs_min,vmax=cs_max,
                             transform=ccrs.PlateCarree())
    else:
        cs=ax.pcolormesh(glon,glat,plt_var,cmap=cs_cmap,rasterized=True,
                         vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())

    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=cbar_extend)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(nm_svar,fontsize=6)
    # Output figure
    ndpi=300
    out_file(work_dir,out_fn,ndpi)


# data plot ========================================================== CHJ =====
def plot_data(plt_var,plt_var_nm,plt_out_txt,zlvl,jedi_type,out_title_base,
              out_fn_base,work_dir,opt_tile):

    var_max=np.nanmax(plt_var)
    var_min=np.nanmin(plt_var)
    logging.info(f''' var_max= {var_max}''')
    logging.info(f''' var_min= {var_min}''')
    var_max05=var_max*0.5
    var_min05=var_min*0.5
    logging.info(f''' var_max05= {var_max05}''')
    logging.info(f''' var_min05= {var_min05}''')

    if plt_var_nm == 'snodl' or plt_var_nm == 'snwdph':
        cmap_range_opt='fixed'
    elif plt_var_nm == 'mask':
        cmap_range_opt='fixed'
    else:
        cmap_range_opt='real'
    cs_cmap='gist_ncar_r'
    if cmap_range_opt=='symmetry':
        n_rnd=0
        tmp_cmp=max(abs(var_max05),abs(var_min05))
        cs_min=round(-tmp_cmp,n_rnd)
        cs_max=round(tmp_cmp,n_rnd)
        cbar_extend='both'
    elif cmap_range_opt=='round':
        n_rnd=0
        cs_min=round(var_min05,n_rnd)
        cs_max=round(var_max05,n_rnd)
        cbar_extend='both'
    elif cmap_range_opt=='real':
        cs_min=var_min
        cs_max=var_max
        cbar_extend='neither'
    elif cmap_range_opt=='fixed':
        if plt_var_nm == "mask":
            cs_cmap=plt.get_cmap('Accent',3)
            cs_min=0
            cs_max=2
            cbar_extend='neither'
        else:
            cs_min=0.0
            cs_max=800.0
            cbar_extend='both'
    else:
        sys.exit('FATAL ERROR: wrong colormap-range flag !!!')

    logging.info(f''' cs_max= {cs_max}''')
    logging.info(f''' cs_min= {cs_min}''')

    if jedi_type == 'soil_moisture':
        out_title=f'''{out_title_base}::{plt_var_nm}::L{zlvl+1}::{plt_out_txt}'''
        out_fn=f'''{out_fn_base}_{plt_var_nm}_z{zlvl+1}_{plt_out_txt}'''
    else:
        out_title=f'''{out_title_base}::{plt_var_nm}::{plt_out_txt}'''
        out_fn=f'''{out_fn_base}_{plt_var_nm}_{plt_out_txt}'''
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_title(out_title, fontsize=6)
    # Call background plot
    back_plot(ax)
    if opt_tile:
        for it in range(num_tiles):
            cs=ax.pcolormesh(glon[it,:,:],glat[it,:,:],plt_var[it,:,:],
                             cmap=cs_cmap,rasterized=True,vmin=cs_min,
                             vmax=cs_max,transform=ccrs.PlateCarree())
    else:
        cs=ax.pcolormesh(glon,glat,plt_var,cmap=cs_cmap,rasterized=True,
                         vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())

    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=cbar_extend)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(plt_var_nm,fontsize=6)
    # Output figure
    ndpi=300
    out_file(work_dir,out_fn,ndpi)


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

