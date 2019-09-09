'''
Script that generates monthly netCDF file from folder of daily netCDF files.

Retrieves total dimension size of new netCDF files.
Iterates through daily files to extract values.
Creates list of values for each variable.
Populates new netCDF file.

'''

import netCDF4
import os
import numpy as np
import datetime
import logging
import sys

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)


yearMonth = 201908
vessel = 'LMEL' #G.O. Sars
#vessel = 'OXYH2' #Nuka Arctica
directory = '/home/mve036/Desktop/201908_LMEL' 
dim_tot = {}

logging.info(f'Creating monthly netCDF file  for {vessel}, {yearMonth}')

logging.debug('Retrieving list of daily netCDF files')
for root, dirs, files in os.walk( directory, topdown=False):
    dim_tot=0
    for name in sorted(files):
        logging.debug(f'reading filename: {name}')
        dataset = netCDF4.Dataset(os.path.join(root, name))
        dim_depth = len(dataset.dimensions['DEPTH'])
        dim_time = len(dataset.dimensions['TIME'])
        dim_lat = len(dataset.dimensions['LATITUDE'])
        dim_lon = len(dataset.dimensions['LONGITUDE'])
        dim_pos = len(dataset.dimensions['POSITION'])
        dim_tot += dim_time

nc_name = 'GL_' + str(yearMonth) + '_TS_TS_'  + vessel + '.nc'
logging.debug(f'Creating new empty monthly file: {nc_name}')
dataset_m = netCDF4.Dataset(nc_name,'w',format='NETCDF4_CLASSIC')

logging.debug('Assigning dimensions')
depth_dim = dataset_m.createDimension('DEPTH',1)
time_dim = dataset_m.createDimension('TIME',dim_tot)
lat_dim = dataset_m.createDimension('LATITUDE',dim_tot)
lon_dim = dataset_m.createDimension('LONGITUDE',dim_tot)
pos_dim = dataset_m.createDimension('POSITION',dim_tot)

logging.debug('Assigning attributes/variables')
variables = dataset.variables.keys()
for var in variables:
    fill = dataset[var]._FillValue
    variable = dataset_m.createVariable(var,dataset[var].dtype,dataset[var].dimensions,fill_value = fill)
    set_attr = {}
    for attr in dataset[var].ncattrs():
        if '_FillValue' in attr:
            continue
        attr_val = dataset[var].getncattr(attr)
        set_attr[attr] = attr_val
    variable.setncatts(set_attr)

logging.debug('populate file')
for root, dirs, files in os.walk(directory, topdown=False):
    start = {}
    end = {}
    for var in dataset.variables.keys():        
        start[var]=0
        for file in sorted(files):
            dataset_d = netCDF4.Dataset(os.path.join(root, file))
            dim_len = len(dataset_d[var].dimensions) 
            if dim_len == 1:
                array = dataset_d[var][:]
                end[var] = start[var] + len(array)
                dataset_m[var][start[var]:end[var]] = array 
                start[var] = end[var]
            elif dim_len == 2:
                array = dataset_d[var][:,:]
                end[var] = start[var] + len(array)
                dataset_m[var][start[var]:end[var]] = array 
                start[var] = end[var]

set_gattr = {}
for gattr in dataset.ncattrs():
    set_gattr[gattr] = dataset.getncattr(gattr)

start_date = (datetime.datetime(1950,1,1,0,0) + datetime.timedelta(min(dataset_m['TIME'][:]))).strftime("%Y-%m-%dT%H:%M:%SZ")
end_date = (datetime.datetime(1950,1,1,0,0) + datetime.timedelta(max(dataset_m['TIME'][:]))).strftime("%Y-%m-%dT%H:%M:%SZ")
set_gattr['geospatial_lat_min'] = min(dataset_m['LATITUDE'][:])
set_gattr['geospatial_lat_max'] = max(dataset_m['LATITUDE'][:])
set_gattr['geospatial_lon_min'] = min(dataset_m['LONGITUDE'][:])
set_gattr['geospatial_lon_max'] = max(dataset_m['LONGITUDE'][:])
set_gattr['time_coverage_start'] = start_date
set_gattr['time_coverage_end'] = end_date
set_gattr['update_interval']='void'
dataset_m.setncatts(set_gattr)
dataset_m.close()

logging.info('nc-file complete')